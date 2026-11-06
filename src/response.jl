# 计算给定流量、荷载等时间过程引起的响应：逐次增加频率采样点数、减半波数间距并记录收敛状态。
# 输入 r×z 是计算采样网格；任意几何查询在 queries.jl 中进行，不隐式重算。
# 源半径固定为输入 sradius；接收点网格只影响数值采样，不改变物理源分布。
# 空间积分由 integration.jl 计算，本文件叠加给定源时间过程的响应。

function response_peaks(data)
    nt,nr,nz,nc = size(data)
    peaks = zeros(nc,nr,nz)
    for i in 1:nc, iz in 1:nz, ir in 1:nr, it in 1:nt
        x = data[it,ir,iz,i]
        isfinite(x) && (peaks[i,ir,iz] = max(peaks[i,ir,iz],abs(x)))
    end
    peaks
end

function response_difference(a,b)
    scale = max.(response_peaks(a),response_peaks(b))
    error = zeros(size(scale))
    for i in axes(scale,1), iz in axes(scale,3), ir in axes(scale,2), it in axes(a,1)
        x,y = a[it,ir,iz,i],b[it,ir,iz,i]
        delta = isinf(x) && x == y ? 0.0 : abs(x-y)
        error[i,ir,iz] = max(error[i,ir,iz],isnan(delta) ? Inf : delta)
    end
    ratio = zeros(size(scale))
    for i in axes(scale,1)
        floor = 128eps(Float64)*maximum(view(scale,i,:,:))
        for iz in axes(scale,3), ir in axes(scale,2)
            den = max(scale[i,ir,iz],floor)
            ratio[i,ir,iz] = iszero(den) ? (iszero(error[i,ir,iz]) ? 0.0 : Inf) : error[i,ir,iz]/den
        end
    end
    maximum(ratio),ratio,scale
end

coefficient(C,i,ir,iz) = (linear=C[1,i,ir,iz],half=C[2,i,ir,iz],
    direct=C[3,i,ir,iz],inverse_half=C[4,i,ir,iz])

function inverse_fields(values,df,alpha,c,C;mechanical_stride=1)
    nr,nz = size(values[1],2),size(values[1],3)
    data = zeros(c.nt,nr,nz,10)
    # Inversion is sequential across fields. Reuse the Fourier buffers,
    # time primitives and exponential weights for each actual frequency grid.
    strides = mechanical_stride == 1 ? (1,) : (1,mechanical_stride)
    transforms = map(strides) do stride
        nfg = length(1:stride:length(values))
        ntg = 2nfg
        dtg = 1/(ntg*(stride*df))
        tg = (0:ntg-1).*dtg
        fs = zeros(ComplexF64,ntg)
        (spectrum=Vector{ComplexF64}(undef,nfg),fs=fs,fa=similar(fs),
            step=Vector{Float64}(undef,ntg),ramp=Vector{Float64}(undef,ntg),
            factor=(stride*df).*exp.(alpha.*tg),plan=plan_bfft!(fs))
    end
    for i in 1:10, iz in 1:nz, ir in 1:nr
        stride = i <= 7 ? mechanical_stride : 1
        workspace = transforms[i <= 7 && mechanical_stride != 1 ? 2 : 1]
        for (k,j) in enumerate(1:stride:length(values))
            workspace.spectrum[k] = values[j][i,ir,iz]
        end
        data[:,ir,iz,i] .= weighted_time_response(workspace.spectrum,stride*df,alpha,c;
            coeff=coefficient(C,i,ir,iz),zero_initial=true,
            fft_plan=workspace.plan,fft_workspace=workspace)
    end
    data
end

function field_spectrum(ws,s,distances,zrs2,nd,tables;atol=nothing,extend=true)
    c,m,g = ws.input,ws.model,ws.grid
    source_ok = Ref(true)
    sfct = source_vector(ws,s;normalization_status=source_ok)
    if separated_volume_receivers(ws,distances)
        # A smooth source-volume Green integral avoids both Bessel oscillations
        # and low-k cancellation. These are already physical p and Darcy flux.
        factor = c.istype == 1 ? s : one(s)
        volume_atol = atol === nothing ? nothing :
            (iszero(factor) ? fill(Inf,size(atol)) : atol ./ abs(factor))
        y,status = volume_field_spectrum(ws,s,distances,sfct;atol=volume_atol)
        y .*= factor
        return y,(tail_error=status.error .* abs(factor),nk=0,kmax=0.0,
            converged=status.converged,source_converged=source_ok[],
            volume_samples=status.npoints,adaptive_samples=0)
    end
    conversion = fill(ComplexF64(c.istype == 1 ? s : 1),10,length(distances),length(c.zr))
    for iz in eachindex(c.zr)
        n = g.nno[g.lzr[iz]]
        chi = m.dm[n]*m.alf[n]/m.qa[n]
        conversion[9:10,:,iz] .*= -chi
    end
    raw_atol = isnothing(atol) ? nothing : [iszero(v) ? Inf : a/abs(v) for (a,v) in zip(atol,conversion)]
    isnothing(raw_atol) || (raw_atol=reshape(raw_atol,size(conversion)))
    volume_fluid = length(m.h) == 1 && c.isurfcon == 0 && c.istype in (1,2,3) &&
        all(r > c.sradius || z < c.zstop || z > c.zsbtm for z in c.zr for r in distances)
    if volume_fluid && raw_atol !== nothing
        # Source-exterior pressure can be far below the cancellation error of
        # an oscillatory Hankel sum. Integrate its decaying Green function in
        # source coordinates and reserve the Hankel integral for deformation.
        raw_atol[8:10,:,:] .= Inf
    end
    table = tables[1]
    correction = nothing
    correction_status = nothing
    adaptive_samples = 0
    volume_direct = nothing
    direct_volume_samples = 0
    if ws.direct_workspace !== nothing && raw_atol !== nothing
        vertical_gap = minimum(max(c.zstop-z,z-c.zsbtm,0.0) for z in c.zr)
        # Outside the source depth interval the direct spectrum already has
        # an exponential tail. Keep its cheap uniform rule on coarse grids.
        if iszero(vertical_gap) || 10/vertical_gap/table.dk >= 2048
            candidate = volume_direct_spectrum(ws,s,distances,sfct,raw_atol ./ 2)
            if candidate !== nothing && candidate[2].converged
                volume_direct = candidate
                direct_volume_samples = candidate[2].npoints
            end
        end
        geometry = correction_geometry(ws,distances)
        # On a coarse grid the existing correction often needs only 256
        # samples. Switch only when resolving its decay would be expensive.
        if geometry.eligible && 10/geometry.gap/table.dk >= 2048
            mapped = mapped_gauss_correction(ws,s,sfct,distances,raw_atol ./ 2)
            if mapped !== nothing && mapped[2].converged
                correction,correction_status = mapped
                adaptive_samples = correction_status.nk
            end
        end
    end
    local y,status
    while true
        if ws.direct_workspace !== nothing && raw_atol !== nothing
            # The homogeneous direct field is cheap; interface/surface
            # corrections decay rapidly with k. Give each half the error budget.
            direct,direct_status = if volume_direct === nothing
                wavenumber_integral(ws,s,sfct,distances,zrs2,nd,table;
                    absolute_tolerance=raw_atol ./ 2,mode=:direct)
            else
                volume_direct
            end
            if correction_status === nothing || !correction_status.converged
                correction,correction_status = wavenumber_integral(ws,s,sfct,distances,zrs2,nd,table;
                    absolute_tolerance=raw_atol ./ 2,mode=:correction)
            end
            y = direct+correction
            status = (converged=direct_status.converged && correction_status.converged,
                nk=max(direct_status.nk,correction_status.nk),
                tail_error=direct_status.tail_error+correction_status.tail_error,
                kmax=max(direct_status.kmax,correction_status.kmax))
        else
            y,status = wavenumber_integral(ws,s,sfct,distances,zrs2,nd,table;
                absolute_tolerance=raw_atol,skip_fluid=volume_fluid)
        end
        (status.converged || !extend || table.nbsj >= 32NBSJMAX) && break
        n = min(2table.nbsj,32NBSJMAX)
        table=bessel_table(distances,table.radius,table.dk;nbsj=n,previous=table)
        # Every enlarged table contains the same samples as its prefix.
        # Store the largest table for subsequent integrations and frequencies.
        tables[1]=table
    end
    y .*= conversion
    tail = status.tail_error .* abs.(conversion)
    volume_ok = true
    volume_samples = direct_volume_samples
    if volume_fluid
        factor = c.istype == 1 ? s : one(s)
        fluid_atol = atol === nothing ? nothing :
            (iszero(factor) ? fill(Inf,3,length(distances),length(c.zr)) : atol[8:10,:,:] ./ abs(factor))
        fluid,fluid_error,volume_ok,volume_samples = volume_fluid_spectrum(ws,s,distances,sfct;atol=fluid_atol)
        y[8:10,:,:] .= fluid .* factor
        tail[8:10,:,:] .= fluid_error .* abs(factor)
    end
    y,(tail_error=tail,nk=status.nk,kmax=status.kmax,
       converged=status.converged && volume_ok,source_converged=source_ok[],
       volume_samples=volume_samples,adaptive_samples=adaptive_samples)
end

function field_asymptotics(ws,distances,zrs2,nd,tables,pilot_scale)
    c,g = ws.input,ws.grid
    nr,nz = length(distances),length(c.zr)
    coarse,fine = zeros(4,10,nr,nz),zeros(4,10,nr,nz)
    # Injection-rate responses have only decaying high-frequency terms.
    # None contributes a direct or growing term to the time reconstruction.
    c.istype == 2 && return coarse,fine,(spatial_converged=true,source_converged=true)
    dt = c.twindow/(c.nt-1)
    base = 8pi/dt
    xs = [1.0,4.0,16.0,64.0,256.0]
    values = Array{Float64}(undef,5,10,nr,nz)
    spatial_ok,source_ok = true,true
    for j in 1:5
        y,status = field_spectrum(ws,ComplexF64(base*xs[j]),distances,zrs2,nd,tables;
            atol=pilot_scale .* (c.accuracy/256),extend=true)
        values[j,:,:,:] .= real.(y)
        spatial_ok &= status.converged
        source_ok &= status.source_converged
    end
    # Only the presence of a sqrt(s) term changes these small fit matrices.
    # Reuse the same LU/column-pivoted QR factors for all receiver components.
    fits = map((false,true)) do half
        powers = half ? [0.5,0.0,-0.5,-1.0] : [0.0,-0.5,-1.0]
        factors = map((1:4,2:5)) do indices
            matrix = hcat([xs[indices] .^ p for p in powers]...)
            length(powers) == 4 ? lu(matrix) : qr(matrix,ColumnNorm())
        end
        (powers,factors)
    end
    for iz in 1:nz, ir in 1:nr, i in 1:10
        l = g.lzr[iz]
        interface = l > 1 && g.nno[l] != g.nno[l-1]
        endpoint = (c.zr[iz] == c.zstop || c.zr[iz] == c.zsbtm) && distances[ir] < c.sradius
        half = i == 9 && c.istype != 2 &&
            ((l == 1 && c.isurfcon == 1) || interface || endpoint)
        # In a homogeneous halfspace the undrained surface pressure of a
        # force disk follows its compact traction profile. Layer reflections
        # can break this property, so restrict the exact zero accordingly.
        if length(ws.model.h) == 1 && c.istype == 0 && c.zstop == c.zsbtm == 0 &&
                l == 1 && distances[ir] >= tables[1].radius
            half = false
        end
        powers,factors = fits[half ? 2 : 1]
        for (C,indices,factor) in ((coarse,1:4,factors[1]),(fine,2:5,factors[2]))
            fit = factor \ values[indices,i,ir,iz]
            for (p,a) in zip(powers,fit)
                p == 0.5 && (C[2,i,ir,iz]=a/sqrt(base))
                p == 0.0 && (C[3,i,ir,iz]=a)
                # Decaying terms improve the fit of the leading coefficients.
                # They remain in the numerical residual: an uncertain fitted
                # s^(-1/2) coefficient would multiply a growing half integral.
            end
        end
    end
    coarse,fine,(spatial_converged=spatial_ok,source_converged=source_ok)
end

function weighted_profile(ws,distances,zrs2,nd,dk,r0,minimum_samples,capacity,image_factor;
                          verbose=true,workspaces=(ws,))
    c = ws.input
    df = 0.5/c.twindow
    alpha = -df*log(image_factor)
    # Positive-time Fourier images can dominate a diffusion signal before
    # it reaches an exterior receiver. Subtracting the t=0 value does not
    # remove their time dependence. In a homogeneous full space, extend the
    # inversion period from 2T to 8T while keeping the same Bromwich line:
    # the first image is then weighted by image_factor^4, rather than
    # image_factor. The physical source history and accuracy are unchanged.
    exterior_diffusion = length(ws.model.h) == 1 && c.isurfcon == 0 && c.istype in (1,2,3) &&
            all(r > c.sradius || z < c.zstop || z > c.zsbtm for z in c.zr for r in distances)
    mechanical_stride = exterior_diffusion ? 4 : 1
    if exterior_diffusion
        df /= 4
    end
    volume_fields = separated_volume_receivers(ws,distances)
    # Geometry is independent of frequency and read-only after construction.
    # Share its lazy cache; kernel state remains private to each worker.
    for own_ws in workspaces
        own_ws === ws && continue
        own_ws.volume_geometry=ws.volume_geometry
        own_ws.volume_geometry_lock=ws.volume_geometry_lock
        reference,pilot_reference = own_ws.direct_workspace,ws.direct_workspace
        if reference !== nothing && pilot_reference !== nothing
            reference.volume_geometry=pilot_reference.volume_geometry
            reference.volume_geometry_lock=pilot_reference.volume_geometry_lock
        end
    end
    tables = [bessel_table(distances,r0,dk;nbsj=volume_fields ? 6 : minimum_samples)]
    p0,pilot_status = field_spectrum(ws,ComplexF64(alpha),distances,zrs2,nd,tables)
    p1,_ = field_spectrum(ws,complex(alpha,2pi*df),distances,zrs2,nd,tables)
    # Propagate spatial errors through the same causal step/ramp responses as
    # the fields, using only loading that can affect the output window.
    largest_dtg = 1/(2min(NFMIN÷2,capacity÷2)*df)
    weighted_scale = max.(abs.(p0).*causal_history_bound(c,ComplexF64(alpha),largest_dtg),
                          abs.(p1).*causal_history_bound(c,complex(alpha,2pi*df),largest_dtg))
    C0,C,asym_status = field_asymptotics(ws,distances,zrs2,nd,tables,max.(abs.(p0),abs.(p1)))
    values = Array{ComplexF64,3}[]
    tail_sum = zeros(Float64,size(p0))
    source_ok = asym_status.source_converged
    raw_failed,nkmax,kmax = 0,0,0.0
    volume_samples,adaptive_samples = 0,0
    # Preserve the deformation's original frequency grid and time resolution.
    # Every fourth sample of the denser diffusion grid lies on that grid.
    n = min(mechanical_stride*(NFMIN÷2),capacity÷2)
    previous = nothing
    freq_error,asym_error = Inf,Inf
    ratios = fill(Inf,size(p0))
    local data
    while true
        first_new=length(values)+1
        statuses=Vector{typeof(pilot_status)}(undef,n-first_new+1)
        resize!(values,n)
        if length(workspaces) == 1
            for lf in first_new:n
                s = complex(alpha,2pi*(lf-1)*df)
                hs = causal_history_bound(c,s,largest_dtg)
                atol = iszero(hs) ? fill(Inf,size(p0)) : weighted_scale .* (c.accuracy/(32capacity*hs))
                values[lf],statuses[lf-first_new+1] = field_spectrum(ws,s,distances,zrs2,nd,tables;atol)
            end
        else
            # Frequencies are independent. Each task owns its kernel state
            # and any enlarged table; the starting Bessel table is read-only.
            # Keep a fixed frequency-to-workspace assignment across grids so
            # source-pressure normalization caches remain useful.
            count=length(workspaces)
            worker_tables=[[tables[1]] for _ in 1:count]
            @sync for worker in 1:count
                Threads.@spawn begin
                    own_ws=workspaces[worker]
                    own_tables=worker_tables[worker]
                    start=first_new+mod(worker-first_new,count)
                    for lf in start:count:n
                        local s = complex(alpha,2pi*(lf-1)*df)
                        local hs = causal_history_bound(c,s,largest_dtg)
                        local atol = iszero(hs) ? fill(Inf,size(p0)) : weighted_scale .* (c.accuracy/(32capacity*hs))
                        values[lf],statuses[lf-first_new+1] = field_spectrum(own_ws,s,distances,zrs2,nd,own_tables;atol)
                    end
                end
            end
            for own_tables in worker_tables
                own_tables[1].nbsj > tables[1].nbsj && (tables[1]=own_tables[1])
            end
        end
        # Accumulate error estimates in increasing frequency order.
        for lf in first_new:n
            s = complex(alpha,2pi*(lf-1)*df)
            hs = causal_history_bound(c,s,largest_dtg)
            status=statuses[lf-first_new+1]
            !iszero(hs) && (tail_sum .+= status.tail_error .* (hs*(lf == 1 ? 1 : 2)))
            source_ok &= status.source_converged
            raw_failed += !status.converged
            nkmax=max(nkmax,status.nk); kmax=max(kmax,status.kmax)
            volume_samples=max(volume_samples,status.volume_samples)
            adaptive_samples=max(adaptive_samples,status.adaptive_samples)
        end
        data=inverse_fields(values,df,alpha,c,C;mechanical_stride)
        if C == C0
            asym_error=0.0
        else
            alternative=inverse_fields(values,df,alpha,c,C0;mechanical_stride)
            asym_error,_,_=response_difference(data,alternative)
        end
        if previous !== nothing
            freq_error,ratios,_ = response_difference(data,previous)
            ((freq_error <= c.accuracy/4 && asym_error <= c.accuracy/4) || n >= capacity) && break
        end
        previous=data
        n=min(2n,capacity)
    end
    peak=response_peaks(data)
    tail_estimate=tail_sum .* df
    tail_ratio=0.0
    for i in 1:10
        floor=128eps(Float64)*maximum(view(peak,i,:,:))
        for iz in axes(peak,3), ir in axes(peak,2)
            p=max(peak[i,ir,iz],floor); e=tail_estimate[i,ir,iz]
            tail_ratio=max(tail_ratio,iszero(p) ? (iszero(e) ? 0.0 : Inf) : e/p)
        end
    end
    if verbose
        spatial = volume_fields ? "source-volume Gauss" : "dk=$dk"
        println("  $spatial; $n frequencies; time refinement error=$freq_error; weighted spatial error=$tail_ratio; asymptotic error=$asym_error")
        flush(stdout)
    end
    data,(nfrequencies=n,frequency_converged=freq_error<=c.accuracy/4 && asym_error<=c.accuracy/4,
        frequency_error=freq_error,component_frequency_error=ntuple(i->maximum(view(ratios,i,:,:)),10),
        tail_estimate=tail_estimate,wavenumber_error=tail_ratio,
        wavenumber_converged=tail_ratio<=c.accuracy/4 && source_ok && asym_status.spatial_converged,
        asymptotic_error=asym_error,asymptotic_spatial_converged=asym_status.spatial_converged,
        source_normalization_converged=source_ok,wavenumber_samples=nkmax,kmax=kmax,
        raw_tail_failures=raw_failed,table_samples=volume_fields ? 0 : maximum(t.nbsj for t in tables),
        volume_samples=volume_samples,adaptive_samples=adaptive_samples)
end

# Rapidly convergent 2F1 series, used only for the small-parameter surface-disk formulas.
function surface_disk_series(a,b,c,m)
    term,total = 1.0,1.0
    for j in 1:64
        term *= (a+j-1)*(b+j-1)*m/((c+j-1)*j)
        updated = total+term
        updated == total && return total
        total = updated
    end
    total
end

"""Surface Boussinesq displacement and radial derivative for a parabolic disk of unit total force."""
function surface_disk_vertical(r,a,mu,nu)
    factor = (1-nu)/(2pi*mu)
    iszero(r) && return factor*8/(3a),0.0
    r == a && return factor*32/(9pi*a),-factor*16/(3pi*a*a)
    m = r < a ? (r/a)^2 : (a/r)^2
    if m < 0.1
        # These forms avoid cancellation between elliptic integrals near m=0.
        if r < a
            h = surface_disk_series(0.5,-1.5,1.0,m)
            g = surface_disk_series(1.5,-0.5,2.0,m)
            return factor*8h/(3a),-factor*4r*g/a^3
        end
        h = surface_disk_series(0.5,0.5,3.0,m)
        g = surface_disk_series(1.5,1.5,4.0,m)
        return factor*h/r,-factor*(h+m*g/6)/(r*r)
    end
    K,E = ellipk(m),ellipe(m)
    if r < a
        integral = (2(2-m)*E-(1-m)*K)/3
        derivative = ((2m-1)*E+(1-m)*K)/(3m)
        return factor*16integral/(3pi*a),-factor*16r*derivative/(pi*a^3)
    end
    integral = ((3m*m-5m+2)*K+(4m-2)*E)/(3m*m)
    derivative = ((2-m)*E-2(1-m)*K)/(m*m)
    return factor*16integral/(3pi*r),-factor*16derivative/(3pi*r*r)
end

function elastic_profile(ws,distances,zrs2,nd,dk,r0,minimum_samples;verbose=true)
    c = ws.input
    surface_disk = length(ws.model.h) == 1 && c.isurfcon > 0 && c.zstop == c.zsbtm == 0
    exact_components = 1:10
    table_samples = 0
    if surface_disk && all(iszero,c.zr)
        y = zeros(ComplexF64,10,length(distances),length(c.zr))
        status = (converged=true,nk=0,kmax=0.0,tail_error=zeros(size(y)))
    else
        tables = [bessel_table(distances,r0,dk;nbsj=minimum_samples)]
        pilot,_ = field_spectrum(ws,0im,distances,zrs2,nd,tables)
        atol = abs.(pilot) .* (c.accuracy/8)
        if surface_disk
            for iz in eachindex(c.zr)
                c.zr[iz] == 0 || continue
                for ir in eachindex(distances), i in exact_components
                    atol[i,ir,iz] = Inf
                end
            end
        end
        y,status = field_spectrum(ws,0im,distances,zrs2,nd,tables;atol)
        table_samples = only(tables).nbsj
    end
    if surface_disk
        # The surface radial displacement follows the enclosed disk force.
        # Its derivatives and the prescribed normal/shear tractions give the
        # surface strains exactly. The Boussinesq disk formulas give uz and tilt.
        la,mu = only(ws.model.la),only(ws.model.mu)
        nu = la/(2(la+mu))
        factor = -(1-2nu)/(2pi*mu*c.sradius^2)
        for iz in eachindex(c.zr)
            c.zr[iz] == 0 || continue
            for (ir,r) in enumerate(distances)
                if r <= c.sradius
                    x = (r/c.sradius)^2
                    ur = factor*r*(1-x/2)
                    err,ett,ezz = factor*(1-3x/2),factor*(1-x/2),2factor*(1-x)
                else
                    ur = factor*c.sradius^2/(2r)
                    err,ett,ezz = -ur/r,ur/r,0.0
                end
                y[2,ir,iz],y[3,ir,iz],y[4,ir,iz],y[5,ir,iz] = ur,ezz,err,ett
                y[1,ir,iz],y[7,ir,iz] = surface_disk_vertical(r,c.sradius,mu,nu)
                for i in (6,8,9,10)
                    y[i,ir,iz] = 0
                end
                for i in exact_components
                    status.tail_error[i,ir,iz] = 0
                end
            end
        end
    end
    values = [source_value(c,t) for t in range(0,c.twindow;length=c.nt)]
    data = [real(y[i,ir,iz])*values[it] for it in 1:c.nt,
        ir in eachindex(distances), iz in eachindex(c.zr), i in 1:10]
    tail_ratio = 0.0
    for i in 1:10
        floor = 128eps(Float64)*maximum(abs,view(y,i,:,:))
        for iz in eachindex(c.zr), ir in eachindex(distances)
            scale = max(abs(y[i,ir,iz]),floor)
            error = status.tail_error[i,ir,iz]
            tail_ratio = max(tail_ratio,iszero(scale) ? (iszero(error) ? 0.0 : Inf) : error/scale)
        end
    end
    verbose && println("  dk=$dk; static elastic response; relative tail=$tail_ratio")
    data,(nfrequencies=0,frequency_converged=true,frequency_error=0.0,
        component_frequency_error=ntuple(_->0.0,10),
        tail_estimate=status.tail_error .* maximum(abs,values),wavenumber_error=tail_ratio,
        wavenumber_converged=status.converged && tail_ratio<=c.accuracy/4,
        asymptotic_error=0.0,asymptotic_spatial_converged=true,
        source_normalization_converged=true,wavenumber_samples=status.nk,kmax=status.kmax,
        raw_tail_failures=0,table_samples=table_samples,volume_samples=0,adaptive_samples=0)
end

"""Compute fields for the supplied source history and refine frequency count and wavenumber spacing."""
function solve_spectra(input::PoroelasticInput,model::LayerModel,grid::LayerGrid;
        verbose=true,max_frequencies=nothing,image_factor=DEFAULT_IMAGE_FACTOR,
        parallel::Bool=Threads.nthreads()>1)
    c=input
    0 < image_factor < 1 || error("image_factor must be between zero and one.")
    capacity=frequency_capacity(c,max_frequencies)
    static_elastic = c.istype == 0 && all(iszero,model.alf)
    ws=KernelWorkspace(model,grid,c)
    workspaces=[ws]
    if parallel && !static_elastic
        for _ in 2:Threads.nthreads()
            push!(workspaces,KernelWorkspace(model,grid,c))
        end
    end
    nr,nz=length(c.r),length(c.zr)
    data=zeros(c.nt,nr,nz,10)
    diagnostics=NamedTuple[]
    zrs2=[min((z-c.zstop)^2,(z-c.zsbtm)^2) for z in c.zr]
    rsdis=[hypot(sqrt(zrs2[iz]),c.r[ir]) for ir in 1:nr,iz in 1:nz]
    positive=filter(>(0),vec(rsdis)); rsmin=isempty(positive) ? 0.0 : minimum(positive)
    zpositive=filter(>(0),sqrt.(zrs2)); zmin=isempty(zpositive) ? maximum(rsdis) : minimum(zpositive)
    r0=c.sradius
    d1=c.sradius+max(rsmin,first(c.r))
    ng=last(c.r)>0 ? max(1,1+trunc(Int,log(last(c.r)/d1)/log(5))) : 1
    beta=last(c.r)>0 ? exp(log(last(c.r)/d1)/ng) : 1.0
    first_receiver=1
    while first_receiver<=nr
        last_receiver=first_receiver
        limit=max(r0,beta*c.r[first_receiver],beta*rsmin)
        for ir in first_receiver+1:nr
            c.r[ir]<=limit && (last_receiver=ir)
        end
        r=c.r[first_receiver:last_receiver]
        volume_fields=separated_volume_receivers(ws,r)
        nd=last(r)==0 ? 0 : clamp(round(Int,log10(last(r)/hypot(r0,zmin)),RoundNearestTiesAway),0,NDMAX)
        # Differential filtering divides by (r^2+zrs2)^nd; at a source-axis
        # endpoint use the direct Hankel integral for the entire receiver group.
        iszero(first(r)) && any(iszero,zrs2) && (nd=0)
        volume_fields && (nd=0)
        dk=pi/(r0+0.5(minimum(rsdis[first_receiver,:])+maximum(rsdis[last_receiver,:])))/4
        samples=bessel_sample_count(r,r0,dk)
        verbose && println("Profile $(length(diagnostics)+1): $(first(r))–$(last(r)) m; ",
            volume_fields ? "source-volume integration" : "source-weighted order $nd")
        coarse,cd=static_elastic ? elastic_profile(ws,r,zrs2,nd,dk,r0,samples;verbose) :
            weighted_profile(ws,r,zrs2,nd,dk,r0,samples,capacity,image_factor;verbose,workspaces)
        sampling_error=Inf
        local fine,fd
        if volume_fields
            # Its Gauss-order error is already included in cd.tail_estimate;
            # repeating it at a smaller dk would be the identical calculation.
            fine,fd=coarse,cd
            sampling_error=0.0
        end
        for refinement in 1:(volume_fields ? 0 : 8)
            dk/=2
            samples=min(2max(samples,cd.table_samples),32NBSJMAX)
            fine,fd=static_elastic ? elastic_profile(ws,r,zrs2,nd,dk,r0,samples;verbose) :
                weighted_profile(ws,r,zrs2,nd,dk,r0,samples,capacity,image_factor;verbose,workspaces)
            sampling_error,_,_=response_difference(fine,coarse)
            if verbose
                println("  wavenumber sampling refinement error=$sampling_error")
                flush(stdout)
            end
            if sampling_error<=c.accuracy/4 && fd.wavenumber_converged
                break
            end
            coarse,cd=fine,fd
        end
        data[:,first_receiver:last_receiver,:,:] .= fine
        push!(diagnostics,(profile=length(diagnostics)+1,first_receiver=first_receiver,last_receiver=last_receiver,
            differential_order=nd,effective_radius=r0,dk=volume_fields ? 0.0 : dk,wavenumber_samples=fd.wavenumber_samples,
            kmax=fd.kmax,table_samples=fd.table_samples,
            source_volume_samples=fd.volume_samples,adaptive_samples=fd.adaptive_samples,
            sampling_converged=sampling_error<=c.accuracy/4,sampling_error=sampling_error,
            wavenumber_converged=fd.wavenumber_converged,wavenumber_error=fd.wavenumber_error,
            source_normalization_converged=fd.source_normalization_converged,image_factor=Float64(image_factor),
            nfrequencies=fd.nfrequencies,requested_frequencies=fd.frequency_converged ? fd.nfrequencies : 2fd.nfrequencies,
            frequency_converged=fd.frequency_converged,frequency_error=fd.frequency_error,
            asymptotic_error=fd.asymptotic_error,asymptotic_spatial_converged=fd.asymptotic_spatial_converged,
            component_frequency_error=fd.component_frequency_error,raw_tail_failures=fd.raw_tail_failures,
            singular_samples=any(isinf,fine),
            spatial_integration=volume_fields ? :source_volume :
                (ws.direct_workspace === nothing ? :hankel : :split_hankel),
            method=static_elastic ? :static_elastic : :exact_history_asymptotic))
        first_receiver=last_receiver+1
    end
    data,diagnostics
end

"""Compute the corresponding source-pressure integral or fluid injection rate.

For an injection-rate input, return mean source pore pressure times source volume
(Pa*m^3 = N*m). For prescribed drained volume change or source-pressure integral,
return fluid injection rate (m^3/s); any Dirac pulse volumes are recorded separately.
The supplied piecewise-linear time history is used without resampling.
`parallel=true` evaluates independent frequencies with separate kernel workspaces.
"""
function equivalent_source(input::PoroelasticInput,model::LayerModel,grid::LayerGrid;
        verbose=true,max_frequencies=nothing,image_factor=DEFAULT_IMAGE_FACTOR,
        parallel::Bool=Threads.nthreads()>1)
    c,m,g=input,model,grid
    c.istype != 0 || error("Equivalent source applies to source types 1–3.")
    0 < image_factor < 1 || error("image_factor must be between zero and one.")
    capacity=frequency_capacity(c,max_frequencies)
    nlayer=g.nno[g.lstop]
    jumps=source_jump_amounts(c)
    if c.istype==1
        a=-m.alf[nlayer]
        out=[a*source_slope(c,t) for t in range(0,c.twindow;length=c.nt)]
        return out,(nfrequencies=0,frequency_converged=true,frequency_error=0.0,wavenumber_converged=true,
            direct_limit_converged=true,image_factor=Float64(image_factor),quantity=:injection_rate,
            impulsive_steps=!isempty(jumps),impulses=[(time=j.time,volume=a*j.amount) for j in jumps],method=:analytic_derivative)
    end
    ws=KernelWorkspace(m,g,c)
    workspaces=[ws]
    if parallel
        for _ in 2:Threads.nthreads()
            push!(workspaces,KernelWorkspace(m,g,c))
        end
    end
    sf=zeros(ComplexF64,6); sf[6]=-1/(2pi*(c.zsbtm-c.zstop))
    df=0.25/c.twindow; alpha=-df*log(image_factor)
    source_ok=Ref(true)
    function transfer(s)
        spv,status=integrated_source_pressure(ws,s,sf)
        source_ok[] &= status.converged
        c.istype==2 ? spv : inv(spv)
    end
    coeff=(linear=0.0,half=0.0,direct=0.0,inverse_half=0.0)
    coeff_coarse=coeff
    linear_error=0.0
    if c.istype==3
        base=8pi/(c.twindow/(c.nt-1))
        xs=[1.0,4.0,16.0,64.0,256.0]
        powers=[1.0,0.5,0.0,-0.5]
        high=real.([transfer(ComplexF64(base*x)) for x in xs])
        function fitted(indices)
            fit=hcat([xs[indices] .^ p for p in powers]...) \ high[indices]
            (linear=fit[1]/base,half=fit[2]/sqrt(base),direct=fit[3],inverse_half=0.0)
        end
        coeff_coarse=fitted(1:4); coeff=fitted(2:5)
        linear_error=abs(coeff.linear-coeff_coarse.linear)/max(abs(coeff.linear),abs(coeff_coarse.linear))
    end
    values=ComplexF64[]
    n=64; previous=nothing; err=Inf; asym_error=Inf
    local out
    while true
        first_new=length(values)+1
        if length(workspaces) == 1
            for lf in first_new:n
                push!(values,transfer(complex(alpha,2pi*(lf-1)*df)))
            end
        else
            resize!(values,n)
            statuses=Vector{Bool}(undef,n-first_new+1)
            count=length(workspaces)
            @sync for worker in 1:count
                Threads.@spawn begin
                    own_ws=workspaces[worker]
                    start=first_new+mod(worker-first_new,count)
                    for lf in start:count:n
                        local s=complex(alpha,2pi*(lf-1)*df)
                        local spv,status=integrated_source_pressure(own_ws,s,sf)
                        values[lf]=c.istype==2 ? spv : inv(spv)
                        statuses[lf-first_new+1]=status.converged
                    end
                end
            end
            source_ok[] &= all(statuses)
        end
        out=weighted_time_response(values,df,alpha,c;coeff,zero_initial=true)
        alternative=weighted_time_response(values,df,alpha,c;coeff=coeff_coarse,zero_initial=true)
        asym_error,_,_=response_difference(reshape(out,:,1,1,1),reshape(alternative,:,1,1,1))
        if previous!==nothing
            err,_,_=response_difference(reshape(out,:,1,1,1),reshape(previous,:,1,1,1))
            ((err<=c.accuracy/4 && asym_error<=c.accuracy/4 && linear_error<=c.accuracy/4) || n>=capacity) && break
        end
        previous=out; n=min(2n,capacity)
    end
    verbose && println("Equivalent source: $n weighted frequencies; refinement error=$err")
    out,(nfrequencies=n,frequency_converged=err<=c.accuracy/4 && asym_error<=c.accuracy/4,
        frequency_error=err,asymptotic_error=asym_error,linear_coefficient_error=linear_error,
        wavenumber_converged=source_ok[],direct_limit_converged=linear_error<=c.accuracy/4 && asym_error<=c.accuracy/4,
        image_factor=Float64(image_factor),
        quantity=c.istype==2 ? :injection_moment : :injection_rate,impulsive_steps=c.istype==3 && !isempty(jumps),
        impulses=[(time=j.time,volume=coeff.linear*j.amount) for j in jumps],singular_samples=any(isinf,out),
        method=:exact_history_asymptotic)
end
