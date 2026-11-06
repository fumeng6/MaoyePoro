# Hankel/Bessel 径向积分；按圆柱源内的孔压体积积分确定给定平均孔压时的响应幅度。
# 记录高波数积分剩余量的估计和采样加密误差；
# 不能把某个频点的相对尾差直接当作全部物理时序的误差界。
const NFMIN = 512
const NBSJMIN = 512
const NBSJMAX = 16384
const NDMAX = 2
const DEFAULT_IMAGE_FACTOR = 1e-3
# Euler--Maclaurin half-axis correction from the first six integrand samples.
const HANKEL_ENDPOINT_WEIGHTS = (138241, -309047, 396502, -291754, 115385, -19087) ./ 60480

function frequency_capacity(input,max_frequencies=nothing)
    isnothing(max_frequencies) && return nextpow(2,max(input.nt,16NFMIN))
    max_frequencies isa Integer && max_frequencies >= NFMIN && ispow2(max_frequencies) ||
        error("max_frequencies must be a power of two, at least $NFMIN.")
    return Int(max_frequencies)
end

"""Hankel factor 8J₂(x)/x² of the normalized parabolic radial source density; x=k*radius."""
disk_filter(x) = abs(x) < 1e-3 ? 1-x*x/12+x^4/384 : 8besselj(2, x)/(x*x)

function source_vector(ws::KernelWorkspace, s; normalization_status=nothing)
    c, m, g = ws.input, ws.model, ws.grid
    sfct = zeros(ComplexF64, 6)
    len = c.zsbtm-c.zstop
    if c.istype == 0
        sfct[2] = 1/(2pi*(len == 0 ? 1 : len))
    elseif c.istype == 1
        sfct[6] = m.alf[g.nno[g.lstop]]/(2pi*len)
    elseif c.istype == 2
        sfct[6] = -1/(2pi*len)
    elseif c.istype == 3
        sfct[6] = -1/(2pi*len)
        # The normalization depends on source geometry, material and s;
        # reuse it when refining the receiver wavenumber spacing dk.
        spv, converged = get!(ws.source_pressure_integrals,ComplexF64(s)) do
            value,status = integrated_source_pressure(ws,s,sfct)
            (value,status.converged)
        end
        isnothing(normalization_status) || (normalization_status[] &= converged)
        sfct ./= spv
    end
    return sfct
end

"""Integrate the pore-pressure response over the cylindrical source volume for source vector sfct."""
function integrated_source_pressure(ws::KernelWorkspace, s, sfct)
    c, m, g = ws.input, ws.model, ws.grid
    n = g.nno[g.lstop]
    len = c.zsbtm-c.zstop
    uniform_fluid = c.isurfcon == 0 && length(m.h) == 1 && len > 0 &&
        !iszero(m.alf[n]) && all(iszero,view(sfct,1:5))
    beta = (1+m.alf[n]*m.qa[n]/(m.la[n]+2m.mu[n]))/m.dm[n]
    chi = m.dm[n]*m.alf[n]/m.qa[n]
    tail_tolerance = min(1e-3, c.accuracy/8)
    last_value = nothing
    sampling_error = Inf
    local value, tail_converged, ik, dk
    for refinement in 0:8
        dk = (2pi/c.sradius/128)/2^refinement
        # Each table belongs to its actual dk; reuse it across frequencies.
        while length(ws.source_radial_weights) <= refinement
            push!(ws.source_radial_weights, NTuple{2,Float64}[])
        end
        weights = ws.source_radial_weights[refinement+1]
        acc, previous = 0.0+0.0im, 0.0+0.0im
        ik, block, tail_converged = 0, 128, false
        while ik < NBSJMAX
            block_length = min(block, NBSJMAX-ik)
            lastk = ik+block_length
            if length(weights) < lastk
                firstk = length(weights)+1
                resize!(weights, lastk)
                for j in firstk:lastk
                    x = (j*dk)*c.sradius
                    weights[j] = (disk_filter(x), besselj1(x))
                end
            end
            for _ in 1:block_length
                ik += 1
                k = ik*dk
                disk, radial_bessel = weights[ik]
                if uniform_fluid
                    # Integrate the full-space diffusion Green function over
                    # the uniform source depth: cp=-sfct[6]*(1+(exp(-qL)-1)/(qL))/(chi*q^2).
                    q = sqrt(ComplexF64(k*k+s*beta))
                    x = q*len
                    factor = if abs(x) < 0.1
                        # [x+exp(-x)-1]/x, without subtracting nearly equal terms.
                        term = x/2
                        total = term
                        for j in 3:12
                            term *= -x/j
                            total += term
                        end
                        total
                    else
                        1+expm1(-x)/x
                    end
                    cp = -sfct[6]*factor/(chi*q*q)
                else
                    kernel!(ws, s, k, sfct)
                    if len == 0
                        cp = ws.ypsv[5, g.lstop]
                    else
                        dmai = differential_matrix!(ws.differential, ws, s, k, n)
                        top = ntuple(j -> ws.ypsv[j, g.lstop], Val(6))
                        bottom = ntuple(j -> ws.ypsv[j, g.lsbtm], Val(6))
                        ytop = (top[1], top[2]-m.alf[n]*top[5], top[3], top[4], top[5],
                            top[6]+m.alf[n]*s*top[1])
                        # Restore total traction using the layer of each endpoint,
                        # then express q_aux in the homogeneous source layer.
                        ybtm = (bottom[1], bottom[2]-m.alf[g.nno[g.lsbtm]]*bottom[5], bottom[3], bottom[4], bottom[5],
                            bottom[6]+m.alf[n]*s*bottom[1])
                        cp = sum(dmai[5,j]*((ybtm[j]-ytop[j])/len+sfct[j]) for j in 1:6)
                    end
                end
                term = cp*disk*radial_bessel
                acc += term
                ik <= length(HANKEL_ENDPOINT_WEIGHTS) && (acc += HANKEL_ENDPOINT_WEIGHTS[ik]*term)
            end
            tail_converged = abs(acc-previous) <= tail_tolerance*abs(acc)
            tail_converged && break
            previous = acc
            block *= 2
        end
        value = acc*(2pi*len*c.sradius*dk)
        if last_value !== nothing
            scale = max(abs(value), abs(last_value))
            sampling_error = iszero(scale) ? 0.0 : abs(value-last_value)/scale
            tail_converged && sampling_error <= c.accuracy/4 &&
                return value, (converged=true, nk=ik, dk=dk, sampling_error=sampling_error)
        end
        last_value = value
    end
    return value, (converged=false, nk=ik, dk=dk, sampling_error=sampling_error)
end

function bessel_sample_count(r, r0, dk)
    for ik in NBSJMIN:NBSJMAX
        if ik*dk*max(r0, first(r)) >= 1000
            return ik
        end
    end
    NBSJMAX
end

function bessel_table(r, r0, dk; nbsj=nothing, previous=nothing)
    isnothing(nbsj) && (nbsj = bessel_sample_count(r, r0, dk))
    nbsj >= 6 || error("At least six wavenumber samples are required for the endpoint correction.")
    disk = Vector{Float64}(undef, nbsj)
    bsj = Array{Float64}(undef, nbsj, 2, length(r))
    # An enlarged table has exactly the same prefix. Copy those samples rather
    # than evaluating their special functions again for every extension.
    first_new = 1
    if previous !== nothing
        previous.dk == dk && previous.radius == r0 && size(previous.bsj,3) == length(r) && previous.nbsj <= nbsj ||
            error("An extended Bessel table must use the same spacing, source radius and receiver grid size.")
        copyto!(disk, 1, previous.disk, 1, previous.nbsj)
        bsj[1:previous.nbsj,:,:] .= previous.bsj
        first_new = previous.nbsj+1
    end
    for ik in first_new:nbsj
        disk[ik] = disk_filter(ik*dk*r0)
    end
    for ir in eachindex(r)
        if iszero(r[ir])
            for ik in first_new:nbsj
                bsj[ik,1,ir] = 1.0
                bsj[ik,2,ir] = r[ir]
            end
        else
            for ik in first_new:nbsj
                x = ik*dk*r[ir]
                bsj[ik,1,ir] = besselj0(x)
                bsj[ik,2,ir] = besselj1(x)
            end
        end
    end
    return (disk=disk, bsj=bsj, dk=dk, nbsj=nbsj, radius=r0)
end

# Internal component 4 is err+ett; component 5 is filled after integration.
function wavenumber_integral(ws, s, sfct, distances, zrs2, nd, table;
        absolute_tolerance=nothing,mode::Symbol=:full,skip_fluid::Bool=false)
    c, g = ws.input, ws.grid
    mode in (:full,:direct,:correction) || throw(ArgumentError("Wavenumber integral mode must be :full, :direct, or :correction."))
    mode == :full || !isnothing(ws.direct_workspace) ||
        throw(ArgumentError("The direct-field decomposition requires a layered or half-space fluid-source model."))
    nd > 0 && any(iszero, distances) && any(iszero, zrs2) &&
        throw(ArgumentError("A receiver on the source axis at a source endpoint requires differential order zero."))
    nr, nz = length(distances), length(g.lzr)
    dk, nbsj = table.dk, table.nbsj
    k0 = max(nd*100dk, 20pi/(c.sradius+last(distances)))
    low = zeros(ComplexF64, 10, nr, nz)
    # A buried force plane has a slowly decaying tail, but its low-k field
    # varies on the much shorter scale set by source/interface distances.
    # Integrate the first six intervals directly; apply the existing endpoint
    # rule at 6dk, where its samples resolve that variation.
    origin_samples = nd == 0 && c.istype == 0 && c.zstop == c.zsbtm && c.zstop > 0 &&
        !all(iszero,ws.model.alf) ? 6 : 0
    origin_error = zeros(Float64,10,nr,nz)
    if origin_samples > 0
        low,origin_error = force_hankel_origin(ws,s,sfct,distances,dk,origin_samples)
    end
    high = zeros(ComplexF64, 10, nr, nz)
    previous = similar(high); fill!(previous, 0)
    current = similar(high)
    cy = zeros(ComplexF64, 10, nz)
    delta = Vector{ComplexF64}(undef, 10)
    tail_error = fill(Inf, 10, nr, nz)
    isnothing(absolute_tolerance) || size(absolute_tolerance) == size(high) ||
        error("The wavenumber absolute tolerance must have shape (10, radius, depth).")
    # Euler--Maclaurin half-axis endpoint correction from the first six
    # actual integrand samples. It contains f(0)/2, f'(0)/12, -f'''(0)/720
    # and f'''''(0)/30240, without assuming a homogeneous/point-source model.
    endpoint_weights = HANKEL_ENDPOINT_WEIGHTS
    ym = zeros(ComplexF64, 10, nd+1, nz)
    y0 = similar(ym); fill!(y0, 0)
    yp = similar(ym); fill!(yp, 0)
    m2 = (0,1,0,0,0,1,1,0,0,1)
    radial = (false,true,false,false,false,true,true,false,false,true)
    nnk, nk, finished = 128, nbsj, false
    for ik in 1:nbsj
        k = ik*dk
        if mode == :direct
            y = direct_kernel!(ws,s,k,sfct)
        else
            y = kernel!(ws,s,k,sfct)
            if mode == :correction
                direct = direct_kernel!(ws,s,k,sfct)
                for j in eachindex(y)
                    y[j] -= direct[j]
                end
            end
        end
        for iz in 1:nz
            cy[1,iz],cy[2,iz],cy[3,iz],cy[4,iz],cy[5,iz] =
                y[1,iz],y[2,iz],y[4,iz],y[5,iz],0
            cy[6,iz],cy[7,iz],cy[8,iz],cy[9,iz],cy[10,iz] =
                y[6,iz],k*y[1,iz],y[3,iz],y[7,iz],k*y[3,iz]
            for i in 1:10
                cy[i,iz] *= table.disk[ik]
            end
        end
        if 1 <= ik-origin_samples <= length(endpoint_weights)
            for iz in 1:nz, ir in 1:nr, i in 1:10
                i == 5 && continue
                b = radial[i] ? -table.bsj[ik,2,ir] : table.bsj[ik,1,ir]
                low[i,ir,iz] += endpoint_weights[ik-origin_samples]*cy[i,iz]*(b*k*dk)
            end
        end
        # The force-plane route has nd=0: the low/high split would cancel
        # algebraically, so sum it directly. Retain other routes' summation.
        if origin_samples == 0
            for iz in 1:nz, i in 1:10
                ym[i,1,iz] = y0[i,1,iz]
                y0[i,1,iz] = yp[i,1,iz]
            end
            if ik <= nd
                for iz in 1:nz, i in 1:10
                    yp[i,1,iz] = 0
                end
            else
                fac = exp(-((k-nd*dk)/k0)^2)
                for iz in 1:nz, i in 1:10
                    yp[i,1,iz] = cy[i,iz]*(1-fac)
                    cy[i,iz] *= fac
                end
            end
        end
        for iz in 1:nz, ir in 1:nr, i in 1:10
            i == 5 && continue
            ik <= origin_samples && continue
            b = radial[i] ? -table.bsj[ik,2,ir] : table.bsj[ik,1,ir]
            low[i,ir,iz] += cy[i,iz]*(b*k*dk)
        end
        for id in 1:min(ik-1,nd)
            kf = k-id*dk
            for iz in 1:nz, i in 1:10
                ym[i,id+1,iz] = y0[i,id+1,iz]
                y0[i,id+1,iz] = yp[i,id+1,iz]
                yp[i,id+1,iz] = y0[i,id,iz]*(zrs2[iz]+m2[i]/kf^2) -
                    (yp[i,id,iz]-ym[i,id,iz])/(2dk*kf) -
                    (yp[i,id,iz]+ym[i,id,iz]-2y0[i,id,iz])/dk^2
            end
        end
        if origin_samples == 0 && ik > nd
            kf = (ik-nd)*dk
            for iz in 1:nz, ir in 1:nr, i in 1:10
                i == 5 && continue
                b = radial[i] ? -table.bsj[ik-nd,2,ir] : table.bsj[ik-nd,1,ir]
                high[i,ir,iz] += yp[i,nd+1,iz]*(b*kf*dk)
            end
        end
        if ik == nnk
            nnk *= 2
            finished = true
            # Compare the complete integral across successive wavenumber
            # cutoffs: the low-pass part plus the scaled differential remainder.
            for iz in 1:nz, ir in 1:nr
                den = (zrs2[iz]+distances[ir]^2)^nd
                for i in 1:10
                    current[i,ir,iz] = low[i,ir,iz]+high[i,ir,iz]/den
                end
            end
            for iz in 1:nz, ir in 1:nr
                for i in 1:10
                    delta[i] = current[i,ir,iz]-previous[i,ir,iz]
                end
                delta[5] = distances[ir] == 0 ? delta[4]/2 : delta[2]/distances[ir]
                delta[4] -= delta[5]
                delta[7] -= 2delta[6]
                for i in 1:10
                    tail_error[i,ir,iz] = abs(delta[i])+origin_error[i,ir,iz]
                end
                if isnothing(absolute_tolerance)
                    for idx in (skip_fluid ? ((1,2),(3,4,6,7)) : ((1,2),(3,4,6,7),(9,10)))
                        change = sqrt(sum(abs2(current[i,ir,iz]-previous[i,ir,iz]) for i in idx))
                        scale = sqrt(sum(abs2(current[i,ir,iz]) for i in idx))
                        finished &= change <= c.accuracy*scale
                    end
                    if !skip_fluid && (g.lzr[iz] != 1 || c.isurfcon != 1)
                        finished &= abs(current[8,ir,iz]-previous[8,ir,iz]) <= c.accuracy*abs(current[8,ir,iz])
                    end
                else
                    finished &= all(tail_error[i,ir,iz] <= absolute_tolerance[i,ir,iz] for i in 1:10)
                end
            end
            if finished
                nk = ik
                break
            end
            previous .= current
        end
    end
    for iz in 1:nz, ir in 1:nr
        den = (zrs2[iz]+distances[ir]^2)^nd
        for i in 1:10
            high[i,ir,iz] = low[i,ir,iz]+high[i,ir,iz]/den
        end
        high[5,ir,iz] = distances[ir] == 0 ? high[4,ir,iz]/2 : high[2,ir,iz]/distances[ir]
        high[4,ir,iz] -= high[5,ir,iz]
        high[7,ir,iz] -= 2high[6,ir,iz]
    end
    return high, (converged=finished, nk=nk, tail_error=tail_error,
                  kmax=nk*dk, available_kmax=nbsj*dk)
end

function spectrum_scales(y)
    nr, nz = size(y,2), size(y,3)
    out = zeros(Float64, 4, nr, nz)
    for iz in 1:nz, ir in 1:nr
        out[:,ir,iz] .= (sqrt(sum(abs2(y[i,ir,iz]) for i in (1,2))),
            sqrt(sum(abs2(y[i,ir,iz]) for i in 3:7)), abs(y[8,ir,iz]),
            sqrt(sum(abs2(y[i,ir,iz]) for i in (9,10))))
    end
    out
end

function spectrum_changed(y, previous, scale, c, g; sampling=false)
    dif = spectrum_scales(y-previous)
    for iz in axes(y,3), ir in axes(y,2)
        for i in (1,2,4)
            dif[i,ir,iz] > c.accuracy*scale[i,ir,iz] && return true
        end
        if g.lzr[iz] != 1 || c.isurfcon != 1
            dif[3,ir,iz] > c.accuracy*scale[3,ir,iz] && return true
        end
        if sampling
            for i in (3,4)
                abs(y[i,ir,iz]-previous[i,ir,iz]) > c.accuracy*abs(y[i,ir,iz]) && return true
            end
        end
    end
    false
end

