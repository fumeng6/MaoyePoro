# Nonoscillatory source-volume quadrature in a homogeneous full space.
# Geometry is independent of s.
function volume_fluid_geometry(ws,order)
    lock(ws.volume_geometry_lock) do
        get!(ws.volume_geometry,order) do
            c = ws.input
            off_diagonal = [j/sqrt(4j*j-1) for j in 1:order-1]
            gauss = eigen(SymTridiagonal(zeros(order),off_diagonal))
            nodes = gauss.values
            weights = 2 .* gauss.vectors[1,:].^2
            v = (nodes .+ 1)./2
            rho = c.sradius .* sqrt.(v)
            source_z = (c.zstop+c.zsbtm)/2 .+ (c.zsbtm-c.zstop)/2 .* nodes
            # v=rho^2/a^2: normalized parabolic radial weights are 2(1-v)dv.
            radial_weights = (1 .- v).*weights
            axial_weights = weights./2
            nphi = 4order
            angles = (2pi/nphi).*(0:nphi-1)
            source_x = rho .* transpose(cos.(angles))
            source_y = rho .* transpose(sin.(angles))
            count = order*order*nphi
            # Cache source nodes only. Receiver-distance tables would multiply
            # memory by the number of receivers and parallel frequency workers.
            x = Matrix{Float64}(undef,count,1)
            y,z,source_weight = similar(x),similar(x),similar(x)
            point = 0
            @inbounds for jz in 1:order, jp in 1:nphi, jv in 1:order
                point += 1
                x[point] = source_x[jv,jp]
                y[point] = source_y[jv,jp]
                z[point] = source_z[jz]
                source_weight[point] = radial_weights[jv]*axial_weights[jz]/nphi
            end
            (x,y,z,source_weight)
        end
    end
end

# Keeping the source at least one characteristic source dimension away from
# every receiver avoids the near-singular point kernels in this tensor quadrature.
function separated_volume_receivers(ws,distances)
    c = ws.input
    length(ws.model.h) == 1 && c.isurfcon == 0 && c.istype in (1,2,3) || return false
    dimension = max(c.sradius,(c.zsbtm-c.zstop)/2)
    all(hypot(max(r-c.sradius,0.0),max(c.zstop-z,z-c.zsbtm,0.0)) >= dimension
        for z in c.zr for r in distances)
end

"""Homogeneous source-material direct field in the raw units of the layered split.

Use the actual model's source vector and surface projection. Divide the reference
Darcy flux by its own -chi; the caller later uses the receiver material's chi.
"""
function volume_direct_spectrum(ws,s,distances,sfct,raw_atol)
    reference = ws.direct_workspace
    reference !== nothing && separated_volume_receivers(reference,distances) || return nothing
    m,c = reference.model,ws.input
    chi = m.dm[1]*m.alf[1]/m.qa[1]
    atol = copy(raw_atol)
    atol[9:10,:,:] .*= abs(chi)
    for (iz,l) in enumerate(ws.grid.lzr)
        c.isurfcon > 0 && l == 1 || continue
        atol[6,:,iz] .= Inf
        if c.isurfcon == 1
            atol[8,:,iz] .= Inf
            atol[10,:,iz] .= Inf
        else
            atol[9,:,iz] .= Inf
        end
    end
    fields,status = volume_field_spectrum(reference,s,distances,sfct;atol)
    fields[9:10,:,:] ./= -chi
    error = status.error
    error[9:10,:,:] ./= abs(chi)
    for (iz,l) in enumerate(ws.grid.lzr)
        c.isurfcon > 0 && l == 1 || continue
        # tlt=duz/dr-2ezr. Before projection duz/dr=ezr, after ezr=0.
        fields[7,:,iz] .*= -1
        fields[6,:,iz] .= 0
        error[6,:,iz] .= 0
        if c.isurfcon == 1
            fields[8,:,iz] .= 0; fields[10,:,iz] .= 0
            error[8,:,iz] .= 0; error[10,:,iz] .= 0
        else
            fields[9,:,iz] .= 0
            error[9,:,iz] .= 0
        end
    end
    fields,(converged=status.converged,nk=0,kmax=0.0,tail_error=error,npoints=status.npoints)
end

# Radial displacement of a point fluid source, including its s=0 limit.
function point_fluid_displacement(x)
    if abs(x) < 0.1
        term,total = 0.5+0im,0.5+0im
        for j in 0:12
            term *= -x*(j+2)/((j+1)*(j+3))
            total += term
        end
        return total
    end
    (-expm1(-x)-x*exp(-x))/(x*x)
end

"""Integrate all ten homogeneous full-space fluid-source fields in source coordinates.

The point pressure is Q*exp(-gamma*R)/(4pi*chi*R). Its irrotational displacement
and strain follow from div(u)=alpha*p/(lambda+2mu). Source geometry and physical
units are identical to the Hankel representation. Receivers must be separated
from the cylinder by at least max(radius,half length); source-interior and near
source receivers retain the depth-integrated Hankel representation.
"""
function volume_field_spectrum(ws,s,distances,sfct;atol=nothing)
    separated_volume_receivers(ws,distances) && all(iszero,view(sfct,1:5)) ||
        throw(ArgumentError("Full source-volume fields require separated receivers and a homogeneous full-space fluid source."))
    c,m = ws.input,ws.model
    xi,alpha = m.la[1]+2m.mu[1],m.alf[1]
    chi = m.dm[1]*alpha/m.qa[1]
    gamma = sqrt(ComplexF64(s*(1+alpha*m.qa[1]/xi)/m.dm[1]))
    rate = -2pi*(c.zsbtm-c.zstop)*sfct[6]
    pressure_factor = rate/(4pi*chi)
    displacement_factor = alpha*pressure_factor/xi
    flux_factor = rate/(4pi)
    shape = (10,length(distances),length(c.zr))
    previous = nothing
    error = fill(Inf,shape)
    local fields,npoints
    for order in (4,8,16,32,64)
        sx,sy,sz,weight = volume_fluid_geometry(ws,order)
        fields = zeros(ComplexF64,shape)
        npoints = length(sx)*length(distances)*length(c.zr)
        @inbounds for iz in eachindex(c.zr),ir in eachindex(distances),j in eachindex(sx)
            dx,dy,dz = distances[ir]-sx[j],-sy[j],c.zr[iz]-sz[j]
            R = sqrt(dx*dx+dy*dy+dz*dz)
            nx,ny,nz = dx/R,dy/R,dz/R
            x = gamma*R
            exponential = exp(-x)
            pressure = pressure_factor*exponential/R
            displacement = displacement_factor*point_fluid_displacement(x)
            isotropic = displacement/R
            directional = alpha*pressure/xi-3isotropic
            flux = flux_factor*exponential*(1+x)/(R*R)
            w = weight[j]
            fields[1,ir,iz] += w*displacement*nz
            fields[2,ir,iz] += w*displacement*nx
            fields[3,ir,iz] += w*(isotropic+directional*nz*nz)
            fields[4,ir,iz] += w*(isotropic+directional*nx*nx)
            fields[5,ir,iz] += w*(isotropic+directional*ny*ny)
            fields[6,ir,iz] += w*directional*nx*nz
            fields[7,ir,iz] -= w*directional*nx*nz
            fields[8,ir,iz] += w*pressure
            fields[9,ir,iz] += w*flux*nz
            fields[10,ir,iz] += w*flux*nx
        end
        for iz in eachindex(c.zr),ir in eachindex(distances)
            if 2c.zr[iz] == c.zstop+c.zsbtm
                for i in (1,6,7,9); fields[i,ir,iz]=0; end
            end
            if iszero(distances[ir])
                for i in (2,6,7,10); fields[i,ir,iz]=0; end
                fields[4,ir,iz] = fields[5,ir,iz] = (fields[4,ir,iz]+fields[5,ir,iz])/2
            end
        end
        if previous !== nothing
            error .= abs.(fields.-previous)
            converged = atol === nothing ?
                all(error .<= (c.accuracy/4).*max.(abs.(fields),abs.(previous))) :
                all(error .<= atol)
            converged && return fields,(converged=true,error=error,npoints=npoints)
        end
        previous = fields
    end
    fields,(converged=false,error=error,npoints=npoints)
end

"""Compute pressure and Darcy flux by integrating the full-space diffusion Green function.

This path requires a homogeneous full space, source type 1, 2 or 3, a pure fluid-source
jump and receivers strictly outside the cylindrical source. Returned fields and
absolute quadrature errors have shape (3,radius,depth), ordered p [Pa], q_z [m/s],
q_r [m/s]. The source vector supplies the equivalent rate Q=-2pi*L*sfct[6].
For type 1 the caller multiplies these fields by s to apply the volume derivative.

Gauss orders 8,16,32,64 are compared successively. With `atol`, each physical
field uses its supplied absolute tolerance; otherwise each component must agree
to accuracy/4 relative to its own amplitude. No weak-field floor is applied.
The final Boolean records convergence; `npoints` counts source evaluations at
the final order across the requested receiver positions.
"""
function volume_fluid_spectrum(ws,s,distances,sfct;atol=nothing)
    c,m = ws.input,ws.model
    length(m.h) == 1 && c.isurfcon == 0 && c.istype in (1,2,3) &&
        all(iszero,view(sfct,1:5)) ||
        throw(ArgumentError("Source-volume fluid quadrature requires a homogeneous full space and a pure fluid source of type 1, 2 or 3."))
    all((r > c.sradius || z < c.zstop || z > c.zsbtm) for z in c.zr for r in distances) ||
        throw(ArgumentError("Source-volume fluid quadrature requires receivers strictly outside the source cylinder."))
    radial_indices = searchsortedfirst.(Ref(c.r),distances)
    all(i <= length(c.r) && c.r[i] == r for (i,r) in zip(radial_indices,distances)) ||
        throw(ArgumentError("Receiver distances must be taken from the workspace input grid."))
    nr,nz = length(distances),length(c.zr)
    shape = (3,nr,nz)
    isnothing(atol) || size(atol) == shape ||
        throw(ArgumentError("Fluid-field absolute tolerances must have shape (3,radius,depth)."))
    xi = m.la[1]+2m.mu[1]
    beta = (1+m.alf[1]*m.qa[1]/xi)/m.dm[1]
    chi = m.dm[1]*m.alf[1]/m.qa[1]
    gamma = sqrt(ComplexF64(s*beta))
    rate = -2pi*(c.zsbtm-c.zstop)*sfct[6]
    pressure_factor = rate/(4pi*chi)
    flux_factor = rate/(4pi)
    previous = nothing
    errors = fill(Inf,shape)
    local fields,npoints
    converged = false
    for order in (8,16,32,64)
        source_x,source_y,source_z,weight = volume_fluid_geometry(ws,order)
        fields = zeros(ComplexF64,shape)
        npoints = length(source_x)*nr*nz
        for iz in 1:nz, ir in 1:nr
            pressure,qz,qr = 0.0+0im,0.0+0im,0.0+0im
            @inbounds for point in eachindex(source_x)
                dx = distances[ir]-source_x[point]
                dz = c.zr[iz]-source_z[point]
                r = sqrt(dx*dx+source_y[point]^2+dz*dz)
                green = weight[point]*exp(-gamma*r)/r
                gradient = green*(gamma+1/r)
                pressure += green
                qz += gradient*dz/r
                qr += gradient*dx/r
            end
            fields[1,ir,iz] = pressure_factor*pressure
            fields[2,ir,iz] = 2c.zr[iz] == c.zstop+c.zsbtm ? 0im : flux_factor*qz
            fields[3,ir,iz] = iszero(distances[ir]) ? 0im : flux_factor*qr
        end
        if previous !== nothing
            errors .= abs.(fields.-previous)
            converged = if isnothing(atol)
                all(errors .<= (c.accuracy/4).*max.(abs.(fields),abs.(previous)))
            else
                all(errors .<= atol)
            end
            converged && return fields,errors,true,npoints
        end
        previous = fields
    end
    fields,errors,converged,npoints
end
