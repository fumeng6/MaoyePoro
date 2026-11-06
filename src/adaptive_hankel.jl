# Adaptive integration is used only for the rapidly decaying difference between
# a layered/half-space field and its homogeneous full-space direct field.
function hankel_gauss_rule(order)
    e = eigen(SymTridiagonal(zeros(order),[j/sqrt(4j*j-1) for j in 1:order-1]))
    (e.values,2 .* e.vectors[1,:].^2)
end
const HANKEL_GAUSS8 = hankel_gauss_rule(8)
const HANKEL_GAUSS16 = hankel_gauss_rule(16)
const HANKEL_GAUSS32 = hankel_gauss_rule(32)

# Raw ordering, before err/ett and tilt are reconstructed by the caller.
function force_hankel_origin(ws,s,sfct,distances,dk,intervals)
    nr,nz = length(distances),length(ws.grid.lzr)
    fine,coarse = zeros(ComplexF64,10,nr,nz),zeros(ComplexF64,10,nr,nz)
    for (rule,out) in ((HANKEL_GAUSS32,fine),(HANKEL_GAUSS16,coarse))
        x,w = rule
        for panel in 1:intervals, j in eachindex(x)
            k = dk*(panel-.5+.5x[j])
            y = kernel!(ws,s,k,sfct)
            weight = .5dk*w[j]*disk_filter(k*ws.input.sradius)*k
            for ir in 1:nr
                b0,b1 = besselj0(k*distances[ir]),-besselj1(k*distances[ir])
                for iz in 1:nz
                    out[1,ir,iz] += weight*b0*y[1,iz]
                    out[2,ir,iz] += weight*b1*y[2,iz]
                    out[3,ir,iz] += weight*b0*y[4,iz]
                    out[4,ir,iz] += weight*b0*y[5,iz]
                    out[6,ir,iz] += weight*b1*y[6,iz]
                    out[7,ir,iz] += weight*b1*k*y[1,iz]
                    out[8,ir,iz] += weight*b0*y[3,iz]
                    out[9,ir,iz] += weight*b0*y[7,iz]
                    out[10,ir,iz] += weight*b1*k*y[3,iz]
                end
            end
        end
    end
    error = abs.(fine-coarse)
    for iz in 1:nz,ir in 1:nr
        ett = iszero(distances[ir]) ? error[4,ir,iz]/2 : error[2,ir,iz]/distances[ir]
        error[5,ir,iz] = ett
        error[4,ir,iz] += ett
        error[7,ir,iz] += 2error[6,ir,iz]
    end
    fine,error
end

function correction_geometry(ws,distances)
    c,m = ws.input,ws.model
    reflectors = Float64[]
    c.isurfcon > 0 && push!(reflectors,0.0)
    tops = cumsum(m.h[1:end-1])
    for n in eachindex(tops)
        any(getfield(m,name)[n] != getfield(m,name)[n+1] for name in (:la,:mu,:alf,:qa,:dm)) &&
            push!(reflectors,tops[n])
    end
    isempty(reflectors) && return (eligible=false,gap=Inf,scale=0.0)
    gap = minimum(min(abs(c.zstop-z),abs(c.zsbtm-z)) for z in reflectors)
    # A source close to an interface can leave a long oscillatory correction.
    eligible = gap > 0 && (maximum(distances)+c.sradius)/gap <= 8
    (eligible=eligible,gap=gap,scale=1/(gap+c.sradius+maximum(distances)))
end

function mapped_hankel_breaks(ws,s,scale)
    c,m = ws.input,ws.model
    ks = Float64[1/(c.zsbtm-c.zstop),1/c.sradius]
    # Bracket diffusion length scales explicitly: agreement of two Gauss rules
    # on a broad interval alone can miss a narrow feature close to k=0.
    for n in eachindex(m.h)
        beta = (1+m.alf[n]*m.qa[n]/(m.la[n]+2m.mu[n]))/m.dm[n]
        kappa = sqrt(abs(s*beta))
        kappa > 0 && append!(ks,(kappa/4,kappa,4kappa))
    end
    u = Float64[0,.25,.5,.75,1]
    append!(u,(k/(scale+k) for k in ks if isfinite(k) && k > 0))
    sort!(unique!(u))
end

function correction_integrand(ws,s,sfct,distances,k)
    y = kernel!(ws,s,k,sfct)
    direct = direct_kernel!(ws,s,k,sfct)
    y .-= direct
    out = zeros(ComplexF64,10,length(distances),length(ws.input.zr))
    disk = disk_filter(k*ws.input.sradius)*k
    @inbounds for iz in axes(out,3),ir in axes(out,2)
        r = distances[ir]
        j0,j1 = besselj0(k*r),-besselj1(k*r)
        out[1,ir,iz] = y[1,iz]*j0*disk
        out[2,ir,iz] = y[2,iz]*j1*disk
        out[3,ir,iz] = y[4,iz]*j0*disk
        out[5,ir,iz] = iszero(r) ? y[5,iz]*j0*disk/2 : out[2,ir,iz]/r
        out[4,ir,iz] = y[5,iz]*j0*disk-out[5,ir,iz]
        out[6,ir,iz] = y[6,iz]*j1*disk
        out[7,ir,iz] = k*y[1,iz]*j1*disk-2out[6,ir,iz]
        out[8,ir,iz] = y[3,iz]*j0*disk
        out[9,ir,iz] = y[7,iz]*j0*disk
        out[10,ir,iz] = k*y[3,iz]*j1*disk
    end
    out
end

function adaptive_hankel_panel(f,a,b,atol,calls,maxcalls,maxdepth,depth=0)
    calls[]+24 <= maxcalls || return zeros(ComplexF64,size(atol)),fill(Inf,size(atol)),false
    fine = zeros(ComplexF64,size(atol))
    coarse = zeros(ComplexF64,size(atol))
    center,half = (a+b)/2,(b-a)/2
    for (rule,out) in ((HANKEL_GAUSS16,fine),(HANKEL_GAUSS8,coarse))
        x,w = rule
        for j in eachindex(x)
            calls[] += 1
            out .+= (half*w[j]).*f(center+half*x[j])
        end
    end
    error = abs.(fine-coarse)
    all(error .<= atol) && return fine,error,true
    depth == maxdepth && return fine,error,false
    mid = (a+b)/2
    yl,el,cl = adaptive_hankel_panel(f,a,mid,atol./2,calls,maxcalls,maxdepth,depth+1)
    yr,er,cr = adaptive_hankel_panel(f,mid,b,atol./2,calls,maxcalls,maxdepth,depth+1)
    yl+yr,el+er,cl&&cr
end

"""Integrate the boundary/interface correction with k=scale*u/(1-u).

Return `nothing` for an unsuitable geometry or Re(s)<=0. Otherwise return the
same raw ten-field ordering as `wavenumber_integral`: the last two fields are
pressure gradients, before conversion to Darcy flux. Two mapping scales and
local 8/16-point Gauss comparisons estimate absolute errors in the supplied
field units. A failed estimate must fall back to the ordinary Hankel integral.
This route does not apply to the algebraically decaying direct field.
"""
function mapped_gauss_correction(ws,s,sfct,distances,atol)
    geometry = correction_geometry(ws,distances)
    ws.direct_workspace !== nothing && geometry.eligible && real(s) > 0 || return nothing
    calls = Ref(0)
    kmax = Ref(0.0)
    results = Array{ComplexF64,3}[]
    errors = Array{Float64,3}[]
    converged = true
    for scale in (geometry.scale,2geometry.scale)
        edges = mapped_hankel_breaks(ws,s,scale)
        value,error = zeros(ComplexF64,size(atol)),zeros(size(atol))
        f = function(u)
            k = scale*u/(1-u)
            kmax[] = max(kmax[],k)
            correction_integrand(ws,s,sfct,distances,k).*(scale/(1-u)^2)
        end
        panel_budget = atol./(4(length(edges)-1))
        for j in 2:length(edges)
            y,e,ok = adaptive_hankel_panel(f,edges[j-1],edges[j],panel_budget,calls,8192,12)
            value .+= y
            error .+= e
            converged &= ok
        end
        push!(results,value)
        push!(errors,error)
    end
    error = max.(errors[1],errors[2]).+abs.(results[1].-results[2])
    converged &= all(error .<= atol)
    results[2],(converged=converged,nk=calls[],kmax=kmax[],tail_error=error)
end
