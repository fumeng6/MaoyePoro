# 计算分层介质中位移、应变和孔压的波数-复频率响应，使用衰减模态与基向量归一化。
# 每个并行任务独占一个可复用工作区；物理场定义和单位见最终结果接口。
# 传播时第六个状态量为 -chi*dp/dz + alpha*s*uz；其中 q_z=-chi*dp/dz 为竖向 Darcy 通量。
# chi 的单位为 m²/(Pa·s)，alpha 为 Biot 系数，s 为复频率（s⁻¹）。
mutable struct KernelWorkspace
    model::LayerModel
    grid::LayerGrid
    input::PoroelasticInput
    ypsv::Matrix{ComplexF64}
    ka::Vector{ComplexF64}
    ck3::Vector{ComplexF64}
    smalls::BitVector
    max_step::Vector{Float64}
    maup::Array{ComplexF64,3}
    maiup::Array{ComplexF64,3}
    malw::Array{ComplexF64,3}
    mailw::Array{ComplexF64,3}
    yup::Array{ComplexF64,3}
    ylw::Array{ComplexF64,3}
    ysup::Array{ComplexF64,3}
    yslw::Array{ComplexF64,3}
    c0::Matrix{ComplexF64}
    c1::Matrix{ComplexF64}
    temp::Matrix{ComplexF64}
    orth::Matrix{ComplexF64}
    result::Matrix{ComplexF64}
    boundary_matrix::Matrix{ComplexF64}
    differential::Matrix{ComplexF64}
    source_rhs::Vector{ComplexF64}
    coeff6::Matrix{ComplexF64}
    rhs6::Vector{ComplexF64}
    scale6::Vector{Float64}
    pivots6::Vector{LinearAlgebra.BlasInt}
    coeff12::Matrix{ComplexF64}
    rhs12::Vector{ComplexF64}
    scale12::Vector{Float64}
    pivots12::Vector{LinearAlgebra.BlasInt}
    source_radial_weights::Vector{Vector{NTuple{2,Float64}}}
    source_pressure_integrals::Dict{ComplexF64,Tuple{ComplexF64,Bool}}
    direct_workspace::Union{Nothing,KernelWorkspace}
    volume_geometry::Dict{Int,NTuple{4,Matrix{Float64}}}
    volume_geometry_lock::ReentrantLock
end

function KernelWorkspace(model::LayerModel, grid::LayerGrid, input::PoroelasticInput)
    if any(iszero, model.alf) && !(all(iszero, model.alf) && input.istype == 0)
        error("Zero coupling is supported only as a wholly drained elastic model with source type 0; mixed zero/nonzero coupling requires a defined degenerate limit.")
    end
    nl, np = length(model.h), length(grid.z)
    max_step = zeros(nl)
    for l in eachindex(grid.hp)
        n = grid.nno[l]
        max_step[n] = max(max_step[n],grid.hp[l])
    end
    matrix_arrays = ntuple(_ -> zeros(ComplexF64, 6, 6, max(1, np-1)), 4)
    vector_arrays = ntuple(_ -> zeros(ComplexF64, 6, 3, np), 4)
    ws = KernelWorkspace(model, grid, input, zeros(ComplexF64, 6, np),
        zeros(ComplexF64, nl), zeros(ComplexF64, nl), falses(nl), max_step,
        matrix_arrays..., vector_arrays..., zeros(ComplexF64, 6, 3),
        zeros(ComplexF64, 6, 3), zeros(ComplexF64, 6, 3),
        zeros(ComplexF64, 3, 3), zeros(ComplexF64, 7, length(grid.lzr)),
        zeros(ComplexF64, 6, 6), zeros(ComplexF64, 6, 6), zeros(ComplexF64, 6),
        zeros(ComplexF64, 6, 6), zeros(ComplexF64, 6), zeros(6), zeros(LinearAlgebra.BlasInt,6),
        zeros(ComplexF64, 12, 12), zeros(ComplexF64, 12), zeros(12), zeros(LinearAlgebra.BlasInt,12),
        Vector{NTuple{2,Float64}}[],
        Dict{ComplexF64,Tuple{ComplexF64,Bool}}(), nothing,
        Dict{Int,NTuple{4,Matrix{Float64}}}(),ReentrantLock())
    if input.istype in (1,2,3) && (input.isurfcon > 0 || nl > 1)
        # The uniform full-space field in the source material contains the
        # slowly decaying high-k response. Interfaces and the surface remain
        # in the full solution and are integrated as its difference from this field.
        n = grid.nno[grid.lstop]
        reference_model = LayerModel([0.0],[model.la[n]],[model.mu[n]],
            [model.alf[n]],[model.qa[n]],[model.dm[n]])
        reference_grid = LayerGrid(copy(grid.z),copy(grid.hp),ones(Int,np),
            copy(grid.lzr),grid.lstop,grid.lsbtm)
        reference_input = PoroelasticInput(input;isurfcon=0)
        ws.direct_workspace = KernelWorkspace(reference_model,reference_grid,reference_input)
    end
    ws
end

function update_modes!(ws, s, k)
    m = ws.model
    for n in eachindex(m.h)
        beta = (1 + m.alf[n]*m.qa[n]/(m.la[n]+2m.mu[n]))/m.dm[n]
        ws.ka[n] = sqrt(ComplexF64(k*k + s*beta))
        # The near-coincident basis contains exp(±(q-k)h). Keep its phase
        # separation small over every propagation step; otherwise use the
        # original q mode with exp(-q*h). The quotient also works at s=0.
        separation = s*beta/(ws.ka[n]+k)
        ws.smalls[n] = abs(separation)/k < 1e-3 && abs(separation)*ws.max_step[n] <= 1
        ws.ck3[n] = ws.smalls[n] ? ComplexF64(k) : ws.ka[n]
    end
end

function layer_matrix!(a, ws, s, k, z, n)
    m = ws.model
    la, mu, alpha, qa, dm = m.la[n], m.mu[n], m.alf[n], m.qa[n], m.dm[n]
    q = ws.ka[n]
    x = k*z
    etu, xi = la+mu+alpha*qa, la+2mu
    xiu = xi+alpha*qa
    cb = dm*xi/xiu
    a[1,1],a[2,1],a[3,1],a[4,1],a[5,1],a[6,1] = (1, 2mu*k, 1, 2mu*k, 0, 0)
    a[1,2],a[2,2],a[3,2],a[4,2],a[5,2],a[6,2] = (1, -2mu*k, -1, 2mu*k, 0, 0)
    a[1,3],a[2,3],a[3,3],a[4,3],a[5,3],a[6,3] = (1+etu*(1-x)/mu, 2etu*(1-x)*k,
        -1-etu*x/mu, -2etu*x*k, -2qa*k, 2alpha*dm*k*k)
    a[1,4],a[2,4],a[3,4],a[4,4],a[5,4],a[6,4] = (1+etu*(1+x)/mu, -2etu*(1+x)*k,
        1-etu*x/mu, 2etu*x*k, 2qa*k, 2alpha*dm*k*k)
    a[1,5],a[2,5],a[3,5],a[4,5],a[5,5],a[6,5] = (alpha, 2alpha*mu*k*k/q, alpha*k/q,
        2alpha*mu*k, xi*s/(cb*q), -s*alpha*xiu/qa)
    a[1,6],a[2,6],a[3,6],a[4,6],a[5,6],a[6,6] = (alpha, -2alpha*mu*k*k/q, -alpha*k/q,
        2alpha*mu*k, -xi*s/(cb*q), -s*alpha*xiu/qa)
    if ws.smalls[n]
        # (q-k)/s = beta/(q+k) also at s=0; expm1 avoids cancellation.
        beta = xiu/(xi*dm)
        delta = beta/(q+k)
        arg = delta*s*z
        phase = exp(arg)
        phases = iszero(s) ? delta*z : expm1(arg)/s
        a[1,5] *= phases
        a[2,5] = -2alpha*mu*k*delta/q + a[2,5]*phases
        a[3,5] = -alpha*delta/q + a[3,5]*phases
        a[4,5] *= phases
        a[5,5] = xi/(cb*q) + a[5,5]*phases
        a[6,5] = -alpha*xiu/qa + a[6,5]*phases
        a[1,6] *= -phases/phase
        a[2,6] = 2alpha*mu*k*delta/q - a[2,6]*phases/phase
        a[3,6] = alpha*delta/q - a[3,6]*phases/phase
        a[4,6] *= -phases/phase
        a[5,6] = -xi/(cb*q) - a[5,6]*phases/phase
        a[6,6] = -alpha*xiu/qa - a[6,6]*phases/phase
    end
    for j in 1:6
        a[6,j] += alpha*s*a[1,j]
    end
    a
end

function layer_matrix(ws,s,k,z,n)
    update_modes!(ws,s,k)
    layer_matrix!(zeros(ComplexF64,6,6),ws,s,k,z,n)
end

function inverse_layer_matrix!(a, ws, s, k, n)
    m = ws.model
    la, mu, alpha, qa, dm = m.la[n], m.mu[n], m.alf[n], m.qa[n], m.dm[n]
    q = ws.ka[n]
    et1, xi = la+mu+qa*alpha, la+2mu
    xi1 = xi+qa*alpha
    cb = dm*xi/xi1
    cv1 = qa*alpha*mu*2cb*k*k
    cv2, cv3 = 2xi*xi1, 4mu*k*xi*xi1
    cv4, cv5 = cv2*s, cv3*s
    a[3,1],a[3,2],a[3,3],a[3,4],a[3,5],a[3,6] = (mu/(2xi1), 1/(4xi1*k), -mu/(2xi1), -1/(4xi1*k), 0, 0)
    a[4,1],a[4,2],a[4,3],a[4,4],a[4,5],a[4,6] = (mu/(2xi1), -1/(4xi1*k), mu/(2xi1), -1/(4xi1*k), 0, 0)
    if ws.smalls[n]
        delta = (xi1/(xi*dm))/(q+k)
        a[1,1],a[1,2],a[1,3],a[1,4],a[1,5],a[1,6] = (0, (cv1*delta/k+xi*mu)/cv3,
            (-cv1*delta/k+xi*et1)/cv2, xi*xi1/cv3,
            alpha*cb*delta/(2xi), 0)
        a[5,1],a[5,2],a[5,3],a[5,4],a[5,5],a[5,6] = (qa*mu*2cb*k*k/cv2, cb*qa*q/cv2,
            -2mu*k*cb*qa*q/cv2, -cb*qa*k/cv2,
            cb*q/(2xi), -qa/(2alpha*xi1))
    else
        a[1,1],a[1,2],a[1,3],a[1,4],a[1,5],a[1,6] = (-cv1/cv4, (-cv1+s*xi*mu)/cv5,
            (cv1+s*xi*et1)/cv4, (cv1+s*xi*xi1)/cv5,
            -alpha*cb*k/(2xi*s), qa/(2xi1*s))
        a[5,1],a[5,2],a[5,3],a[5,4],a[5,5],a[5,6] = (qa*mu*2cb*k*k/cv4, cb*qa*q/cv4,
            -2mu*k*cb*qa*q/cv4, -cb*qa*k/cv4,
            cb*q/(2xi*s), -qa/(2alpha*xi1*s))
    end
    for j in 1:6
        parity = j in (2,3,5) ? -1 : 1
        a[2,j] = parity*a[1,j]
        a[6,j] = parity*a[5,j]
    end
    for i in 1:6
        a[i,1] -= alpha*s*a[i,6]
    end
    a
end

function inverse_layer_matrix(ws,s,k,n)
    update_modes!(ws,s,k)
    inverse_layer_matrix!(zeros(ComplexF64,6,6),ws,s,k,n)
end

function differential_matrix!(a, ws, s, k, n)
    m = ws.model
    la, mu, alpha, qa, dm = m.la[n],m.mu[n],m.alf[n],m.qa[n],m.dm[n]
    lau, xi = la+alpha*qa, la+2mu
    xiu = lau+2mu
    cv = xiu*s+xi*dm*k*k
    fill!(a,0)
    a[1,2],a[1,3] = 1/(mu*k*k),-1/k
    a[2,1] = 4mu*((lau+mu)*s+(la+mu)*dm*k*k)/cv
    a[2,4],a[2,6] = -(lau*s+la*dm*k*k)/(k*cv), 2mu*qa/cv
    a[3,1],a[3,4],a[3,6] = (lau*s+la*dm*k*k)/(k*cv),(s+dm*k*k)/(k*k*cv),qa/(k*cv)
    a[4,2] = 1/k
    a[5,1],a[5,4],a[5,6] = -2mu*qa*s/cv,qa*s/(k*cv),-xi*qa/(alpha*cv)
    a[6,5] = -dm*alpha/qa
    for j in 1:6
        a[6,j] += alpha*s*a[1,j]
    end
    for i in 1:6
        a[i,1] -= alpha*s*a[i,6]
    end
    a
end

differential_matrix(ws,s,k,n) = differential_matrix!(zeros(ComplexF64,6,6),ws,s,k,n)

# Solve the source-matching system after scaling each row by its largest
# coefficient. Matrix, right-hand-side and pivot storage are reused across k.
function solve_scaled!(coef,rhs,scale,pivots)
    @inbounds for i in axes(coef,1)
        value=abs(coef[i,1])
        for j in 2:size(coef,2)
            value=max(value,abs(coef[i,j]))
        end
        scale[i]=value
    end
    @inbounds for j in axes(coef,2), i in axes(coef,1)
        coef[i,j] /= scale[i]
    end
    @inbounds for i in eachindex(rhs)
        rhs[i] /= scale[i]
    end
    @static if VERSION >= v"1.11"
        # Julia >=1.11 accepts the workspace's pivot array for LAPACK LU.
        _,_,info = LAPACK.getrf!(coef,pivots)
        info == 0 || throw(SingularException(info))
        LAPACK.getrs!('N',coef,pivots,rhs)
    else
        ldiv!(lu!(coef),rhs)
    end
end

# Invert the selected 3x3 modal block to normalize that block of the propagated
# basis to the identity matrix.
function modal_inverse!(b, c, rows)
    i,j,k = rows
    a11,a12,a13 = c[i,1],c[i,2],c[i,3]
    a21,a22,a23 = c[j,1],c[j,2],c[j,3]
    a31,a32,a33 = c[k,1],c[k,2],c[k,3]
    d = a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31)
    iszero(d) && error("Singular poroelastic propagation basis")
    b[1,1],b[1,2],b[1,3] = (a22*a33-a23*a32)/d,(a13*a32-a12*a33)/d,(a12*a23-a13*a22)/d
    b[2,1],b[2,2],b[2,3] = (a23*a31-a21*a33)/d,(a11*a33-a13*a31)/d,(a13*a21-a11*a23)/d
    b[3,1],b[3,2],b[3,3] = (a21*a32-a22*a31)/d,(a12*a31-a11*a32)/d,(a11*a22-a12*a21)/d
    b
end

# These fixed six-state products are small enough that direct sums avoid
# the overhead of a BLAS call. Inputs and output use separate workspace storage.
@inline function mul6x6_6x3!(c,a,b)
    @inbounds for j in 1:3, i in 1:6
        value = a[i,1]*b[1,j]
        for l in 2:6
            value += a[i,l]*b[l,j]
        end
        c[i,j] = value
    end
    c
end

@inline function mul6x3_3x3!(c,a,b)
    @inbounds for j in 1:3, i in 1:6
        c[i,j] = a[i,1]*b[1,j]+a[i,2]*b[2,j]+a[i,3]*b[3,j]
    end
    c
end

function propagate!(y, ws, first, last, s, k)
    first == last && return y
    g, m = ws.grid, ws.model
    step = first < last ? 1 : -1
    for l in (first+step):step:last
        interval = step == 1 ? l-1 : l
        n, h = g.nno[interval],g.hp[interval]
        wave = (ComplexF64(exp(-k*h)),ComplexF64(exp(-k*h)),exp(-ws.ck3[n]*h))
        ai = step == 1 ? ws.maiup : ws.mailw
        a = step == 1 ? ws.maup : ws.malw
        # q_aux=q_z+alpha*s*u_z is a layer-dependent coordinate. Match
        # physical q_z across interfaces, while storing each node's lower trace.
        if step == -1 && g.nno[l+1] != n
            ws.temp .= view(y,:,:,l+1)
            delta = m.alf[n]-m.alf[g.nno[l+1]]
            for j in 1:3
                ws.temp[6,j] += delta*s*ws.temp[1,j]
            end
            mul6x6_6x3!(ws.c0, view(ai,:,:,interval), ws.temp)
        else
            mul6x6_6x3!(ws.c0, view(ai,:,:,interval), view(y,:,:,l-step))
        end
        rows = step == 1 ? (1,3,5) : (2,4,6)
        modal_inverse!(ws.orth,ws.c0,rows)
        mul6x3_3x3!(ws.c1,ws.c0,ws.orth)
        for depth in first:step:(l-step)
            mul6x3_3x3!(ws.temp,view(y,:,:,depth),ws.orth)
            for j in 1:3, i in 1:6
                y[i,j,depth] = ws.temp[i,j]*wave[j]
            end
        end
        for j in 1:3, i in 1:3
            ws.c1[rows[i],j] = i == j ? 1 : 0
            opposite = step == 1 ? 2i : 2i-1
            ws.c1[opposite,j] *= wave[i]*wave[j]
        end
        mul6x6_6x3!(view(y,:,:,l),view(a,:,:,interval),ws.c1)
        if step == 1 && g.nno[l] != n
            delta = m.alf[g.nno[l]]-m.alf[n]
            for j in 1:3
                y[6,j,l] += delta*s*y[1,j,l]
            end
        end
    end
    y
end

function solve_coupled_modes!(ws,s,k,sfct)
    g,m,p = ws.grid,ws.model,ws.input
    top,bottom,np = g.lstop,g.lsbtm,length(g.z)
    for l in 1:(bottom-1)
        layer_matrix!(view(ws.maup,:,:,l),ws,s,k,g.hp[l],g.nno[l])
        inverse_layer_matrix!(view(ws.maiup,:,:,l),ws,s,k,g.nno[l])
    end
    for l in top:(np-1)
        layer_matrix!(view(ws.malw,:,:,l),ws,s,k,-g.hp[l],g.nno[l])
        inverse_layer_matrix!(view(ws.mailw,:,:,l),ws,s,k,g.nno[l])
    end
    fill!(ws.yup,0)
    fill!(ws.ylw,0)
    if p.isurfcon == 0
        a = layer_matrix!(ws.boundary_matrix,ws,s,k,0.0,g.nno[1])
        for j in 1:3, i in 1:6
            ws.yup[i,j,1] = a[i,2j-1]
        end
    elseif p.isurfcon == 1
        ws.yup[1,1,1],ws.yup[3,2,1],ws.yup[6,3,1] = 1,1,1
    elseif p.isurfcon == 2
        ws.yup[1,1,1],ws.yup[3,2,1],ws.yup[5,3,1] = 1,1,1
        ws.yup[6,1,1] = m.alf[g.nno[1]]*s
    else
        error("Surface condition must be 0, 1, or 2")
    end
    propagate!(ws.yup,ws,1,top,s,k)
    a = layer_matrix!(ws.boundary_matrix,ws,s,k,0.0,g.nno[np])
    for j in 1:3, i in 1:6
        ws.ylw[i,j,np] = a[i,2j]
    end
    propagate!(ws.ylw,ws,np,bottom,s,k)
    if top == bottom
        coef,b = ws.coeff6,ws.rhs6
        for j in 1:3, i in 1:6
            coef[i,j] = ws.yup[i,j,top]
            coef[i,j+3] = -ws.ylw[i,j,top]
        end
        copyto!(b,sfct)
        solve_scaled!(coef,b,ws.scale6,ws.pivots6)
        for l in 1:(top-1)
            mul!(view(ws.ypsv,:,l),view(ws.yup,:,:,l),view(b,1:3))
        end
        for l in (top+1):np
            mul!(view(ws.ypsv,:,l),view(ws.ylw,:,:,l),view(b,4:6))
        end
        if top == 1 && p.istype == 0 && p.isurfcon > 0
            # A force disk on the free surface prescribes boundary traction;
            # evaluate the receiver using the material-side limit.
            mul!(view(ws.ypsv,:,top),view(ws.ylw,:,:,top),view(b,4:6))
        else
            for i in 1:6
                value = (ws.yup[i,1,top]*b[1]+ws.ylw[i,1,top]*b[4])/2
                for j in 2:3
                    value += (ws.yup[i,j,top]*b[j]+ws.ylw[i,j,top]*b[j+3])/2
                end
                ws.ypsv[i,top] = value
            end
        end
    else
        copyto!(view(ws.ysup,:,:,top),view(ws.yup,:,:,top))
        propagate!(ws.ysup,ws,top,bottom,s,k)
        copyto!(view(ws.yslw,:,:,bottom),view(ws.ylw,:,:,bottom))
        propagate!(ws.yslw,ws,bottom,top,s,k)
        b = ws.source_rhs
        mul!(b,differential_matrix!(ws.differential,ws,s,k,g.nno[top]),sfct)
        coef,bl = ws.coeff12,ws.rhs12
        fill!(coef,0)
        for j in 1:3, i in 1:6
            coef[i,j] = ws.yup[i,j,top]
            coef[i,j+3] = -ws.ysup[i,j,top]
            coef[i,j+6] = -ws.yslw[i,j,top]
            coef[i+6,j+3] = -ws.ysup[i,j,bottom]
            coef[i+6,j+6] = -ws.yslw[i,j,bottom]
            coef[i+6,j+9] = ws.ylw[i,j,bottom]
        end
        for i in 1:6
            bl[i]=b[i]
            bl[i+6]=b[i]
        end
        # At a source endpoint on a material interface, both propagated bases
        # use the lower trace. Express the source-layer particular solution in
        # that same q_aux=q_z+alpha*s*u_z coordinate at the lower endpoint.
        bl[12] += (m.alf[g.nno[bottom]]-m.alf[g.nno[top]])*s*b[1]
        solve_scaled!(coef,bl,ws.scale12,ws.pivots12)
        for l in 1:top
            mul!(view(ws.ypsv,:,l),view(ws.yup,:,:,l),view(bl,1:3))
        end
        for l in bottom:np
            mul!(view(ws.ypsv,:,l),view(ws.ylw,:,:,l),view(bl,10:12))
        end
        for l in (top+1):(bottom-1), i in 1:6
            value = bl[4]*ws.ysup[i,1,l]+bl[7]*ws.yslw[i,1,l]
            for j in 2:3
                value += bl[j+3]*ws.ysup[i,j,l]+bl[j+6]*ws.yslw[i,j,l]
            end
            ws.ypsv[i,l] = b[i]+value
        end
    end
    for l in 1:np
        alpha = m.alf[g.nno[l]]
        ws.ypsv[2,l] += alpha*ws.ypsv[5,l]
        ws.ypsv[6,l] -= alpha*s*ws.ypsv[1,l]
    end
    ws.ypsv
end

# At νu=ν the six-state poroelastic basis is singular (χ=α=0). Mechanical
# loading in this limit is computed with four elastic states and zero fluid response.
function elastic_matrix!(a,ws,k,z,n)
    m = ws.model
    mu, et = m.mu[n],m.la[n]+m.mu[n]
    x=k*z
    a[1,1],a[2,1],a[3,1],a[4,1] = 1,2mu*k,1,2mu*k
    a[1,2],a[2,2],a[3,2],a[4,2] = 1,-2mu*k,-1,2mu*k
    a[1,3],a[2,3],a[3,3],a[4,3] = 1+et*(1-x)/mu,2et*(1-x)*k,-1-et*x/mu,-2et*x*k
    a[1,4],a[2,4],a[3,4],a[4,4] = 1+et*(1+x)/mu,-2et*(1+x)*k,1-et*x/mu,2et*x*k
    a
end

elastic_matrix(ws,k,z,n) = elastic_matrix!(zeros(ComplexF64,4,4),ws,k,z,n)

function elastic_propagate!(y,ws,first,last,k)
    first == last && return y
    g=ws.grid
    step=first<last ? 1 : -1
    a=view(ws.boundary_matrix,1:4,1:4)
    c0,c1,temp=view(ws.c0,1:4,1:2),view(ws.c1,1:4,1:2),view(ws.temp,1:4,1:2)
    for l in (first+step):step:last
        interval=step==1 ? l-1 : l
        n,h=g.nno[interval],g.hp[interval]
        ai=LinearAlgebra.inv!(lu!(elastic_matrix!(a,ws,k,0.0,n)))
        mul!(c0,ai,view(y,:,:,l-step))
        rows=step==1 ? (1:2:3) : (2:2:4)
        orth=inv(c0[rows,:])
        mul!(c1,c0,orth)
        wave=exp(-k*h)
        for depth in first:step:l-step
            mul!(temp,view(y,:,:,depth),orth)
            for j in 1:2, i in 1:4
                y[i,j,depth] = temp[i,j]*wave
            end
        end
        c1[rows[1],1],c1[rows[1],2],c1[rows[2],1],c1[rows[2],2] = 1,0,0,1
        opposite=step==1 ? (2,4) : (1,3)
        for j in 1:2, i in opposite
            c1[i,j] *= wave^2
        end
        mul!(view(y,:,:,l),elastic_matrix!(a,ws,k,step*h,n),c1)
    end
    y
end

function solve_elastic_modes!(ws,k,sfct)
    g=ws.grid
    top,bottom,np=g.lstop,g.lsbtm,length(g.z)
    # The drained limit uses four states and two modes, stored in the leading
    # blocks of the coupled workspace.
    up=view(ws.yup,1:4,1:2,:)
    low=view(ws.ylw,1:4,1:2,:)
    fill!(up,0)
    fill!(low,0)
    a=view(ws.boundary_matrix,1:4,1:4)
    if ws.input.isurfcon==0
        elastic_matrix!(a,ws,k,0.0,g.nno[1])
        for j in 1:2, i in 1:4
            up[i,j,1] = a[i,2j-1]
        end
    else
        up[1,1,1],up[3,2,1]=1,1
    end
    elastic_propagate!(up,ws,1,top,k)
    elastic_matrix!(a,ws,k,0.0,g.nno[np])
    for j in 1:2, i in 1:4
        low[i,j,np] = a[i,2j]
    end
    elastic_propagate!(low,ws,np,bottom,k)
    fill!(ws.ypsv,0)
    if top==bottom
        coef,b=view(ws.coeff6,1:4,1:4),view(ws.rhs6,1:4)
        for j in 1:2, i in 1:4
            coef[i,j] = up[i,j,top]
            coef[i,j+2] = -low[i,j,top]
        end
        copyto!(b,view(sfct,1:4))
        solve_scaled!(coef,b,view(ws.scale6,1:4),view(ws.pivots6,1:4))
        for l in 1:top-1
            mul!(view(ws.ypsv,1:4,l),view(up,:,:,l),view(b,1:2))
        end
        for l in top+1:np
            mul!(view(ws.ypsv,1:4,l),view(low,:,:,l),view(b,3:4))
        end
        if top==1 && ws.input.isurfcon>0
            mul!(view(ws.ypsv,1:4,top),view(low,:,:,top),view(b,3:4))
        else
            upper,lower=view(ws.temp,1:4,1),view(ws.temp,1:4,2)
            mul!(upper,view(up,:,:,top),view(b,1:2))
            mul!(lower,view(low,:,:,top),view(b,3:4))
            for i in 1:4
                ws.ypsv[i,top] = (upper[i]+lower[i])/2
            end
        end
    else
        sup,slw=view(ws.ysup,1:4,1:2,:),view(ws.yslw,1:4,1:2,:)
        copyto!(view(sup,:,:,top),view(up,:,:,top))
        elastic_propagate!(sup,ws,top,bottom,k)
        copyto!(view(slw,:,:,bottom),view(low,:,:,bottom))
        elastic_propagate!(slw,ws,bottom,top,k)
        # A^{-1} times the constant source gives the in-cylinder particular
        # solution. A is the first-order elastic depth-system matrix.
        n=g.nno[top]
        la,mu=ws.model.la[n],ws.model.mu[n]
        xi=la+2mu
        di=view(ws.differential,1:4,1:4)
        fill!(di,0)
        di[1,2],di[1,3] = 1/(mu*k*k),-1/k
        di[2,1],di[2,4] = 4mu*(la+mu)/xi,-la/(k*xi)
        di[3,1],di[3,4],di[4,2] = la/(k*xi),1/(k*k*xi),1/k
        b=view(ws.source_rhs,1:4)
        mul!(b,di,view(sfct,1:4))
        coef,bl=view(ws.coeff12,1:8,1:8),view(ws.rhs12,1:8)
        fill!(coef,0)
        for j in 1:2, i in 1:4
            coef[i,j] = up[i,j,top]
            coef[i,j+2] = -sup[i,j,top]
            coef[i,j+4] = -slw[i,j,top]
            coef[i+4,j+2] = -sup[i,j,bottom]
            coef[i+4,j+4] = -slw[i,j,bottom]
            coef[i+4,j+6] = low[i,j,bottom]
        end
        for i in 1:4
            bl[i],bl[i+4] = b[i],b[i]
        end
        solve_scaled!(coef,bl,view(ws.scale12,1:8),view(ws.pivots12,1:8))
        for l in 1:top
            mul!(view(ws.ypsv,1:4,l),view(up,:,:,l),view(bl,1:2))
        end
        for l in bottom:np
            mul!(view(ws.ypsv,1:4,l),view(low,:,:,l),view(bl,7:8))
        end
        for l in top+1:bottom-1
            upper,lower=view(ws.temp,1:4,1),view(ws.temp,1:4,2)
            mul!(upper,view(sup,:,:,l),view(bl,3:4))
            mul!(lower,view(slw,:,:,l),view(bl,5:6))
            for i in 1:4
                ws.ypsv[i,l] = b[i]+upper[i]+lower[i]
            end
        end
    end
    ws.ypsv
end

# exp(-t*a)-exp(-t*b), retaining small depth differences.
decay_difference(t,a,b) = a <= b ? -exp(-t*a)*expm1(-t*(b-a)) :
    exp(-t*b)*expm1(-t*(a-b))

# [exp(-q*r)-exp(-k*r)]/(q-k), including its coincident-root limit.
decay_divided_difference(k,separation,r) = iszero(separation) ? -r*exp(-k*r) :
    exp(-k*r)*expm1(-separation*r)/separation

"""Integrate the homogeneous full-space fluid-source Green function over source depth.

For H_t(z)=integral(exp(-t*abs(z-zs))/(2t), zs=top:bottom), pore pressure
is J*H_q and displacement is (alpha/xi)*grad(Phi), where
Phi=J*(H_q-H_k)/(q^2-k^2), J=-sfct[6]/chi and xi=lambda+2mu.
Stable divided differences retain the zero-frequency limit without subtracting
nearly equal roots. Stored states use matrix stress and physical Darcy flux,
as after solve_coupled_modes!; receiver fields are written directly to retain
small strains that would otherwise cancel in stress-to-strain conversion.
"""
function solve_fullspace_fluid!(ws,s,k,sfct)
    c,m,g = ws.input,ws.model,ws.grid
    la,mu,alpha,qa,dm = m.la[1],m.mu[1],m.alf[1],m.qa[1],m.dm[1]
    xi = la+2mu
    chi = dm*alpha/qa
    q = ws.ka[1]
    separation = s*((1+alpha*qa/xi)/dm)/(q+k)
    amplitude = -sfct[6]/chi
    elastic = alpha*amplitude/(2xi)
    receiver = 1
    for l in eachindex(g.z)
        z = g.z[l]
        if z < c.zstop || z > c.zsbtm
            above = z < c.zstop
            a,b = above ? (c.zstop-z,c.zsbtm-z) : (z-c.zsbtm,z-c.zstop)
            direction = above ? 1 : -1
            nk,nq = decay_difference(k,a,b),decay_difference(q,a,b)
            mk,mq = direction*nk,direction*nq
            da = decay_divided_difference(k,separation,a)
            db = decay_divided_difference(k,separation,b)
            dn = da-db
            dmz = direction*dn
        else
            a,b = z-c.zstop,c.zsbtm-z
            nk = -expm1(-k*a)-expm1(-k*b)
            nq = -expm1(-q*a)-expm1(-q*b)
            mk,mq = decay_difference(k,a,b),decay_difference(q,a,b)
            da = decay_divided_difference(k,separation,a)
            db = decay_divided_difference(k,separation,b)
            dn,dmz = -da-db,da-db
        end
        divided_n = dn/(q+k)
        uz = elastic/(q*(q+k))*(dmz-mk/k)
        ur = elastic/(q*q)*(k*divided_n-nk/k)
        pressure = amplitude*nq/(2q*q)
        ezz = elastic*divided_n
        radial_strain = -k*ur
        ezr = k*uz
        pressure_gradient = amplitude*mq/(2q)
        ws.ypsv[1,l],ws.ypsv[2,l],ws.ypsv[3,l] = uz,xi*ezz-la*k*ur,ur
        ws.ypsv[4,l],ws.ypsv[5,l],ws.ypsv[6,l] = 2mu*ezr,pressure,-chi*pressure_gradient
        while receiver <= length(g.lzr) && g.lzr[receiver] == l
            ws.result[1,receiver],ws.result[2,receiver],ws.result[3,receiver] = uz,ur,pressure
            ws.result[4,receiver],ws.result[5,receiver],ws.result[6,receiver] = ezz,radial_strain,ezr
            ws.result[7,receiver] = pressure_gradient
            receiver += 1
        end
    end
    ws.result
end

# Once every propagation factor has underflowed, a buried force plane is
# isolated from all interfaces and other depths. Evaluate its homogeneous
# full-space solution directly instead of solving a disconnected modal system.
# These expressions also preserve its exact odd-in-depth zeros.
function isolated_force_plane!(ws,s,k,sfct)
    c,g,m = ws.input,ws.grid,ws.model
    c.istype == 0 && g.lstop == g.lsbtm && g.lstop > 1 && real(s) >= 0 || return nothing
    l = g.lstop
    g.nno[l-1] == g.nno[l] || return nothing
    all(iszero(sfct[j]) for j in (1,3,4,5,6)) || return nothing
    # exp(-k*h) rounds to zero below half the smallest subnormal Float64.
    k*minimum(view(g.hp,1:length(g.hp)-1)) > 746 || return nothing
    n = g.nno[l]
    mu,xi = m.mu[n],m.la[n]+2m.mu[n]
    alpha,qa,dm = m.alf[n],m.qa[n],m.dm[n]
    q = sqrt(ComplexF64(k*k+s*(1+alpha*qa/xi)/dm))
    a = alpha*qa*s/(xi*dm*(k+q)^2)
    fill!(ws.result,0)
    for (iz,lr) in enumerate(g.lzr)
        lr == l || continue
        ws.result[1,iz] = sfct[2]*(1/mu+(1-a)/xi)/(4k)
        ws.result[6,iz] = sfct[2]*(1-a)/(4xi)
        ws.result[7,iz] = sfct[2]*qa*s/(2xi*dm*(k+q))
    end
    ws.result
end

"""Poroelastic Laplace-wavenumber response; output is reused on the next call.

Rows: uz, ur, pore pressure, ezz, err+ett, ezr, dp/dz.
`sfct` is the six-component source jump in the transformed propagation state.
"""
function kernel!(ws::KernelWorkspace,s,k,sfct)
    k > 0 || throw(ArgumentError("Wavenumber k must be positive"))
    m = ws.model
    if all(iszero,m.alf)
        iszero(sfct[5]) && iszero(sfct[6]) || error("The drained elastic reference accepts mechanical loading only.")
        solve_elastic_modes!(ws,k,sfct)
    else
        isolated = isolated_force_plane!(ws,s,k,sfct)
        isolated === nothing || return isolated
        update_modes!(ws,s,k)
        if ws.input.isurfcon == 0 && length(m.h) == 1 && ws.input.istype in (1,2,3) &&
                all(iszero,view(sfct,1:5))
            return solve_fullspace_fluid!(ws,s,k,sfct)
        end
        solve_coupled_modes!(ws,s,k,sfct)
        if iszero(s) && iszero(sfct[5]) && iszero(sfct[6])
            # A static mechanical force cannot sustain excess pore pressure in
            # a source-free, bounded diffusion solution. Avoid cancellation noise.
            ws.ypsv[5:6,:] .= 0
        end
    end
    for (i,l) in enumerate(ws.grid.lzr)
        n = ws.grid.nno[l]
        la,mu = m.la[n],m.mu[n]
        chi = m.dm[n]*m.alf[n]/m.qa[n]
        y = view(ws.ypsv,:,l)
        ws.result[1,i],ws.result[2,i],ws.result[3,i],ws.result[4,i],
        ws.result[5,i],ws.result[6,i],ws.result[7,i] = (y[1],y[3],y[5],(y[2]+la*k*y[3])/(la+2mu),
            -k*y[3],y[4]/(2mu),iszero(chi) ? 0im : -y[6]/chi)
        # Enforce known exact zeros before relative integration criteria are
        # evaluated; leave the raw state available for residual inspection.
        c = ws.input
        if l == 1 && c.isurfcon > 0
            ws.result[6,i] = 0
            c.isurfcon == 1 && (ws.result[3,i] = 0)
            c.isurfcon == 2 && (ws.result[7,i] = 0)
        end
        if c.isurfcon == 0 && length(m.h) == 1 && c.istype in (1,2,3) &&
           2c.zr[i] == c.zstop+c.zsbtm
            ws.result[1,i] = 0
            ws.result[6,i] = 0
            ws.result[7,i] = 0
        end
    end
    ws.result
end

"""Uniform source-material reference, with the full solution's exact surface zeros.

The same sfct is used in both fields; source pressure is normalized only by the
actual model. Projecting these known zeros makes each integral preserve them.
"""
function direct_kernel!(ws::KernelWorkspace,s,k,sfct)
    reference = ws.direct_workspace
    isnothing(reference) && throw(ArgumentError("The direct-field decomposition requires a layered or half-space fluid-source model."))
    y = kernel!(reference,s,k,sfct)
    if ws.input.isurfcon > 0
        for (iz,l) in enumerate(ws.grid.lzr)
            l == 1 || continue
            y[6,iz] = 0
            ws.input.isurfcon == 1 && (y[3,iz] = 0)
            ws.input.isurfcon == 2 && (y[7,iz] = 0)
        end
    end
    y
end
