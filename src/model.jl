# 材料参数采用 SI；本文件负责物理输入校验、参数转换及沿深度梯度离散。
# LayerModel 最后一层 h=0 表示该层向下无限延伸。
function validate_input(input::PoroelasticInput)
    c = input
    all(isfinite,(c.zstop,c.zsbtm,c.sradius,c.twindow)) ||
        error("Source depths, source radius and time window must be finite.")
    0 <= c.istype <= 3 || error("Source type must be 0 (vertical total force, N), 1 (prescribed drained volume change, m^3; induced correction only), 2 (fluid injection rate, m^3/s), or 3 (mean source pressure * source volume, Pa*m^3).")
    0 <= c.zstop <= c.zsbtm && c.sradius > 0 || error("Invalid source depth/radius.")
    c.istype == 0 || c.zsbtm > c.zstop || error("Source types 1–3 require a finite cylinder.")
    0 <= c.isurfcon <= 2 || error("Surface condition must be 0, 1 or 2.")
    c.nt >= 2 && c.twindow > 0 || error("Require nt ≥ 2 and a positive time window.")
    0 < c.accuracy < 1 || error("Wavenumber accuracy must be between 0 and 1.")
    for (name, x) in (("depth", c.zr), ("distance", c.r))
        !isempty(x) && all(isfinite, x) && minimum(x) >= 0 && issorted(x) ||
            error("Receiver $name samples must be finite, nonnegative and ordered.")
    end
    length(c.ts0) == length(c.sinj0) && !isempty(c.ts0) || error("Empty or unequal source history.")
    first(c.ts0) == 0 && issorted(c.ts0) && all(isfinite, c.ts0) && all(isfinite, c.sinj0) ||
        error("Source history must start at t=0, with nondecreasing finite times/values.")
    length(c.timesn) == length(c.filesn) && all(t -> 0 <= t <= c.twindow, c.timesn) ||
        error("Snapshot times must lie in the output window.")
    nodes = c.model_nodes
    size(nodes, 2) == 6 && size(nodes, 1) >= 1 || error("Model columns: depth, μ, ν, νu, B, D.")
    all(isfinite, nodes) && nodes[1, 1] == 0 && issorted(nodes[:, 1]) || error("Model depths must start at zero and be ordered.")
    for row in eachrow(nodes)
        z, mu, nu, nuu, b, d = row
        z >= 0 && mu > 0 && 0 < nu <= nuu < 0.5 && 0 < b <= 1 && d > 0 ||
            error("Require μ>0, 0<ν≤νu<0.5, 0<B≤1 and D>0.")
        alpha = 3*(nuu-nu)/((1-2nu)*(1+nuu)*b)
        alpha <= 1+8eps(Float64) ||
            error("The Poisson ratios and Skempton B imply a Biot coefficient alpha=$alpha > 1; use consistent material parameters.")
    end
    return c
end

"""Convert physical input parameters and discretize gradients into homogeneous sublayers."""
function build_model(input::PoroelasticInput)
    validate_input(input)
    a = input.model_nodes
    z, mu, nu, nuu, b, d = eachcol(a)
    la = mu .* nu ./ (0.5 .- nu)
    alf = 1.5 .* (nuu .- nu) ./ ((0.5 .- nu) .* (1 .+ nuu) .* b)
    qa = (1 .+ nuu) .* mu .* b ./ (0.5 .- nuu) ./ 3
    dm = d .* (1 .- nuu) .* (0.5 .- nu) ./ ((1 .- nu) .* (0.5 .- nuu))
    values = hcat(la, mu, alf, qa, dm)
    h = Float64[]
    rows = Vector{Vector{Float64}}()
    for i in 1:length(z)-1
        dz = z[i+1] - z[i]
        dz == 0 && continue   # repeated depth gives the two sides of an interface
        v1, v2 = values[i, :], values[i+1, :]
        rel = [v1[j] + v2[j] == 0 ? 0.0 : 2abs(v2[j]-v1[j]) / (v2[j]+v1[j]) for j in 1:5]
        n = round(Int, max(1.0, maximum(rel ./ [0.05, 0.05, 0.05, 0.05, 0.25])), RoundNearestTiesAway)
        for j in 1:n
            push!(h, dz/n)
            push!(rows, v1 + (v2-v1) * ((j-0.5)/n))
        end
    end
    push!(h, 0.0)
    push!(rows, values[end, :])
    p = reduce(vcat, permutedims.(rows))
    # 插值后的每个子层也须对应物理可用的 Skempton 系数。
    for row in eachrow(p)
        la, mu, alpha, qa, _ = row
        b = qa/(la+2mu/3+alpha*qa)
        b <= 1+8eps(Float64) ||
            error("Interpolating the material parameters gives a sublayer with Skempton B=$b > 1; refine the material nodes or specify homogeneous layers.")
    end
    return LayerModel(h, (p[:, j] for j in 1:5)...)
end

function build_grid(input::PoroelasticInput, model::LayerModel)
    tops = [0.0; cumsum(model.h[1:end-1])]
    z = sort!(unique([tops; input.zr; input.zstop; input.zsbtm]))
    hp = [diff(z); 0.0]
    nno = searchsortedlast.(Ref(tops), z)
    lzr = searchsortedfirst.(Ref(z), input.zr)
    lstop, lsbtm = searchsortedfirst(z, input.zstop), searchsortedfirst(z, input.zsbtm)
    if lstop < lsbtm && nno[lstop] != nno[lsbtm-1]
        error("The entire source screen must lie within one homogeneous sublayer.")
    end
    return LayerGrid(z, hp, nno, lzr, lstop, lsbtm)
end
