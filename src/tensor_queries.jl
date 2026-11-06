"""
    strain_tensor(source; r, z, azimuth_deg=0.0, times=nothing, tspan=nothing, method=:exact)
    strain_tensor(source; x, y, z, times=nothing, tspan=nothing, method=:exact)

按物理坐标和时间查询应变，并把局部柱坐标张量转到北、东、下（NED）基底。
`source` 可为完整计算结果或 JLD2 路径；后者只读取所需的四个应变分量及查询元数据。
空间坐标是逐点配对的同形数组，而非新的距离×深度笛卡尔积。`x` 为北向距离、
`y` 为东向距离，`z` 向下为正，单位均为 m。用 `r,z` 时，`azimuth_deg` 为从北
顺时针的径向方位角，单位为度；轴线上径向方向不唯一，约定取 0°。

返回的 `data` 维度为 `3 × 3 × time × point_shape...`，前两轴依次为
`(:north, :east, :down)`。应变无量纲；`ezr` 是张量剪应变，不作工程剪应变的二倍换算。
精确查询与显式插值遵循 `sample_field` 的覆盖范围、材料边界及加载变化时刻限制，张量旋转不会增加空间信息。
所需分量若有 `Inf` 或 `NaN`，会明确报错，避免旋转时出现 `0*Inf` 并生成虚假分量。

例如，`out.data[:, :, it, ip]` 是第 `ip` 个配对点、第 `it` 个所选时刻的 NED
应变张量。曲面法向应变还需用该点的单位法向量 `n` 计算 `dot(n, epsilon*n)`；
曲面上的局部径向应变 `err` 一般不能直接当作曲面法向应变。
"""
function strain_tensor(source; kwargs...)
    context, raw = _query_fields(source, (:err, :ett, :ezz, :ezr))
    plan = _query_plan(context; kwargs...)
    sampled = map(values -> _sample_values(values, plan), raw)
    all(values -> all(isfinite, values), sampled) ||
        throw(ArgumentError("Cannot rotate nonfinite strain samples. Query finite saved times or use a physically finite source rise time."))

    nt, np = length(plan.time), prod(plan.point_shape)
    err, ett, ezz, ezr = map(values -> reshape(values, nt, np), sampled)
    angles = plan.azimuth_deg isa Real ? fill(Float64(plan.azimuth_deg), np) : vec(Float64.(plan.azimuth_deg))
    radii = plan.r isa Real ? fill(Float64(plan.r), np) : vec(plan.r)
    # 只规定轴线上的坐标基底，不改动或平均原始应变分量。
    for ip in 1:np
        iszero(radii[ip]) && (angles[ip] = 0.0)
    end

    data = Array{Float64}(undef, 3, 3, nt, plan.point_shape...)
    tensor = reshape(data, 3, 3, nt, np)
    for ip in 1:np
        c, s = cosd(angles[ip]), sind(angles[ip])
        for it in 1:nt
            rr, tt, zz, rz = err[it,ip], ett[it,ip], ezz[it,ip], ezr[it,ip]
            tensor[1,1,it,ip] = c*c*rr + s*s*tt
            tensor[2,2,it,ip] = s*s*rr + c*c*tt
            tensor[3,3,it,ip] = zz
            tensor[1,2,it,ip] = tensor[2,1,it,ip] = c*s*(rr-tt)
            tensor[1,3,it,ip] = tensor[3,1,it,ip] = c*rz
            tensor[2,3,it,ip] = tensor[3,2,it,ip] = s*rz
        end
    end
    # 返回实际采用的角度，轴线上的 0° 约定也在元数据中明确保留。
    azimuth = isempty(plan.point_shape) ? first(angles) : reshape(angles, plan.point_shape...)
    metadata = merge(_query_metadata(plan), (azimuth_deg=azimuth,))
    (; metadata..., data, unit="1", basis=:ned, components=(:north, :east, :down))
end
