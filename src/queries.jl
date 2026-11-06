# 查询已保存的轴对称场，使用结果中记录的计算坐标、材料模型和源半径。
# 查询坐标逐点配对并按 Julia 广播规则组合；与输入中 r×z 的全组合网格不同。

function _history_break_times(input)
    times, values = input.ts0, input.sinj0
    breaks = Float64[]
    left_slope = 0.0
    i = 1
    while i <= length(times) && times[i] <= input.twindow
        j = i
        while j < length(times) && times[j+1] == times[i]
            j += 1
        end
        jump = values[j] - (i == 1 ? 0.0 : values[i])
        right_slope = j == length(times) ? 0.0 : (values[j+1]-values[j])/(times[j+1]-times[j])
        # 只禁止跨实际跳变或斜率转折插值；共线冗余节点不影响插值。
        # 窗口末端的连续斜坡采用左导数，窗外斜率改变不影响本窗口。
        if !iszero(jump) || (times[i] < input.twindow && left_slope != right_slope)
            push!(breaks, times[i])
        end
        left_slope = right_slope
        i = j+1
    end
    breaks
end

function _check_query_time_interval(input, lower, upper, breaks=_history_break_times(input))
    # 上端点包含在检查内：节点后的响应不能被插回该节点之前。
    i = searchsortedlast(breaks, lower)+1
    if i <= length(breaks) && breaks[i] <= upper
        error("Time interpolation crosses a source-history change at $(breaks[i]) s. Use saved times or recompute a finer time grid.")
    end
end

function _axis_bracket(axis, value, method, label)
    isfinite(value) || error("Query $label must be finite.")
    first(axis) <= value <= last(axis) || error("Query $label=$value is outside the saved range [$(first(axis)), $(last(axis))]; extrapolation is not supported.")
    i = searchsortedfirst(axis, value)
    if i <= length(axis) && axis[i] == value
        return (lo=i, hi=i, weight=0.0)
    end
    method == :linear || error("Query $label=$value is not a saved node. Use method=:linear explicitly for interpolation, or compute this coordinate.")
    (lo=i-1, hi=i, weight=(value-axis[i-1])/(axis[i]-axis[i-1]))
end

function _coordinate_triplet(a, b, c)
    all(v -> v isa Real || v isa AbstractArray, (a,b,c)) ||
        error("Query coordinates must be real scalars or arrays; use arrays for point sets.")
    # 按相同广播规则分别生成三个坐标数组，使其形状一致。
    (broadcast((x,y,z) -> Float64(x), a,b,c),
     broadcast((x,y,z) -> Float64(y), a,b,c),
     broadcast((x,y,z) -> Float64(z), a,b,c))
end

_flat_coordinates(value::Real) = [Float64(value)]
_flat_coordinates(value) = vec(value)

function _query_times(context, times, tspan)
    times === nothing || tspan === nothing || error("Use times or tspan, not both.")
    if tspan !== nothing
        length(tspan) == 2 || error("tspan must contain a start and an end time in seconds.")
        lower, upper = Float64.(tspan)
        first(context.time) <= lower <= upper <= last(context.time) || error("tspan must lie inside the saved time window.")
        selected = context.time[findall(t -> lower <= t <= upper, context.time)]
        isempty(selected) && error("tspan contains no saved time samples. Specify times and an explicit interpolation method if needed.")
        return selected
    end
    selected = times === nothing ? copy(context.time) : times isa Real ? [Float64(times)] : vec(Float64.(collect(times)))
    !isempty(selected) || error("The requested time list is empty.")
    selected
end

function _check_spatial_cell(context, rb, zb, effective_radius, layer_tops)
    c = context.input
    rl, rh = c.r[rb.lo], c.r[rb.hi]
    zl, zh = c.zr[zb.lo], c.zr[zb.hi]
    if zb.lo != zb.hi
        # 界面节点保存从下侧接近界面时的极限值，不能用它插值上层应变；梯度离散出的子层也在内。
        searchsortedlast(layer_tops, zl) == searchsortedlast(layer_tops, zh) ||
            error("Spatial interpolation crosses a material sublayer boundary. Query saved nodes or compute samples on the required side.")
        if rl < effective_radius && any(z -> zl <= z <= zh, (c.zstop,c.zsbtm))
            error("Spatial interpolation touches a source-support boundary. Compute this location instead of interpolating across the source.")
        end
    end
    if rb.lo != rb.hi && rl <= effective_radius <= rh && zh >= c.zstop && zl <= c.zsbtm
        error("Spatial interpolation touches the radial source-support boundary. Compute this location instead.")
    end
end

function _query_plan(context; r=nothing, z=nothing, x=nothing, y=nothing,
                     azimuth_deg=nothing, times=nothing, tspan=nothing, method=:exact)
    method = Symbol(method)
    method in (:exact,:linear) || error("Query method must be :exact or :linear.")
    z === nothing && error("Provide query depth z in metres, positive downward.")
    if r !== nothing
        x === nothing && y === nothing || error("Use r or x/y coordinates, not both.")
        R,Z,A = _coordinate_triplet(r,z,azimuth_deg === nothing ? 0.0 : azimuth_deg)
        A = ifelse.(iszero.(R), 0.0, A)
        X,Y = R.*cosd.(A), R.*sind.(A)
    else
        x !== nothing && y !== nothing || error("Provide r/z or x/y/z query coordinates.")
        azimuth_deg === nothing || error("x/y already determine azimuth; do not also supply azimuth_deg.")
        X,Y,Z = _coordinate_triplet(x,y,z)
        R = hypot.(X,Y)
        A = atan.(Y,X).*(180/pi)
        A = ifelse.(iszero.(R), 0.0, A)
    end
    all(v -> all(isfinite,_flat_coordinates(v)), (R,Z,A,X,Y)) || error("Query coordinates and azimuth must be finite.")
    all(>=(0),_flat_coordinates(R)) && all(>=(0),_flat_coordinates(Z)) || error("Query radius and depth must be nonnegative.")
    shape = R isa Real ? () : size(R)
    requested_times = _query_times(context,times,tspan)
    rb = [_axis_bracket(context.input.r,v,method,"radius [m]") for v in _flat_coordinates(R)]
    zb = [_axis_bracket(context.input.zr,v,method,"depth [m]") for v in _flat_coordinates(Z)]
    tb = [_axis_bracket(context.time,t,method,"time [s]") for t in requested_times]
    time_interpolated = any(b -> b.lo != b.hi,tb)
    if time_interpolated
        breaks = _history_break_times(context.input)
        for b in tb
            b.lo == b.hi || _check_query_time_interval(context.input,context.time[b.lo],context.time[b.hi],breaks)
        end
    end
    spatial_interpolated = any(b -> b.lo != b.hi,rb) || any(b -> b.lo != b.hi,zb)
    effective_radius = isempty(context.diagnostics) ? context.input.sradius : first(context.diagnostics).effective_radius
    if spatial_interpolated
        layer_tops = [0.0; cumsum(context.model.h[1:end-1])]
        for (radial,depth) in zip(rb,zb)
            _check_spatial_cell(context,radial,depth,effective_radius,layer_tops)
        end
    end
    (time=requested_times,r=R,z=Z,x=X,y=Y,azimuth_deg=A,point_shape=shape,
        rb=rb,zb=zb,tb=tb,method=method,spatial_interpolated=spatial_interpolated,
        time_interpolated=time_interpolated,effective_radius=effective_radius,
        source_diagnostics=context.diagnostics,
        saved_shape=(length(context.time),length(context.input.r),length(context.input.zr)))
end

function _sample_values(values, plan)
    size(values) == plan.saved_shape || error("Field dimensions disagree with the saved coordinate axes.")
    eltype(values) == Float64 || error("Field queries require the saved Float64 precision.")
    data = Matrix{Float64}(undef,length(plan.time),length(plan.rb))
    for ip in eachindex(plan.rb), it in eachindex(plan.tb)
        rb,zb,tb = plan.rb[ip],plan.zb[ip],plan.tb[it]
        if rb.lo == rb.hi && zb.lo == zb.hi && tb.lo == tb.hi
            # 精确节点直接复制，保留 Inf、NaN 及有符号零，避免 0*Inf。
            data[it,ip] = values[tb.lo,rb.lo,zb.lo]
            continue
        end
        value = 0.0
        for (jt,wt) in ((tb.lo,1-tb.weight),(tb.hi,tb.weight)),
            (jr,wr) in ((rb.lo,1-rb.weight),(rb.hi,rb.weight)),
            (jz,wz) in ((zb.lo,1-zb.weight),(zb.hi,zb.weight))
            weight = wt*wr*wz
            iszero(weight) && continue
            sample = values[jt,jr,jz]
            isfinite(sample) || error("Cannot interpolate a requested field across nonfinite samples.")
            value += weight*sample
        end
        data[it,ip] = value
    end
    reshape(data,(length(plan.time),plan.point_shape...))
end

function _query_metadata(plan)
    (time=plan.time,r=plan.r,z=plan.z,x=plan.x,y=plan.y,azimuth_deg=plan.azimuth_deg,
        method=plan.method,spatial_interpolated=plan.spatial_interpolated,
        time_interpolated=plan.time_interpolated,effective_radius=plan.effective_radius,
        source_diagnostics=plan.source_diagnostics)
end

_query_fields(source, components) = (source, map(c -> field(source,c),components))

function _query_fields(path::AbstractString, components)
    # 查询归档时仅加载所需分量，并读取插值所需的材料/源历史元数据。
    JLD2.jldopen(path,"r") do file
        check_result_archive(file)
        input = PoroelasticInput(;file["input"]...)
        m = file["model"]
        model = LayerModel(m.h,m.la,m.mu,m.alf,m.qa,m.dm)
        isequal(file["radii_m"],input.r) && isequal(file["depths_m"],input.zr) || error("Saved query coordinates disagree with the input.")
        context = (input=input,model=model,time=file["time"],diagnostics=file["diagnostics"])
        context,map(c -> file["fields/$(c)"],components)
    end
end

"""
    sample_field(result_or_path, component; r, z, times=nothing, tspan=nothing, method=:exact)
    sample_field(result_or_path, component; x, y, z, ...)

按物理坐标查询已计算的位移、应变、孔压等量。坐标用米，时间用秒；x=北、y=东、z=下，原点在源轴线。
r/z 或 x/y/z 使用 Julia 广播规则逐点配对，同形矩阵可表示剖面或曲面。
返回 data 的顺序为 time × 查询坐标形状；单点返回时间向量，times 为标量时仍保留时间维。

默认 :exact 只取已保存坐标和时刻。:linear 显式进行线性插值，不外推、不跨材料层或
注水/加载区域的边界，也不跨源随时间的跳变或斜率转折；所用端点不能含 Inf/NaN。
插值误差不包括在原求解诊断中。
tspan=(起点,终点) 仅筛选区间内已有时间样点，不生成端点，也不是时间平均。

十个场量分量采用局部柱坐标基底；固定方向应变请用 strain_tensor。
r/z 查询可另给 azimuth_deg，方位从北顺时针；轴线上固定为 0°。
此函数从已有结果取值；指定 JLD2 路径时只读取所需的场量分量和查询元数据。
"""
function sample_field(source, component; kwargs...)
    ic = component_index(component)
    context,values = _query_fields(source,(COMPONENTS[ic],))
    plan = _query_plan(context;kwargs...)
    (;_query_metadata(plan)...,component=COMPONENTS[ic],unit=RESULT_UNITS[ic],
        basis=:cylindrical,data=_sample_values(only(values),plan))
end
