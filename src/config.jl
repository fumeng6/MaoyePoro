const SOURCE_KINDS = ("force", "matrix_volume", "injection_rate", "injection_moment")
const SOURCE_DESCRIPTIONS = ("竖向总力 [N]", "已知排水体积应变在源区的积分 [m^3]；只计算孔压与附加变形，不含原排水弹性场",
    "每秒注入的流体体积，负值表示抽水 [m^3/s]", "源区平均孔压 × 源区体积 [Pa*m^3 = N*m]")
const SURFACE_CONDITIONS = ("fullspace", "drained", "sealed")
const MODEL_NODE_KEYS = ("depth_m", "shear_modulus_pa", "poisson",
    "undrained_poisson", "skempton", "diffusivity_m2_s")

function check_config_keys(table, allowed, section)
    table isa AbstractDict || throw(ArgumentError("$section must be a TOML table."))
    unknown = setdiff(collect(keys(table)), allowed)
    isempty(unknown) || throw(ArgumentError("Unknown $section key(s): $(join(sort!(unknown), ", "))."))
    table
end

function config_code(value, names, label)
    index = findfirst(==(value), names)
    isnothing(index) && throw(ArgumentError("$label must be one of $(join(names, ", ")); got $(repr(value))."))
    index-1
end

"""Read named simulation parameters from a TOML file.

Histories and model nodes retain their original order, including repeated
times/depths and future history nodes. A nonzero initial value is a step at t=0;
the reader never inserts history nodes. All numerical values use SI units.
"""
function read_input(path::AbstractString)
    lowercase(splitext(path)[2]) == ".toml" ||
        throw(ArgumentError("Simulation input must be a .toml file with named parameters."))
    read_toml_input(path)
end

history_path(input_path, filename) = isabspath(filename) ? normpath(filename) :
    normpath(joinpath(dirname(abspath(input_path)), filename))

# Read time_s and value from a two-column CSV; time is in seconds and source.kind sets the value unit.
function read_source_history(path)
    times, values = Float64[], Float64[]
    open(path) do io
        header_seen = false
        for (line_number, raw) in enumerate(eachline(io))
            line_number == 1 && (raw = replace(raw, '\ufeff' => ""))
            line = strip(raw)
            (isempty(line) || startswith(line, '#')) && continue
            tokens = map(split(line, ','; keepempty=true)) do token
                value = strip(token)
                startswith(value, '"') && endswith(value, '"') && length(value) >= 2 ?
                    chop(value; head=1, tail=1) : value
            end
            if !header_seen
                tokens == ["time_s", "value"] ||
                    throw(ArgumentError("$path:$line_number: source history header must be time_s,value."))
                header_seen = true
                continue
            end
            length(tokens) == 2 || throw(ArgumentError("$path:$line_number: source history requires exactly two columns."))
            t, v = map(token -> tryparse(Float64, token), tokens)
            t !== nothing && v !== nothing && isfinite(t) && isfinite(v) ||
                throw(ArgumentError("$path:$line_number: source history requires two finite numeric values."))
            push!(times, t); push!(values, v)
        end
        header_seen && !isempty(times) || throw(ArgumentError("$path: source history contains no data rows."))
    end
    times, values
end

function source_history(source, input_path)
    if haskey(source, "history_file")
        !haskey(source, "times_s") && !haskey(source, "values") ||
            throw(ArgumentError("Use source.history_file or source.times_s and source.values, not both."))
        filename = source["history_file"]
        filename isa AbstractString || throw(ArgumentError("source.history_file must be a file path."))
        return read_source_history(history_path(input_path, filename))
    end
    haskey(source, "times_s") && haskey(source, "values") ||
        throw(ArgumentError("Provide source.history_file or both source.times_s and source.values."))
    Float64.(source["times_s"]), Float64.(source["values"])
end

function read_toml_input(path)
    config = TOML.parsefile(path)
    check_config_keys(config, ("schema_version", "source", "receivers", "time",
        "medium", "numerics", "output"), "root")
    version = get(config, "schema_version", 2)
    version isa Integer && version == 2 || throw(ArgumentError("Unsupported simulation input schema_version: $version; use schema_version = 2 with named parameter tables."))
    source = check_config_keys(config["source"],
        ("kind", "top_m", "bottom_m", "radius_m", "times_s", "values", "history_file"), "source")
    receivers = check_config_keys(config["receivers"], ("radii_m", "depths_m"), "receivers")
    time = check_config_keys(config["time"], ("window_s", "samples"), "time")
    medium = check_config_keys(config["medium"], ("surface", "nodes"), "medium")
    numerics = check_config_keys(get(config, "numerics", Dict()), ("accuracy",), "numerics")
    output = check_config_keys(get(config, "output", Dict()), ("fields", "snapshot_times_s",
        "snapshot_files", "compute_equivalent", "field_files", "equivalent_filename", "diagnostics"), "output")
    diagnostics = get(output, "diagnostics", false)
    diagnostics isa Bool || throw(ArgumentError("output.diagnostics must be true or false."))

    istype = config_code(source["kind"], SOURCE_KINDS, "source.kind")
    ts0, sinj0 = source_history(source, path)
    isurfcon = config_code(medium["surface"], SURFACE_CONDITIONS, "medium.surface")
    node_records = medium["nodes"]
    node_records isa AbstractVector || throw(ArgumentError("Use [[medium.nodes]] for the ordered model nodes."))
    nodes = Matrix{Float64}(undef, length(node_records), 6)
    for (i, record) in enumerate(node_records)
        check_config_keys(record, MODEL_NODE_KEYS, "medium.nodes[$i]")
        for j in 1:6
            nodes[i,j] = record[MODEL_NODE_KEYS[j]]
        end
    end
    fields = String.(get(output, "fields", string.(COMPONENTS)))
    unknown_fields = setdiff(fields, string.(COMPONENTS))
    isempty(unknown_fields) || throw(ArgumentError("Unknown output fields: $(join(unknown_fields, ", "))."))
    selected = ntuple(i -> string(COMPONENTS[i]) in fields, 10)
    field_files = check_config_keys(get(output, "field_files", Dict()),
        string.(COMPONENTS), "output.field_files")
    filenames = ntuple(i -> String(get(field_files, string(COMPONENTS[i]), "$(COMPONENTS[i]).dat")), 10)
    timesn = Float64.(get(output, "snapshot_times_s", Float64[]))
    filesn = String.(get(output, "snapshot_files",
        ["snapshot_$(lpad(string(i),4,'0')).dat" for i in eachindex(timesn)]))
    input = PoroelasticInput(
        zstop=Float64(source["top_m"]), zsbtm=Float64(source["bottom_m"]),
        sradius=Float64(source["radius_m"]), istype=istype,
        ts0=ts0, sinj0=sinj0,
        zr=Float64.(receivers["depths_m"]), r=Float64.(receivers["radii_m"]),
        twindow=Float64(time["window_s"]), nt=Int(time["samples"]),
        accuracy=Float64(get(numerics, "accuracy", 0.05)),
        selected=selected, filenames=filenames, timesn=timesn, filesn=filesn,
        compute_equivalent=Bool(get(output, "compute_equivalent", false)),
        equivalent_filename=String(get(output, "equivalent_filename", "equivalent_source.dat")),
        diagnostics=diagnostics,
        isurfcon=isurfcon, model_nodes=nodes)
    validate_input(input)
end

"""Write all input parameters as named TOML records without changing their values.

With `history_file="history.csv"`, write the source history to that CSV and
reference it from TOML. A relative path is relative to the destination TOML.
The default keeps the history inline.
"""
function write_input(path::AbstractString, input::PoroelasticInput; history_file=nothing)
    lowercase(splitext(path)[2]) == ".toml" || throw(ArgumentError("write_input requires a .toml destination."))
    c = validate_input(input)
    config = Dict(
        "schema_version" => 2,
        "source" => Dict("kind" => SOURCE_KINDS[c.istype+1], "top_m" => c.zstop,
            "bottom_m" => c.zsbtm, "radius_m" => c.sradius, "times_s" => c.ts0, "values" => c.sinj0),
        "receivers" => Dict("radii_m" => c.r, "depths_m" => c.zr),
        "time" => Dict("window_s" => c.twindow, "samples" => c.nt),
        "medium" => Dict("surface" => SURFACE_CONDITIONS[c.isurfcon+1],
            "nodes" => [Dict(MODEL_NODE_KEYS[j] => c.model_nodes[i,j] for j in 1:6)
                for i in axes(c.model_nodes,1)]),
        "numerics" => Dict("accuracy" => c.accuracy),
        "output" => Dict("fields" => [string(COMPONENTS[i]) for i in 1:10 if c.selected[i]],
            "field_files" => Dict(string(COMPONENTS[i]) => c.filenames[i] for i in 1:10), "snapshot_times_s" => c.timesn,
            "snapshot_files" => c.filesn, "compute_equivalent" => c.compute_equivalent,
            "equivalent_filename" => c.equivalent_filename, "diagnostics" => c.diagnostics))
    if history_file !== nothing
        history_file isa AbstractString || throw(ArgumentError("history_file must be a file path."))
        destination = history_path(path, history_file)
        same_path = Sys.iswindows() ? lowercase(destination) == lowercase(abspath(path)) : destination == abspath(path)
        same_path && throw(ArgumentError("The source history CSV and TOML must have different paths."))
        config["source"]["history_file"] = String(history_file)
        delete!(config["source"], "times_s")
        delete!(config["source"], "values")
        mkpath(dirname(destination))
        open(destination, "w") do io
            println(io, "time_s,value")
            for (t, value) in zip(c.ts0, c.sinj0)
                @printf(io, "%.16e,%.16e\n", t, value)
            end
        end
    end
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# MaoyePoro 输入：SI 单位，距离 m、时间 s、压力/模量 Pa，深度向下为正。")
        println(io, "# source.kind 选择物理过程，values 的含义和单位如下：")
        for (kind, description) in zip(SOURCE_KINDS, SOURCE_DESCRIPTIONS)
            println(io, "#   ", kind, "：", description)
        end
        println(io, "# source 的 top_m、bottom_m 和 radius_m 给出源区顶深度、底深度和半径，单位 m；源强度沿径向呈抛物分布。")
        println(io, "# times_s 与 values 一一对应，首时刻为 0；非零首值表示零时刻阶跃，相邻不同时刻之间线性变化。")
        println(io, "# 同一时刻的相邻值表示阶跃前后的状态；history_file 是相对本 TOML 的 CSV 路径，也可填写绝对路径。")
        println(io, "# receivers 的 radii_m 与 depths_m 作全部组合，给出每个计算位置的径向距离和深度，单位 m。")
        println(io, "# time.samples 是从 0 到 time.window_s 等间隔保存的时间点数，包含首末点。")
        println(io, "# medium.surface：fullspace 无自由面；drained 表面孔压为零；sealed 表面法向 Darcy 通量为零。")
        println(io, "# medium.nodes 给出材料随深度的变化：深度 m、剪切模量 Pa、排水和不排水泊松比、Skempton B、水力扩散系数 m²/s。")
        println(io, "# 同一深度的相邻材料节点分别表示界面上下两侧；不同深度间的材料参数按梯度离散为均质子层。")
        println(io, "# numerics.accuracy 控制波数和频率采样加密时允许的相对变化，简短收敛状态见结果 status。")
        println(io, "# output.fields 选择要导出的文本时序分量；JLD2 保存全部十个场量分量。快照时刻单位 s，须在输出窗口内。")
        println(io, "# output.diagnostics=true 或命令行 --diagnostics 单独保存 diagnostics.txt，不改变计算精度。")
        TOML.print(io, config; sorted=true)
    end
    abspath(path)
end
