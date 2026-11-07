# JLD2 保存 Float64 场量、输入参数、材料模型、计算网格、加载历史和简短计算状态。
# 各场量分量分别存储，查询一个分量时只需读取对应数组及其坐标和单位。
const RESULT_UNITS = ("m", "m", "1", "1", "1", "1", "rad", "Pa", "m/s", "m/s")
const RESULT_AXIS_ORDER = ("time", "radius", "depth", "component")

# Save only plain Julia records, so archives do not depend on serialized
# MaoyePoro struct definitions.
result_record(value) = NamedTuple{fieldnames(typeof(value))}(
    Tuple(getfield(value, name) for name in fieldnames(typeof(value))))

"""Save a complete result in one lossless JLD2 file; return its absolute path.

All ten fields are saved regardless of the text-output selection. Each field
has its own time × radius × depth dataset, allowing `read_field` to load it
without loading the other fields. `compress=true` enables lossless compression.
Detailed diagnostics are written separately by `write_outputs(...; diagnostics=true)`.
"""
function save_result(path::AbstractString, result; compress::Bool=false)
    input, data = result.input, result.data
    size(data) == (length(result.time), length(input.r), length(input.zr), length(COMPONENTS)) ||
        error("Result data must have shape (time, radius, depth, component).")
    eltype(data) == Float64 || error("Result archives require Float64 fields.")
    destination = abspath(path)
    mkpath(dirname(destination))
    JLD2.jldopen(destination, "w"; compress) do file
        file["schema_version"] = 2
        file["package_version"] = string(pkgversion(@__MODULE__))
        file["julia_version"] = string(VERSION)
        file["axis_order"] = collect(RESULT_AXIS_ORDER)
        file["field_axis_order"] = collect(RESULT_AXIS_ORDER[1:3])
        file["components"] = string.(collect(COMPONENTS))
        file["units"] = collect(RESULT_UNITS)
        file["coordinate_units"] = (time="s", radius="m", depth="m")
        file["source_value_unit"] = ("N", "m^3", "m^3/s", "N*m")[input.istype+1]
        file["depth_positive"] = "down"
        file["time"] = result.time
        file["radii_m"] = input.r
        file["depths_m"] = input.zr
        file["input"] = result_record(input)
        file["model"] = result_record(result.model)
        file["grid"] = result_record(result.grid)
        file["status"] = result.status
        file["equivalent"] = result.equivalent
        file["equivalent_impulses"] = result.equivalent_impulses
        file["settings"] = hasproperty(result, :settings) ? result.settings : nothing
        for (ic, name) in enumerate(COMPONENTS)
            file["fields/$(name)"] = data[:, :, :, ic]
        end
    end
    destination
end

function check_result_archive(file)
    file["schema_version"] == 2 || error("Unsupported MaoyePoro result schema.")
    file["components"] == string.(collect(COMPONENTS)) || error("Unexpected result component order.")
    file["axis_order"] == collect(RESULT_AXIS_ORDER) || error("Unexpected result axis order.")
    file["field_axis_order"] == collect(RESULT_AXIS_ORDER[1:3]) || error("Unexpected field axis order.")
    file["units"] == collect(RESULT_UNITS) || error("Unexpected result units.")
end

"""Load physical results, input, model/grid and brief numerical status from JLD2.

Detailed diagnostics are separate from the result archive. The returned
`diagnostics` is empty and `equivalent_diagnostics` is nothing.
"""
function load_result(path::AbstractString)
    JLD2.jldopen(path, "r") do file
        check_result_archive(file)
        input = PoroelasticInput(; file["input"]...)
        m, g = file["model"], file["grid"]
        model = LayerModel(m.h, m.la, m.mu, m.alf, m.qa, m.dm)
        grid = LayerGrid(g.z, g.hp, g.nno, g.lzr, g.lstop, g.lsbtm)
        time = file["time"]
        isequal(file["radii_m"], input.r) && isequal(file["depths_m"], input.zr) ||
            error("Result coordinates disagree with the saved input.")
        shape = (length(time), length(input.r), length(input.zr))
        data = Array{Float64}(undef, shape..., length(COMPONENTS))
        for (ic, name) in enumerate(COMPONENTS)
            values = file["fields/$(name)"]
            size(values) == shape && eltype(values) == Float64 ||
                error("Unexpected dimensions or precision for field $name.")
            copyto!(view(data, :, :, :, ic), values)
        end
        (; input, model, grid, time, data, status=file["status"], diagnostics=NamedTuple[],
            equivalent=file["equivalent"], equivalent_impulses=file["equivalent_impulses"], equivalent_diagnostics=nothing,
            settings=file["settings"])
    end
end

"""Read one archived field and its coordinates/unit without loading all ten fields."""
function read_field(path::AbstractString, component)
    ic = component_index(component)
    JLD2.jldopen(path, "r") do file
        check_result_archive(file)
        time, r, z = file["time"], file["radii_m"], file["depths_m"]
        data = file["fields/$(COMPONENTS[ic])"]
        size(data) == (length(time), length(r), length(z)) && eltype(data) == Float64 ||
            error("Unexpected dimensions or precision for field $(COMPONENTS[ic]).")
        (; time, r, z, component=COMPONENTS[ic], unit=RESULT_UNITS[ic], data)
    end
end
