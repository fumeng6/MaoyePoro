"""Homogeneous sublayers; the final zero thickness denotes the lower half-space."""
struct LayerModel
    h::Vector{Float64}
    la::Vector{Float64}
    mu::Vector{Float64}
    alf::Vector{Float64}
    qa::Vector{Float64}
    dm::Vector{Float64}
end

"""Depths at material interfaces, source endpoints and receivers; interface values use the limit from below."""
struct LayerGrid
    z::Vector{Float64}
    hp::Vector{Float64}
    nno::Vector{Int}
    lzr::Vector{Int}
    lstop::Int
    lsbtm::Int
end

"""孔隙弹性计算输入，全部采用 SI，深度向下为正。

istype=0：地表或地下施加的竖向总力，N。
istype=1：给定排水体积应变的体积积分，m³；只计算它引起的孔压与附加变形，不含原排水弹性场。
istype=2：每秒注入的流体体积，m³/s；负值表示抽水。
istype=3：源区平均孔压乘以源区体积，Pa·m³ = N·m。
ts0/sinj0 记录上述物理量随时间的变化；重复时刻表示阶跃，其他节点间线性变化。
r 和 zr 是内部计算轴，求解全部 r×zr 组合，不是逐点配对的坐标表。
model_nodes 每行六列：深度 m、剪切模量 Pa、排水/不排水泊松比、Skempton B、扩散系数 m²/s。

已有结果按任意形状的配对坐标取样请用 sample_field；北、东、下（NED）坐标下的应变用 strain_tensor。
参数扫描用 PoroelasticInput(input; kwargs...) 复制参数，避免原地改动已有结果的输入记录。
"""
Base.@kwdef struct PoroelasticInput
    # 注水、施力或体积变化的区域：沿深度均匀、径向抛物分布。整个井段须在一个实际均质子层内。
    zstop::Float64
    zsbtm::Float64
    sradius::Float64
    istype::Int
    # 时间单位 s；值的单位由 istype 决定。保留重复时刻和窗外节点，不自动抽稀。
    ts0::Vector{Float64}
    sinj0::Vector{Float64}
    # 计算轴：径向距离、接收深度均为 m，非负有序；二者作全部组合。
    zr::Vector{Float64}
    r::Vector{Float64}
    # 输出时间窗 s、等距时间样点数。nt 不等于源时间记录的行数，也不直接指定内部计算的频率数。
    twindow::Float64
    nt::Int
    accuracy::Float64 = 0.05
    # 以下选项控制文件输出和等效源计算。selected 选择要导出的文本时序分量，求解仍计算全部十个分量。
    selected::NTuple{10,Bool} = ntuple(_ -> true, 10)
    filenames::NTuple{10,String} = ("uz.dat", "ur.dat", "ezz.dat", "err.dat", "ett.dat",
                                   "ezr.dat", "tlt.dat", "pp.dat", "dvz.dat", "dvr.dat")
    timesn::Vector{Float64} = Float64[]
    filesn::Vector{String} = String[]
    compute_equivalent::Bool = false
    equivalent_filename::String = "equivalent_source.dat"
    # 顶部：0 全空间；1 排水自由面；2 封闭自由面。材料重复深度表示界面两侧。
    isurfcon::Int = 1
    model_nodes::Matrix{Float64}
end

"""Copy an input with keyword changes, keeping its array parameters independent."""
function PoroelasticInput(input::PoroelasticInput; kwargs...)
    original = NamedTuple{fieldnames(PoroelasticInput)}(
        Tuple(getfield(input, name) for name in fieldnames(PoroelasticInput)))
    parameters = map(x -> x isa AbstractArray ? copy(x) : x, merge(original, (;kwargs...)))
    validate_input(PoroelasticInput(;parameters...))
end

const COMPONENTS = (:uz, :ur, :ezz, :err, :ett, :ezr, :tlt, :pp, :dvz, :dvr)
