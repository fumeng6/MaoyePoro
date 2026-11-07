const PROGRAM_VERSION = v"1.0.0"

const _CLI_INPUT_TEMPLATE = raw"""
# MaoyePoro 待填写输入配置；填写完成前不能计算。
# 标注“必填”的空字符串 "" 和空数组 [] 须填写。
# 数值字段：将 "" 整体替换为数字，去掉引号；kind、surface 仍填写带引号的字符串。
# 全部采用 SI：距离 m、时间 s、压力及模量 Pa；深度 z 向下为正。
# 填写后先运行：myporo case.toml --check（将 case.toml 换成此文件名）。

schema_version = 2 # 输入文件结构版本，请保留。

[source]
# 必填：选择物理过程；values 的含义及单位随 kind 改变：
# "injection_rate"：每秒注入的流体体积，m³/s；正值注水，负值抽水。
# "force"：竖向总力，N；正值向下。
# "matrix_volume"：排水体积应变在源区的积分，m³；只计算孔压与附加变形，不含原排水弹性场。
# "injection_moment"：源区平均孔压 × 源区体积，Pa·m³；不是井壁某点的压力。
kind = ""

# 必填：源区顶、底深度及物理半径，m；整个源区须位于同一实际均质子层内。
# 0 <= top_m <= bottom_m，radius_m > 0；除 force 外，bottom_m 必须大于 top_m。
# 径向分布权重为抛物形，轴线上最大、外缘为零。
top_m = ""
bottom_m = ""
radius_m = ""

# 必填：时间 s 与源值一一对应，数组非空且等长；values 单位见 kind。
# 首时刻为 0，时间非递减；不同时间节点间线性变化，重复时刻表示阶跃，不可去重。
times_s = []
values = []
# 长历史可改用下面的 CSV：删除 times_s、values 两行后再取消注释。
# CSV 首行为 time_s,value；相对路径以本 TOML 所在目录为准。
# history_file = "history.csv"

[receivers]
# 必填：径向距离与深度，m；两个数组均非空、非负且有序。
# 程序计算全部“半径 × 深度”组合，不是逐项配对的坐标。
radii_m = []
depths_m = []

[time]
window_s = "" # 必填：从 0 计算至此时刻，s，须大于 0。
samples = "" # 必填：均匀保存的时间点数，整数 >= 2，含首末点；不是源历史节点数。

[medium]
# 必填：按实际边界选择：
# "fullspace"：全空间，无自由表面；上方按首层材料延伸，z=0 只是深度参考。
# "drained"：排水地表，孔压为零。
# "sealed"：封闭地表，法向 Darcy 通量为零。
surface = ""

# 首材料节点深度必须为 0；单节点表示均匀介质，末节点以下无限延伸。
# 分层时复制整个 [[medium.nodes]] 段并按深度排列。
# 突变界面须重复同一深度，分别给出上下两侧材料；不同深度间的物性变化解释为梯度。
[[medium.nodes]]
depth_m = 0.0
shear_modulus_pa = "" # 必填：剪切模量 μ，Pa，须大于 0。
poisson = "" # 必填：排水泊松比 ν，无量纲。
undrained_poisson = "" # 必填：不排水泊松比 νu，无量纲；通常耦合条件为 0 < ν < νu < 0.5。
skempton = "" # 必填：Skempton B，无量纲，0 < B <= 1；与泊松比须满足 Biot 系数 <= 1。
diffusivity_m2_s = "" # 必填：水力扩散系数 D，m²/s，须大于 0；不是渗透率。

[numerics]
# 可保留默认值：控制波数积分与拉普拉斯逆变换的数值收敛判据，须在 0 与 1 之间。
# 不是每个场量、每个样点的相对误差上限；计算后查看结果 status。
accuracy = 0.05

[output]
# 默认保存全部十类场量的 result.jld2；命令行 --format text 或 both 可导出文本。
# 开启后单独保存 diagnostics.txt 并显示详细进度；不改变计算和收敛检查。
diagnostics = false
# 可选：fields 只选择文本时序分量，不改变计算和 JLD2 保存的分量。
# fields = ["uz", "pp"]
# 可选：快照时刻，s，须在 0 到 window_s 内；空数组表示不保存快照。
snapshot_times_s = []
"""

function cli_main(args)
    if "--version" in args
        println("MaoyePoro ",PROGRAM_VERSION," (myporo)")
        return
    end
    if isempty(args) || "--help" in args || "-h" in args
        println("""
        MaoyePoro: 注水、抽水和竖向荷载引起的孔压扩散与地层变形
        Usage: myporo case.toml [options] [--julia-args -t4]
          or: myporo --template [case.toml]
          or: julia --threads=4 --project=. maoyeporo.jl case.toml [options]
          --template [FILE]     Write an input template (default: case.toml); never overwrite
          --output DIR, -o DIR   Output directory
          --format FORMAT       jld2, text, or both (default: jld2)
          --max-frequencies N   Frequency capacity, a power of two >= 512
          --image-factor X      Suppress time-period aliasing (default: 0.001)
          --serial              Disable frequency parallelism
          --diagnostics         Save diagnostics.txt and show detailed progress
          --quiet               Suppress progress messages
          --check               Read and check input/model/grid without solving
          --version             Show program version
          --help, -h            Show this help
        """)
        return
    end
    if first(args) == "--template"
        length(args) <= 2 || error("Usage: myporo --template [case.toml]; do not combine with solver options.")
        path = length(args) == 2 ? args[2] : "case.toml"
        !startswith(path, "-") && lowercase(splitext(path)[2]) == ".toml" ||
            error("Template destination must be a .toml file; use myporo --template [case.toml].")
        (ispath(path) || islink(path)) && error("模板目标已存在，拒绝覆盖：$path")
        write(path, _CLI_INPUT_TEMPLATE)
        println("已写入待填写模板：", abspath(path))
        println("请先填写标注“必填”的参数；数值字段去掉空字符串的引号。填写后使用 --check 检查输入。")
        return
    end
    options = Dict{Symbol,Any}()
    positional = String[]
    check_only = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("--output", "-o", "--format", "--max-frequencies", "--image-factor")
            i += 1
            i <= length(args) || error("Missing value after $arg.")
            value = args[i]
            startswith(value, "-") && error("Missing value after $arg; got option $value.")
            if arg in ("--output", "-o")
                options[:output_dir] = value
            elseif arg == "--format"
                Symbol(value) in (:jld2, :text, :both) || error("Output format must be jld2, text, or both.")
                options[:format] = Symbol(value)
            elseif arg == "--max-frequencies"
                options[:max_frequencies] = parse(Int, value)
            else
                options[:image_factor] = parse(Float64, value)
            end
        elseif arg == "--serial"
            options[:parallel] = false
        elseif arg == "--diagnostics"
            options[:diagnostics] = true
        elseif arg == "--quiet"
            options[:verbose] = false
        elseif arg == "--check"
            check_only = true
        elseif startswith(arg, "-")
            error("Unknown option $arg; use --help for available options.")
        else
            push!(positional, arg)
        end
        i += 1
    end
    length(positional) == 1 || error("Supply one .toml input file; use named options such as --output DIR.")
    path = first(positional)
    if check_only
        c = read_input(path)
        model = build_model(c)
        grid = build_grid(c, model)
        KernelWorkspace(model, grid, c) # Check the supported coupling limit without solving.
        capacity = frequency_capacity(c, get(options, :max_frequencies, nothing))
        image_factor = get(options, :image_factor, DEFAULT_IMAGE_FACTOR)
        0 < image_factor < 1 || error("image_factor must be between zero and one.")
        println("Input/model/grid checked: ", abspath(path))
        println("Source: ", SOURCE_KINDS[c.istype+1], "; surface: ", SURFACE_CONDITIONS[c.isurfcon+1])
        println("Source meaning / 物理含义: ", SOURCE_DESCRIPTIONS[c.istype+1])
        println("Screen [m]: ", (c.zstop, c.zsbtm), "; radius [m]: ", c.sradius, "; history nodes: ", length(c.ts0))
        println("Receiver radii [m]: ", c.r, "; depths [m]: ", c.zr)
        println("Window [s]: ", c.twindow, "; time samples: ", c.nt, "; accuracy: ", c.accuracy)
        println("Homogeneous layers: ", length(model.h), "; depth levels: ", length(grid.z))
        println("Frequency capacity: ", capacity, "; image factor: ", image_factor)
        println("No field calculation was performed; convergence has not been assessed.")
        return
    end
    result = run_simulation(path; options...)
    println("Numerical criteria: fields=", result.status.field_converged,
        "; equivalent source=", something(result.status.equivalent_converged, "not requested"),
        "; singular samples=", result.status.singular_samples)
    for path in result.paths
        println("Wrote ", path)
    end
end

"""Run the Coulomb-stress command line, shared by the Julia script and executable."""
function coulomb_cli(args)
    if "--version" in args
        println("MaoyePoro ",PROGRAM_VERSION," (mycfs)")
        return
    end
    if isempty(args) || args == ["--help"] || args == ["-h"]
        println("Usage: mycfs stress.toml [--help | --version]")
        println("   or: julia --project=. coulomb.jl stress.toml")
        return
    end
    length(args) == 1 || error("Supply exactly one Coulomb stress TOML file; use --help for usage.")
    run_coulomb(only(args))
    nothing
end

"""PackageCompiler simulation entry point; runtime options follow --julia-args."""
function julia_main()::Cint
    try
        BLAS.set_num_threads(1)
        cli_main(ARGS)
        return 0
    catch exception
        showerror(stderr, exception, catch_backtrace())
        println(stderr)
        return 1
    end
end

"""PackageCompiler Coulomb-stress entry point."""
function coulomb_main()::Cint
    try
        BLAS.set_num_threads(1)
        coulomb_cli(ARGS)
        return 0
    catch exception
        showerror(stderr, exception, catch_backtrace())
        println(stderr)
        return 1
    end
end
