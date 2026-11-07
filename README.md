# MaoyePoro

计算井段注水、抽水和竖向加载引起的孔压扩散、地层变形与渗流。

MaoyePoro V1.0 采用轴对称、小变形、线性、准静态 Biot 孔隙弹性模型。输入源几何与时间历史、材料参数和接收位置，输出位移、应变、倾斜、孔压与 Darcy 通量，并提供等效流体源和库仑应力后处理。

[下载可执行程序](https://github.com/fumeng6/maoyeporo/releases/latest) · [用户手册](docs/用户手册.md) · [原理](docs/MaoyePoro原理.md) · [输入输出格式](docs/输入输出格式.md)

## 研究用途

| 问题 | 输入与计算 |
| --- | --- |
| 注水、抽水与停井 | 给定井段流量历史，计算孔压、变形与渗流的时间响应 |
| 地表或地下加载 | 给定竖向总力，计算耦合孔压与变形 |
| 给定排水体积变化 | 计算该体积变化引起的孔压和附加变形 |
| 给定源区平均孔压 | 输入平均孔压乘源体积，计算场量和等效注入率 |
| 断层应力分析 | 从应变、孔压和接收层材料计算库仑应力变化 |

介质可以均匀、水平分层或沿深度渐变。顶部可选全空间、自由排水面和封闭自由面。源为有限圆柱或竖向力盘，径向采用抛物分布，整个源区位于同一实际均质子层内。

## 安装与运行

从 Releases 下载对应平台的完整程序包。

| 平台 | 安装命令 | 直接运行 |
| --- | --- | --- |
| Windows 10/11 x64 | 双击 `install.cmd`，重新打开终端 | `.\bin\myporo.exe --help` |
| Linux x86_64，glibc 2.34+ | `sh install.sh`，将 `~/.local/bin` 加入 PATH | `./bin/myporo --help` |

运行时和计算依赖随包提供，保留完整解压目录。终端命令为 `myporo` 与 `mycfs`。安装和源码构建步骤见[用户手册](docs/用户手册.md#installation)。

从恒定注水输入开始。

```sh
myporo cases/constant_injection.toml --check
myporo cases/constant_injection.toml --output results/constant_injection --format both --julia-args -t4
```

该问题在均匀全空间中设置半径 100 m、深度 900–1100 m 的注入井段，零时刻注入率阶跃至 0.01 m³/s。计算六个位置、65 个时刻和十个场分量。`--check` 读取并检查输入，求解命令输出完整归档和文本结果。

`cases/` 提供三份合成输入。

| 文件 | 物理过程 |
| --- | --- |
| [constant_injection.toml](cases/constant_injection.toml) | 有限井段恒定注水 |
| [injection_cycle.toml](cases/injection_cycle.toml) | 延迟注水与停井后的持续扩散 |
| [layered_injection.toml](cases/layered_injection.toml) | 两层介质中的井段注水 |

## 使用自己的输入

生成待填写配置。

```sh
myporo --template my_case.toml
```

模板包含物理含义和 SI 单位。填写源区、时间历史、材料、顶部条件和接收位置，随后运行 `myporo my_case.toml --check`。长历史可由 `source.history_file` 引用两列 `time_s,value` 的 CSV，相对路径以 TOML 所在目录为准。

长度为 m，时间为 s，模量与孔压为 Pa，流量为 m³/s。深度和垂向位移向下为正，正流量表示注水，负流量表示抽水。计算采用两个接收轴的全部“半径 × 深度”组合。材料突变用相同深度的上下两侧节点表示。

## 读取结果

默认只保存 `result.jld2`，包含全部十场、坐标、单位、输入、实际模型、求解设置和简短计算状态。`--format text` 导出文本及简短的 `simulation_info.txt`，`--format both` 同时保存两种格式。

需要检查详细数值误差时，使用 `myporo case.toml --diagnostics` 或在 TOML 的 `[output]` 中设置 `diagnostics = true`，单独生成 `diagnostics.txt` 并显示详细进度。默认终端只显示简短计算状态、保存路径和必要警告；`--quiet` 可关闭诊断模式下的进度输出。

| 字段 | 物理量 | 单位 |
| --- | --- | --- |
| `uz`、`ur` | 垂向与径向位移 | m |
| `ezz`、`err`、`ett`、`ezr` | 正应变和张量剪应变 | 1 |
| `tlt` | 径向倾斜 | rad |
| `pp` | 孔压变化 | Pa |
| `dvz`、`dvr` | 垂向与径向 Darcy 通量 | m/s |

计算完成后查看 JLD2 读回结果的 `status` 或文本输出的 `simulation_info.txt`；未满足数值判据时，开启诊断查看 `diagnostics.txt`。开关只控制详细信息的保存与展示，收敛检查始终执行。对关注的场量作采样加密比较，弱扩散信号同时检查绝对误差和时间周期延拓影响。诊断定义见[输入输出格式](docs/输入输出格式.md#数值诊断)。

## Julia 接口

推荐 Julia 1.13.2。在源码根目录准备依赖。

```sh
git clone https://github.com/fumeng6/maoyeporo.git
cd maoyeporo
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --threads=4 --project=.
```

```julia
using MaoyePoro, LinearAlgebra
BLAS.set_num_threads(1)
result = run_simulation("cases/constant_injection.toml";
    output_dir="results/constant_injection", format=:both)
q = sample_field(result, :pp;
    r=300.0, z=500.0, tspan=(20000.0,40000.0))
q.time
q.data
```

`solve` 返回内存结果，`run_simulation` 求解并保存。`load_result` 读取完整归档，`read_field` 选择读取一个场，`sample_field` 按坐标查询，`strain_tensor` 返回北、东、下坐标中的应变张量。`coulomb_from_strain` 使用接收层 Biot 系数构造总应力，并计算断层上的库仑应力变化。

## 项目文件

`src/` 为计算源码，`cases/` 为三份合成输入，`build/` 为构建与安装脚本。用户文档集中在 `docs/` 的三份文件中。

- [用户手册](docs/用户手册.md)说明安装、建模、计算、查询和排错。
- [MaoyePoro 原理](docs/MaoyePoro原理.md)说明物理方程、源与边界、数值处理和方法来源。
- [输入输出格式](docs/输入输出格式.md)说明配置、归档、文本交换和诊断字段。

## 来源、许可与引用

分层求解以 Rongjiang Wang 的 [POEL_2024](https://github.com/RongjiangWang/POEL_2024) 为基础，沿用部分保留原作者署名与版权声明。代码与方法对应关系见[原理文档](docs/MaoyePoro原理.md#方法来源署名与许可)。本项目采用 [GPL-3.0](LICENSE)。科研使用时记录程序版本、源码提交和实际输入，并按采用的模型与方法引用下列论文。

- Wang, R. & Kümpel, H.-J. (2003). *Poroelasticity: Efficient modeling of strongly coupled, slow deformation processes in a multilayered half-space*. Geophysics, 68, 705–717.
- Rudnicki, J. W. (1986). *Fluid mass sources and point forces in linear elastic diffusive solids*. Mechanics of Materials, 5, 383–393.

问题反馈可附最小合成输入、实际报错及有关诊断。[Issues](https://github.com/fumeng6/maoyeporo/issues)
