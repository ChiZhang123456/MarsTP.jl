# O₂⁺ Monte Carlo forward tracing, 500 km source and cubic detector

此示例从 500 km 电离层边界释放 O₂⁺，保存每个粒子的权重、位置、速度，并由探头内的驻留时间建立三维速度分布。理论推导见 [monte_carlo.md](monte_carlo.md)。本目录提供实际运行脚本与 PNG，不包含约 48.3 GB 的原始轨迹或 MHD 输入。

## 图

![Velocity-integrated PSD](probe_psd_projections.png)

两个投影分别积分掉 vz、vy，单位为 **s² m⁻⁵**。三维速度网格各轴从 −500 到 500 km/s，宽度 **5 km/s**，共 **200³** 个 bin；绘图 xlim、ylim 均为 **[−300, 300] km/s**。使用 turbo、对数色标和 Arial，PNG 图底部不加说明文字。灰色表示没有抽样贡献，无平滑。

三维网格以 COO 稀疏格式保存非零 bin，未保存的 bin 为本次抽样估计中的零值。先计算完整三维 PSD，再对完整的被省略速度轴求和乘以 Δv = 5000 m/s，得到二维投影。显示窗口仍为 ±300 km/s，三个探头的全部信号都在此窗口内。这里不再输出零速度切片。

两个新增探头示例：[**(0,0,2) Rm**](probe_0_0_2/README.md)、[**(−1.5,0,1) Rm**](probe_m1p5_0_1/README.md)。它们使用同一套完整已保存轨迹，不重新采样或积分。

![5000 trajectories in XZ, XY and YZ](trajectories_5000.png)

从正率且实际传播时间大于 0 的粒子中，用 NumPy `default_rng(20260906)` 蓄水池抽样等概率选出 **5000** 条，不优先挑选探头命中者。颜色表示终止原因：`inner` 返回 500 km，`outer` 到达 4 Rm，`time_limit` 达到 500 s。这是轨迹示意图，线条没有按物理权重加粗或重采样。

绘图读取已保存的实际轨迹，每条最多显示 300 个均匀索引点，保留首尾点；此降采样仅用于显示，积分、探头统计和原始数据均保持 0.1 s。选择 ID 保存在 [trajectories_5000.selection.json](trajectories_5000.selection.json)。

背景使用 `py_space_zc.maven.bs_mpb` 和 `plot_mars`。XZ、XY 中的虚线、点线分别为库内 BS、MPB 经验 X–ρ 曲线，假定 +X 朝向太阳，仅作参考，不是从本次 MHD 提取的边界。YZ 不画这种二维曲线。圆盘是火星，外侧灰圆是 500 km 释放面，洋红方框是探头投影。**数据仍标为原生 MHD 笛卡尔坐标，尚未核实其 MSO/MSE 来源**，因此不能据参考曲线判断边界符合程度。

## 模型和单位

- O₂⁺，质量 5.352390155808×10⁻²⁶ kg，电荷 +e；静态 MHD 总 E、B，非相对论 Boris 积分。
- Rm = 3390 km；释放面及吸收内边界为离地 500 km，外边界为火心距离 4 Rm。
- 原生球面 71×143 = 10,153 个角度单元，各 100 次完整 Maxwellian proposal draws，总计 1,015,300 次；局部 n、Ui、Ti 在各面积中点读取。
- proposal 温度 Ts = 4Ti，各格 Julia RNG 为 `Xoshiro(20260906 + cell_id)`，不依赖线程顺序。
- dt = 0.1 s，最大飞行年龄 500 s；探头中心 (1,0,2) Rm，边长 0.2 Rm = 678 km。
- 不加入电离层外体积产生、碰撞、复合、重力或粒子反馈。

直接调用 [sample_maxwellian_source](../../../src/tracing/monte_carlo_weight.jl)，采用 reservoir 的单向穿面率权重：

```text
w_i = g(v_i; Ui,Ti) / gs(v_i; Ui,4Ti)                [dimensionless]
Q_i = rate_weight_s1 = n A max(v_i · er, 0) w_i / N   [s^-1]
N = 100, including inward draws
source_density_weight_m3 = n w_i / sum_cell(w)        [m^-3, diagnostic only]
```

向内样本 Q=0，保存初态但不传播，不重新抽样补足向外样本。Q 不做自归一化。CSV 的 `source_flux_m2_s=n|Ui|` 仅是旧模型对照诊断，**不用于此运行的 Q**；因此本运行不是固定 n|Ui|A 总率模型。

```text
f3d[b] = sum_i Q_i * dwell_time_i_in_bin / (V * dvx*dvy*dvz) [s^3 m^-6]
n = sum_b f3d[b] * dvx*dvy*dvz                             [m^-3]
fxy = integral f3d dvz; fxz = integral f3d dvy              [s^2 m^-5]
face_flux_m2_s = Q_i / A_face                             [m^-2 s^-1]
```

内部位置、速度分别为 m、m/s，速度 bin 体积用 SI 单位；图轴才转为 km/s。驻留时间估计器不再除以总运行时间或粒子数。每步先求轨迹段与立方体交集，再按线性插值速度穿越 bin 的位置切分驻留时间。位置与速度使用同步端点。

注意：本示例 `--dv-kms 5` 表示 bin **宽度**；包内 `forward_psd` 的 `vgrid` 表示每轴 bin **数量**，等效设置为 `vlim=(-500,500), vgrid=200, velocity_unit=:km_s`。

## 文件与复现

在仓库根目录运行，Julia 使用本仓库 Project.toml / Manifest.toml。原运行使用 Julia 1.12.6、MarsTP 提交 `58e078fa4d853a0307df0e9c87b9b44528835dc7`。Python 需 NumPy、Matplotlib、h5py；轨迹背景还需要用户本地 `py_space_zc`。Windows 已验证环境为 `C:\Users\Win\.conda\envs\mars\python.exe`。

输入：`data/mars_fields_spherical_from_dat.vts`，SHA256 为 `fa92bd82fe16975ad0d50f4e40ace344e9d41389d6024976d423e6126bc7a5c8`。输入需自行提供，脚本核对文件实际径向网格，不用存在约 353 m 偏差的旧解析径向轴代替它。

新探头提取脚本另需已安装的 Numba。三维稀疏存储并不改变 bin 定义或归一化，例如直接读取发布的 NPZ 重建 XY 投影：

```python
import numpy as np
with np.load('probe_psd_sparse.npz') as a:
    ix, iy, iz = a['indices_xyz'].T
    fxy = np.zeros(tuple(a['shape_xyz'][:2]))
    np.add.at(fxy, (ix, iy), a['f3d_s3_m6'] * float(a['dv_ms']))
    np.testing.assert_allclose(fxy, a['fxy_s2_m5'])
```

```powershell
$example = 'examples/forward_tracing/monte_carlo_forward_tracing'
$run = 'outputs/new_rate100_run'
$py = 'C:\Users\Win\.conda\envs\mars\python.exe'
$env:MC_GIT_COMMIT = git rev-parse HEAD
$env:MC_GIT_STATUS = (git status --short) -join "`n"
julia --startup-file=no --compiled-modules=existing --threads=12 --project=. "$example/monte_carlo_shell.jl" $run 1 500 0.1 100 reservoir_maxwellian_rate false
& $py "$example/analyze_monte_carlo.py" $run --dv-kms 5 --vmax-kms 500 --plot-limit-kms 300 --output-dir "$run/analysis_5kms"
& $py "$example/plot_trajectories.py" $run --output "$run/trajectories_5000.png" --count 5000 --seed 20260906
```

模拟输出目录必须不存在。模拟位置参数依次是输出目录、cell stride、最大飞行年龄、dt、每格抽样数、源模型、JLD2 压缩开关。先用较大 stride、较短时间进行小样本试运行。若 SpacePy 默认配置目录不可写，可将进程环境变量 `SPACEPY` 指向可写目录，它会在该目录下创建 `.spacepy`。若 `py_space_zc` 导入时停在 Numba 缓存初始化，将 `NUMBA_CACHE_DIR` 指向可写缓存目录即可；无需修改已安装的库。

| 文件 | 内容 |
| --- | --- |
| [monte_carlo_shell.jl](monte_carlo_shell.jl) | 源采样、Boris 传播、边界、探头相交与流式保存 |
| [analyze_monte_carlo.py](analyze_monte_carlo.py) | 绘图助手及稠密测试参考，主入口转入稀疏分析 |
| [analyze_probe.py](analyze_probe.py) | 完整三维稀疏 PSD、全速度积分投影、密度一致性验证 |
| [reprobe_saved.py](reprobe_saved.py) | 扫描全部轨迹，一次提取三个探头的驻留段 |
| [synchronize_reprobe.jl](synchronize_reprobe.jl) | 根据原 Boris drift 速度同步探头交点速度 |
| [plot_trajectories.py](plot_trajectories.py) | 固定种子 5000 条轨迹三平面 PNG |
| [qa_monte_carlo.py](qa_monte_carlo.py) | 批次轨迹、初末状态、粒子组计数检查 |
| [refine_monte_carlo.jl](refine_monte_carlo.jl) | 选定子集以 0.1、0.05、0.025 s 重算 |
| [test_monte_carlo.jl](test_monte_carlo.jl) | 几何、采样归一化、解析传播测试 |
| [test_monte_carlo_analysis.py](test_monte_carlo_analysis.py) | 速度 bin 穿越、密度单位及相关样本检查 |

原始输出包含 `particles.csv`、`source_cells.csv`、`probe_residence.csv`、`probe_crossings.csv`、`metadata.toml`、`completion.toml` 和 `trajectories_*.jld2`。JLD2 各 `p<ID>/state` 在 Julia 中为 7×N，h5py 中为 N×7，列为 t,x,y,z,vx,vy,vz，单位 s、m、m/s。权重也保存在各组中。分析另存 `probe_psd_sparse.npz`、`analysis_summary.json` 和粒子贡献 CSV；稀疏文件保存零基 `indices_xyz`、对应 `f3d_s3_m6`、完整速度边界、200³ 形状及两张二维投影。无需构造稠密 200³ 数组。目录内的小型 NPZ 可用于读取已发布的 PSD，原始大型轨迹不上传。

```powershell
julia --startup-file=no --compiled-modules=existing --project=. "$example/test_monte_carlo.jl"
& $py "$example/test_monte_carlo_analysis.py"
& $py "$example/test_reprobe.py"
# 可选：对已有完整运行进行子集时间步细化，再核验轨迹文件。
julia --startup-file=no --compiled-modules=existing --project=. "$example/refine_monte_carlo.jl" $run
& $py "$example/qa_monte_carlo.py" $run
```

## 本次结果与验证范围

图来自 `mc500_rate100_full_20260906a`，未为本次改图重新积分。正率粒子 506,259 个，零率 509,041 个；正率终止类型为 inner 315,992、outer 108,153、time_limit 82,114，无数值失败。12 线程积分及写盘用时 589.5 s，原始轨迹约 48.3 GB。

探头有 416 个独立粒子、834 次穿面事件、11,151 段驻留记录，驻留权重有效样本数 119.29。三维 PSD 积分密度 **30,966.9 m⁻³ = 0.0309669 cm⁻³**。三维网格 ±500 km/s，图示 ±300 km/s 包含本次全部探头速度贡献。源率 2.23077×10²⁴ s⁻¹，与解析局部 Maxwellian 全表面积分相差 −1.32 个 Monte Carlo 标准误差。

整理后通过 226 项 Julia 测试、4 项 Python 测试，以及完整源率逐粒子重建、每格密度权重和 PSD 积分密度检查。此前抽查 992 个批次首尾共 1984 条完整轨迹，并核验全部粒子组数。231 条子集的相邻时间步细化使探头密度改变约 −0.00107%、−0.00218%，终止类型不变。

这是有限样本、有限探头体积、最大飞行年龄 500 s 的结果，尚未证明稳态或完整 PSD 收敛。达到时限的源率占 18.13%。整体最大做功闭合残差 0.173 eV，320 条轨迹相对残差超过 1%（分母以 1 eV 为下限）；探头贡献粒子的最大相对残差为 2.46×10⁻⁷。原始结果均保留。

## 已保存轨迹上的三个探头

以下命令只读取原始运行输出及同一 MHD 场，不重新传播粒子。输出目录必须不存在：

```powershell
$reprobe = 'outputs/new_three_probes'
& $py "$example/reprobe_saved.py" $run $reprobe
julia --startup-file=no --compiled-modules=existing --project=. "$example/synchronize_reprobe.jl" $reprobe
foreach ($probeName in @('probe_1_0_2','probe_0_0_2','probe_m1p5_0_1')) {
    & $py "$example/analyze_probe.py" "$reprobe/$probeName" --output-dir "$reprobe/$probeName/analysis"
}
```

提取器检查全部 1,015,300 个粒子组。交点的位置、时间由保存的逐步轨迹精确裁剪得到；drift 速度由该步位移除以实际步时长重建，再调用与原追踪相同的 `TP.update_velocity` 计算交点同步速度，并验证 MHD 输入 SHA256。

原探头回放的 11,151 段驻留记录与原记录相比：时间、位置、权重相同，密度相对差为 0，最大速度差 2.45×10⁻⁸ m/s。新增 4 项检查覆盖立方体交会、平行漏过、稀疏与稠密 PSD 等价及速度越界拒绝。

| 探头位置 (Rm) | 密度 (cm⁻³) | 独立命中粒子 | 驻留权重有效样本数 | 驻留段 |
| --- | ---: | ---: | ---: | ---: |
| (1,0,2) | 0.0309669 | 416 | 119.29 | 11,151 |
| (0,0,2) | 0.001002168 | 115 | 9.71 | 3,328 |
| (−1.5,0,1) | 0.1177217 | 2,011 | 353.71 | 124,374 |

三个探头均验证 `sum(f3d)*dv³ = sum(fxy)*dv² = sum(fxz)*dv² = sum(Q*tau)/V`，所有速度积分覆盖完整 ±500 km/s。不同图采用各自的对数色标，颜色不能直接跨图比较。特别是 (0,0,2) 的有效样本数较低，5 km/s 图上的细结构仍受 Monte Carlo 噪声影响。
