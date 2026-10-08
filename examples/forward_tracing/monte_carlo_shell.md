# 500 km O₂⁺ Monte Carlo 面源与有限体积探头

## 模型与默认配置

该示例模拟稳恒面源释放的独立轨迹。每条轨迹的时间是释放后的飞行年龄，不是一个瞬时粒子云的统一观测时刻。

- O₂⁺，质量和正电荷取自 TestParticle，非相对论洛伦兹运动。
- 火星半径 3390 km，释放面为离地 500 km，吸收内边界为离地 200 km，外边界为火心距离 4 Rm。
- 静态 MHD 总电场与磁场，不添加体积产生、碰撞、复合、重力或粒子反馈。
- 固定 Boris 步长 0.1 s，最大飞行年龄 500 s。达到时限单独记为 `time_limit`，不声称稳态收敛。
- 探头中心 (1, 0, 2) Rm，立方体**边长** 0.2 Rm，即 678 km。采用文件原生笛卡尔轴，不从未确认的元数据推断 MSO/MSE 标签。
- 使用原生角度节点之间的 71 × 143 个球面单元，每格 10 个粒子，总共 101,530 个。经度末端不重复，极区使用有限面积扇区。
- 面内位置按面积均匀随机抽样，余纬的 cos 和经度分别均匀抽样。各格使用该单元面积中点的 MHD 矩作为分片常数源。
- 每粒子随机流为 `Xoshiro(seed + particle_id)`，不依赖线程调度。默认 seed=20260906。

VTK 的实际坐标以 Rm 为尺度。新示例将 `norm(Points)` 乘项目 Rm，直接读取实际径向轴，并检查全部网格点是否与球坐标张量网格一致。旧加载器按 200 km 内端点重建的径向轴与原文件相差约 353 m，此示例在自身运行中修正插值轴，不修改原文件或现有回溯代码。E/B 的原始笛卡尔分量经 MarsTP 转为球分量供现有插值器使用，返回积分器的是笛卡尔矢量；O₂⁺ 的 Ui 按笛卡尔标量分量分别插值。

## 源权重

按用户指定沿用 `F=n|Ui|`，单位 m⁻² s⁻¹。该量是面产生率处方，不是球面法向净通量。

球面单元面积：

```text
A = (Rm + 500 km)^2 [cos(theta0) - cos(theta1)] (phi1 - phi0)
```

参考 MarsASPEN 的 `monte_carlo_weight.jl`，在每格使用自归一化重要性抽样。物理速度 PDF 是截取 `v·er>0` 后归一化的漂移 Maxwellian `g+`。抽样 PDF `gs+` 的温度为 `Ts = sampling_temperature_factor * Ti`，默认因子 4。两者漂移 Ui 相同。

```text
logw_i = log(g+(vi)/gs+(vi))
Wi = F A exp(logw_i) / sum_cell exp(logw)
sum_cell Wi = F A                 [particles s^-1]
```

代码使用 log-sum-exp，避免下溢造成归一化失败。向外截取的归一化概率也包含在比值中，特别是 Ui 径向分量向内时。`Ts=Ti` 时比值为 1，各格等权。**不能在已经按物理 Maxwellian 抽样时再次额外乘一个 Maxwellian。** 本权重是注入率，不是密度，也不是一次瞬时释放的粒子数。若另行定义脉冲释放时长 Δt_emit，脉冲粒子数权重才是 Wi Δt_emit。

每格仅 10 个抽样，自归一化重要性抽样有有限样本偏差。输出有效粒子数作为集中度诊断，不将其视为包含空间分层效应的严格误差估计。增大 Ts 可能改善尾部覆盖，同时降低有效粒子数。

## 探头与 PSD

每个积分步用线段与立方体的几何交集计算准确的步内进入和离开时间，避免越过探头但两个端点都在外面时漏计。这里的“准确”指 Boris 位置线段模型，真实曲线仍有时间步误差。记录的速度同步到对应位置时刻，不保存为半时间步速度。

对于每个速度 bin：

```text
N_bin = sum_i Wi * residence_time_i_in_bin          [particles]
f3d = N_bin / (V_probe * dvx * dvy * dvz)           [s^3 m^-6]
n_probe = integral f3d d^3v                        [m^-3]
```

这里不再除以 500 s。它是持续注入下、飞行年龄不超过 500 s 的体积平均分布贡献。保存每段驻留的两端速度，分析脚本假设段内速度线性变化，并按速度 bin 边界切分驻留时间。

默认速度 bin 宽度为 10 km/s，范围按全部探头端点速度自动扩展，不静默丢掉越界粒子。内部所有 bin 宽度先转为 m/s。图中速度轴显示 km/s，但 PSD 数值保留 SI 单位。

- `fxy = integral f3d dvz` 与 `fxz = integral f3d dvy`，单位 s² m⁻⁵，输出主要两面板图。
- 另输出省略速度分量中心 bin 的 slab 平均值图，单位 s³ m⁻⁶。它不是严格零宽度数学切片。
- 灰色 bin 表示无抽样贡献，不代表物理上严格为零。不平滑，不插值填洞。

每次穿面事件另外保存 `Wi / A_face`，单位 m⁻² s⁻¹。对给定面的进入事件求和得到该面的面积平均单向进入通量；离开另计。不能把六面的标量和当成三维矢量通量。

体积平均数通量矢量另外按 `sum Wi ∫v dt / V_probe` 计算，单位 m⁻² s⁻¹。

## 文件

- `metadata.toml`：全部配置、SI 单位、源 SHA256、软件版本、Git 信息与网格检查。
- `source_cells.csv`：每格面积、MHD n/T/Ui 与 F。
- `particles.csv`：所有粒子的 ID、格子、注入率权重、初末位置和速度、log importance、终止原因、驻留时间、做功闭合诊断。
- `trajectories_*.jld2`：分批压缩保存全部 0.1 s 状态及精确步内终止点。每组 `p<ID>/state` 在 Julia 中为 7×N，行依次为 t,x,y,z,vx,vy,vz。每组另存 `weight_s1` 和 `cell_id`。位置 m，速度 m/s，时间 s。
- `probe_residence.csv`：全部步内驻留段，含权重、两端位置与速度，供重新建立任意三维速度 bin。
- `probe_crossings.csv`：穿面时刻、权重、面 flux、速度与位置。面编号 ±1=x±、±2=y±、±3=z±；direction=1 进入、-1 离开。
- `probe_psd.npz`：3D PSD、2D 积分图、中心 slab 图、各速度 bin 的独立粒子数和有效粒子数。
- `analysis_summary.json`：密度、数通量、统计量、时限权重比例与能量诊断。
- `completion.toml`：仅在全批次成功后生成。没有此文件的运行不应解释为完整结果。
- 运行目录内保存 Project、Manifest 和模拟脚本快照。拒绝覆盖已有运行目录。

Python `h5py` 可读取 JLD2 中这些纯数值数组；`state` 在 h5py 中通常呈 N×7：

```python
import h5py
with h5py.File('trajectories_00001.jld2') as f:
    state = f['p1/state'][:]  # N x 7, verify shape before use
    weight_s1 = f['p1/weight_s1'][()]
```

## 运行

在仓库根目录，PowerShell：

```powershell
$env:MC_GIT_COMMIT = git rev-parse HEAD
$env:MC_GIT_STATUS = (git status --short) -join "`n"
julia --startup-file=no --compiled-modules=existing --threads=4 --project=. examples/forward_tracing/test_monte_carlo.jl
julia --startup-file=no --compiled-modules=existing --threads=4 --project=. examples/forward_tracing/monte_carlo_shell.jl outputs/mc500_new 1 500 0.1
& 'C:\Users\Win\.conda\envs\mars\python.exe' examples/forward_tracing/analyze_monte_carlo.py outputs/mc500_new --dv-kms 10
```

命令行最后三个数是格子抽样 stride、最大飞行年龄、步长。stride=1 才是完整球面；stride>1 仅供小样本测试，不将稀疏格子的结果重标定成完整球面结果。

更改探头、抽样温度或每格粒子数可通过 Julia 的 `Config` 设置。统计收敛、步长收敛与最大飞行年龄收敛是不同问题，能量闭合本身不能证明 PSD 收敛。

## 2026-09-06 首轮完整结果

以下为旧 500 km 吸收边界的历史运行，不能直接代表当前 200 km 配置。完整运行目录：`outputs/mc500_full_20260906c`。此前无 completion 文件的目录为保留的输入检查或中断试跑，不用于结果解释。

- 全部 101,530 个粒子完成，63,243 个返回内边界，21,763 个到达外边界，16,524 个达到 500 s 时限。
- 源总注入率 2.15423×10²⁴ s⁻¹，达到时限的权重占 15.66%。
- 探头命中 100 个独立粒子，200 次穿面事件，2,633 个步内驻留段。有效粒子数 19.76。
- 有限飞行年龄下的体积平均密度 0.0386717 cm⁻³。
- 默认 10 km/s bin，完整三维速度范围为各分量 ±285 km/s。另有显示全部非零 bin 的放大图 `probe_psd_projections_zoom.png`。
- 全部轨迹压缩后约 8.41 GB；逐个检查了粒子组总数，每批首尾完整轨迹的时刻、权重、初末状态、有限值与边界范围均检查通过，共 794 条。
- Julia 项目现有 227 项测试、新增模拟 217 项测试和 Python 分析 4 项测试通过。
- 全部 100 个探头贡献者加 32 个分散选取的对照粒子，分别以 0.1、0.05、0.025 s 重算。相邻步长的探头密度贡献变化为 0.000477% 与 0.000246%，该子集的终止类型没有改变。此检查不能排除原先未命中粒子在较细步长下命中，不能替代完整集合的 PSD 收敛。
- 全部粒子中有 71 个做功闭合误差超过 1%（分母下限 1 eV），最大绝对残差 0.133 eV。最差 16 个案例重算后，最大绝对残差在 0.05 s 时为 0.0452 eV，0.025 s 时为 0.00845 eV。全部原始 0.1 s 结果保留，没有选择性替换。探头贡献粒子的最大相对闭合残差为 2.53×10⁻⁷。

上述结果是首轮 Monte Carlo 估计，尚未验证粒子数、抽样温度、探头尺寸或飞行年龄收敛。详细数值见运行目录中的 `analysis_summary.json`、`qa_summary.json` 与两份 refinement CSV。
