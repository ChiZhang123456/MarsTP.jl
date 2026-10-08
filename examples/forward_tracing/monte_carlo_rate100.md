# O₂⁺：每源面单元 100 次抽样，rate weight 与探头三维 PSD

此运行按 [detector_3d_psd.md](detector_3d_psd.md) 的驻留时间估计器，并直接调用 [源采样函数](../../src/tracing/monte_carlo_weight.jl) 的 `sample_maxwellian_source`。它与上一轮 `n|U|` 面产生率模型不同，不能把两轮差异仅解释为增加粒子数后的统计变化。

## 配置

- O₂⁺，静态 MHD 总 E、B，非相对论，无体积产生、碰撞、复合、重力或反馈。
- 释放面与吸收内边界：离地 500 km。Rm=3390 km，外边界为 4 Rm。
- 71×143=10,153 个原生球面角度单元，每个单元 **100 次 proposal draws**，共 1,015,300 次。
- 单元面积为 `r² (cos(theta0)-cos(theta1)) (phi1-phi0)`。以各单元面积中点作为局部源 patch，读取其 MHD n、Ui、Ti，所有 100 个初始位置在该点，遵循 `maxwellian_source.jl` 的局部 patch 方法。
- 采样完整、未截断的三维漂移 Maxwellian，Ts=4Ti。每格 RNG 为 `Xoshiro(20260906+cell_id)`，与线程调度无关。
- 正向 Boris，dt=0.1 s，最大飞行年龄 500 s；内边界、外边界及时限分别记录。
- 探头中心 (1,0,2) Rm，立方体边长 0.2 Rm，即 678 km。位置和矢量使用原生 MHD 笛卡尔轴，不冒称已核实 MSO/MSE。
- 12 线程，全部实际传播轨迹逐步保存同步位置、速度及时间，双精度，无压缩。零率样本保存初始状态和零率标记，避免无物理贡献的积分。

## 权重与单位

重要性权重：

```text
w_i = g(v_i; Ui,Ti) / gs(v_i; Ui,4Ti)            [dimensionless]
```

保留 MarsASPEN 风格的源密度权重作为独立诊断量：

```text
source_density_weight_m3 = n * w_i / sum_cell(w) [m^-3]
```

**实际输运使用**：

```text
Q_i = rate_weight_s1 = n*A*max(dot(v_i,er),0)*w_i/N  [s^-1]
N = 100, including inward draws
```

向内样本 Q=0，保留为 `zero_rate`，不重新抽样补足 100 个向外粒子，也不按正率粒子数重新归一化。Q 不做自归一化，因此每格总率不强制等于 n|Ui|A。CSV 中 `source_flux_m2_s` 仍记录 n|Ui|，仅作与旧模型对照的 MHD 诊断，**没有用于生成 Q**。

三维 PSD：

```text
occupancy[b] = sum_i Q_i * residence_time_i_in_bin     [particles]
f3d[b] = occupancy[b] / (V_detector * dvx*dvy*dvz)     [s^3 m^-6]
n_detector = sum_b f3d[b] * dvx*dvy*dvz               [m^-3]
```

不再除以总积分时间或粒子数。空间轨迹段先与探头立方体相交，随后按线性速度变化与速度 bin 边界进一步切分。每步端点速度与位置同步，不使用 Boris 半步速度分 bin。

穿面事件另外保存 `face_flux_m2_s=Q_i/A_face`，单位 m⁻² s⁻¹。按各面的进入或离开方向分别求和得到面积平均通量。

## 图与速度窗口

三维速度 bin 宽度为 10 km/s，bin 范围自动包含所有探头速度，不截断计算出的密度。

用户指定所有图的 xlim、ylim 为 **[-200,200] km/s**。这是显示窗口，范围外的轨迹、驻留时间和三维 PSD 都保留。`analysis_summary.json` 单独报告 XY、XZ 投影落在显示窗口内的密度比例。

- `probe_psd_projections.png/pdf`：分别积分掉 vz、vy，单位 s² m⁻⁵。
- `probe_psd_slices.png/pdf`：省略速度分量在中心 bin 的薄片平均，单位 s³ m⁻⁶。
- `probe_psd.npz`：完整三维 PSD、速度 edges、投影、中心薄片、各 bin 独立粒子数和有效粒子数。

灰色表示没有抽样贡献，不代表物理上严格为零。若信号完全在显示窗口外，图内明确注明。

## 输出与复现

完整运行：`outputs/mc500_rate100_full_20260906a`。预试验为 `outputs/mc500_rate100_pilot_20260906a`，只有 80 个源面单元，不代表完整球面。

- `source_cells.csv`：面积、源 n/T/Ui 及 n|Ui| 诊断。
- `particles.csv`：所有 proposal draws 的初末状态、`rate_weight_s1`、`source_density_weight_m3`、重要性比、终止原因及探头驻留时间。
- `trajectories_*.jld2`：每个 `p<ID>` 组的 `state`、`rate_weight_s1`、`source_density_weight_m3`、`cell_id`。state 在 Julia 中为 7×N，在 h5py 中为 N×7，列依次 t,x,y,z,vx,vy,vz，单位 s、m、m/s。
- `probe_residence.csv` 与 `probe_crossings.csv`：保存完整探头记录，可独立重新分 bin。
- `metadata.toml`：参数、输入 SHA256、Git、软件版本和实际 VTK 网格检查。仍直接使用文件实际径向网格，避免旧解析轴与文件相差约 353 m 的问题。
- `completion.toml`：只在全部批次成功后写出。
- `analysis_summary.json`、`qa_summary.json` 与 refinement CSV：单位、统计及数值检查。

PowerShell，在仓库根目录：

```powershell
$env:MC_GIT_COMMIT = git rev-parse HEAD
$env:MC_GIT_STATUS = (git status --short) -join "`n"
julia --startup-file=no --compiled-modules=existing --threads=12 --project=. examples/forward_tracing/monte_carlo_shell.jl outputs/new_rate100_run 1 500 0.1 100 reservoir_maxwellian_rate false
& 'C:\Users\Win\.conda\envs\mars\python.exe' examples/forward_tracing/analyze_monte_carlo.py outputs/new_rate100_run --plot-limit-kms 200
```

脚本拒绝覆盖已有目录。后六项依次为 cell stride、最大飞行年龄、dt、每格 proposal 数、源模型、是否压缩轨迹。

## 验证范围

除已有解析运动和几何测试外，新增均匀零漂移 Maxwellian 球面测试，将抽样向外率与 `n*A_total*sigma/sqrt(2*pi)` 比较，同时验证源面面积和密度权重。真实 MHD 的解析局部 Maxwellian 单向通量也在分析时与全源 Monte Carlo 粒子率比较，并记录 Monte Carlo 标准误差。

实际结果仍为最大飞行年龄 500 s 的分布贡献，不能自动称为稳态收敛解。改变源模型后，上一轮密度和粒子数不构成严格的同源统计收敛对照。独立源样本数、时间步、源面空间离散、速度 bin 和最大飞行年龄需要分别评估。

## 本轮已完成结果（2026-09-06）

- 1,015,300 次 proposal draws；506,259 个正粒子率样本传播，509,041 个零率样本保留初态。
- 正率轨迹：315,992 个返回 500 km，108,153 个到达 4 Rm，82,114 个达到 500 s 时限，无数值失败。
- 积分与写盘用时 589.5 s（12 线程），轨迹文件合计 48,255,689,492 bytes，约 48.3 GB。
- 源率为 2.23077×10²⁴ s⁻¹，解析局部 Maxwellian 向外率的全表面积分为 2.24031×10²⁴ s⁻¹。Monte Carlo 标准误差 7.23088×10²¹ s⁻¹，两者相差 -1.32 个标准误差。
- 探头命中 416 个独立源样本，834 次穿面事件，11,151 个步内驻留段。少量重复穿越按驻留时间正常累计，有效粒子数按独立源样本合并后为 119.29。
- 完整速度范围内的探头密度为 30,966.9 m⁻³，即 0.0309669 cm⁻³。
- 完整三维速度 bin 范围为各轴 ±265 km/s，宽度 10 km/s，数组 53×53×53。显示范围按要求全部为 ±200 km/s。
- XY 图显示完整密度贡献，XZ 图的非零贡献全部在其显示窗口外。范围外数据未丢弃，也未用显示窗口归一化总密度。
- 达到时限的源率权重占 18.13%，尚未进行最大飞行年龄收敛。
- 抽查全部 992 个批次的首尾完整轨迹，共 1,984 条，权重、时间、初末位置速度、有限值和边界范围通过；全部 1,015,300 个粒子组计数通过。
- 200 个最大探头密度贡献者加 31 个正率空间对照，共 231 个粒子，以 0.1、0.05、0.025 s 重算。子集探头密度的相邻变化约 -0.00107% 和 -0.00218%，终止类型未改变。这不是完整集合 PSD 收敛测试。
- 整体最大绝对做功闭合残差 0.173 eV。320 条轨迹的相对残差超过 1%（分母采用 max(|ΔK|,|W|,1 eV)）。探头贡献粒子的最大相对残差为 2.46×10⁻⁷。原始 0.1 s 结果全部保留，未选择性替换轨迹。
- 新模拟测试 226 项、Python 分析测试 4 项通过；另完成全部逐粒子源率重建、每格源密度份额归一化和三维 PSD 积分密度一致性检查。
