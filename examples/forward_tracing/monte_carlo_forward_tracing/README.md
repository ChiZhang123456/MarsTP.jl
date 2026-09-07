# O₂⁺：Monte Carlo 采样、forward tracing 与探头 PSD

本文合并原 `monte_carlo.md` 的理论推导与运行说明，按源采样、轨迹传播、探头统计的顺序介绍实际实现。本目录只提交代码、说明与 PNG，运行生成的 JSON、NPZ 和大型轨迹文件保存在本地。

## 1. Monte Carlo 源采样

Monte Carlo 用有限的计算粒子表示连续的速度分布。为增加高速尾部的样本，本例从较热的 proposal 分布抽样，再用重要性权重恢复目标分布。计算粒子并不代表相同的物理密度或粒子率。

### 1.1 二维教学示例

![Maxwellian and Monte Carlo density and rate weights](monte_carlo_sampling.png)

[绘图代码](plot_monte_carlo_sampling.py) 使用 O₂⁺，密度 5 cm⁻³（5×10⁶ m⁻³），bulk velocity 为 (−10,0) km/s，温度 10 eV，固定 vz=0，抽取 100,000 个二维速度样本，随机种子为 20260907。proposal 温度为 40 eV。为定义单位为 s⁻¹ 的粒子率，额外指定源面面积 **1 m²**、单位法向 **−X**。

这里是归一化的二维教学模型，不是把三维分布在 vz=0 处直接截取。其解析函数也等于各向同性三维 Maxwellian 对 vz 积分后的二维边缘分布；实际电离层源仍抽样三个速度分量。

三个 panel 均使用 turbo：左侧是解析二维分布，单位 s² m⁻⁵；中间是每个样本的 `density_weight`，单位 m⁻³；右侧是每个样本的 `flux_weight`，单位 s⁻¹。后两幅是按单粒子权重着色的散点图，不是速度 bin 内的 PSD。灰色点表示穿面率为零。

### 1.2 Maxwellian 与重要性抽样

先定义温度参数 T（eV，表示 kBT 的能量值）、换算常数 κ=1.602176634×10⁻¹⁹ J/eV、热能 εT（J）、粒子质量 m=5.352390155808×10⁻²⁶ kg，以及单分量热速度标准差 σ（m/s）：

$$
\epsilon_T=\kappa T,\qquad \sigma=\sqrt{\epsilon_T/m}.
$$

本例 σ=5.471184 km/s。令 d 为无量纲速度维数，教学例 d=2、实际源 d=3；速度 v、漂移速度 U 均为 m/s，g_d 为归一化概率密度，单位 (s/m)^d；n 为 m⁻³，f_d 的单位为 m⁻³(s/m)^d：

$$
g_d(\mathbf v)=\frac{\exp[-|\mathbf v-\mathbf U|^2/(2\sigma^2)]}{(2\pi\sigma^2)^{d/2}},\qquad f_d=n g_d.
$$

令 c=4 为无量纲温度倍率，T_s 为 proposal 温度（eV），σ_s 为 proposal 标准差（m/s）：

$$
T_s=cT,\qquad \sigma_s^2=c\sigma^2.
$$

粒子速度 v_i 从 g_s 抽样。g_d 与 g_s 的单位相同，重要性权重 w_i 无量纲；v_i、U、σ 均用 m/s，c、d 无量纲：

$$
w_i=\frac{g_d(\mathbf v_i)}{g_s(\mathbf v_i)}
=c^{d/2}\exp\left[-\frac{|\mathbf v_i-\mathbf U|^2}{2\sigma^2}\left(1-\frac1c\right)\right].
$$

因此二维示例的系数是 c，实际三维实现是 c^(3/2)。实现见 [sample_maxwellian_source、maxwellian_importance_weight_3d](../../../src/tracing/monte_carlo_weight.jl)。

### 1.3 两种单粒子权重

令 N 为无量纲抽样次数，n 为 m⁻³，w_i 无量纲，密度权重 W_i^(n) 为 m⁻³：

$$
W_i^{(n)}=n\frac{w_i}{\sum_{j=1}^{N}w_j},\qquad \sum_i W_i^{(n)}=n.
$$

这是 `particle_density_weight` 使用的自归一化密度权重，仅作为源密度表示和诊断。它不是 forward tracing 所需的源率权重。

令 A 为源面面积（m²），单位法向 n̂ 无量纲，v_i、v_n,i 为 m/s；n 为 m⁻³，N、w_i 无量纲，Q_i 为粒子率（s⁻¹，粒子数视为计数）：

$$
v_{n,i}=\mathbf v_i\cdot\hat{\mathbf n},\qquad
Q_i=\frac{nA}{N}\max(v_{n,i},0)w_i.
$$

Q_i 对应库中的 `rate_weights_s`，也是教学图的 `flux_weight`。严格说它是源面总粒子率权重，而单位面积通量的单位应为 m⁻² s⁻¹。完整 proposal 中朝内的粒子保留零率，不拒绝后重抽；分母 N 包含全部抽样次数。率权重不使用随机的密度自归一化分母。

用于校验总率时，令 U_n=U·n̂（m/s）、σ（m/s）、n（m⁻³）、A（m²）；φ 和 Φ 分别为无量纲标准正态密度与累积分布，Q_exact 为 s⁻¹：

$$
Q_{\rm exact}=nA\left[\sigma\phi(U_n/\sigma)+U_n\Phi(U_n/\sigma)\right].
$$

教学例的密度权重总和恰为 5×10⁶ m⁻³；粒子率总和为 5.06360×10¹⁰ s⁻¹，解析值为 5.03641×10¹⁰ s⁻¹，相差约 1.33 个 Monte Carlo 标准误差。

### 1.4 实际 500 km 电离层源

实际计算在 500 km 高度从 MHD 读取每个面积单元的局部 n、U_i、T_i，不使用上述固定教学参数。源面也是吸收内边界，球面外法向为径向。令 r_s 为源面火心半径（m）、θ 为余纬（rad）、ϕ 为经度（rad）、A_cell 为 m²：

$$
A_{\rm cell}=r_s^2(\cos\theta_0-\cos\theta_1)(\phi_1-\phi_0).
$$

71×143=10,153 个角度单元，每格 100 次完整三维抽样，共 1,015,300 个样本；每格使用本地密度、面积和 N=100 计算 Q_i。proposal 温度为局部温度的 4 倍，随机数为 `Xoshiro(20260906 + cell_id)`。几何位置使用实际 VTS 径向坐标。

当前方案使用上述单向 Maxwellian 穿面率。早期 n|U| 面产生率只作为诊断，不再作为当前 Q_i 的生成处方。电离层之外不加入体积产生率。

## 2. Forward tracing

### 2.1 方程与实际函数

位置 x 用 m、速度 v 用 m/s、时间 t 用 s、质量 m 用 kg、电荷 q 用 C、电场 E 用 V/m、磁场 B 用 T。O₂⁺ 的 q=+1.602176634×10⁻¹⁹ C，质量与上文一致。非相对论运动方程为：

$$
\frac{d\mathbf x}{dt}=\mathbf v,\qquad
\frac{d\mathbf v}{dt}=\frac qm\left(\mathbf E+\mathbf v\times\mathbf B\right).
$$

本例入口为 [ShellMonteCarlo.run_monte_carlo](monte_carlo_shell.jl)，内部由 `ShellMonteCarlo.trace_particle` 推进粒子，调用 `MarsTP.TP.boris_velocity_update` 完成 Boris 速度更新。使用静态 MHD 总电场和磁场，不能将此处电场简单理解为仅有 −U×B。半步速度用于积分，保存状态时转换为与位置同一时刻的速度。

| 参数 | 设置 |
| --- | --- |
| 火星半径 | Rm=3390 km |
| 释放高度、吸收内边界高度 | 500 km |
| 外边界火心距离 | 4 Rm |
| 时间步、最大飞行年龄 | 0.1 s、500 s |
| 忽略过程 | 体积产生、碰撞、复合、重力、反馈 |
| 终止类型 | inner、outer、time_limit、zero_rate；数值失败报错 |

无损传播期间每条轨迹的 Q_i 保持不变。此处所述是示例实际调用链，公共 `trace_forward` 并非这个源面模拟脚本的直接入口。

### 2.2 5000 条轨迹

![5000 trajectories in XZ, XY and YZ](trajectories_5000.png)

[plot_trajectories.py](plot_trajectories.py) 从正率且传播时间大于零的粒子中等概率选出 5000 条，NumPy 随机种子为 20260906；不是专门选择探头命中粒子，也未按 Q_i 加粗轨迹。每条轨迹最多显示 300 个均匀索引点，保留首尾，此显示降采样不改变原始积分或探头统计。

背景使用 `py_space_zc.maven.bs_mpb` 和 `plot_mars`。XZ、XY 中的 BS、MPB 是假定 +X 向阳的经验参考曲线，不是本次 MHD 提取的边界；YZ 不画这种 X–ρ 曲线。圆盘为火星，灰圆为 500 km 源面，洋红方框为探头投影。数据仍标为原生 MHD 笛卡尔坐标，尚未核实 MSO/MSE 来源。

## 3. 探头的三维 PSD 与二维 VDF

### 3.1 粒子率乘驻留时间

令 dN 为粒子计数，x 用 m、v 用 m/s、n 用 m⁻³；三维速度分布 f 的单位为 s³ m⁻⁶：

$$
dN=f(\mathbf x,\mathbf v)d^3x\,d^3v,\qquad n(\mathbf x)=\int f(\mathbf x,\mathbf v)d^3v.
$$

立方体边长 L=0.2 Rm=678,000 m，空间体积 V_D 用 m³；各速度 bin 宽度 Δv_x、Δv_y、Δv_z 用 m/s，速度空间体积 Δ³v_b 用 m³ s⁻³：

$$
V_D=L^3,\qquad \Delta^3v_b=\Delta v_x\Delta v_y\Delta v_z.
$$

令 τ_i,D,b 为驻留时间（s），T_i 为该轨迹的保存时长（s），x_i 用 m、v_i 用 m/s，空间区域 D 与速度 bin B_b 的指示函数为无量纲：

$$
\tau_{i,D,b}=\int_0^{T_i}\mathbf1_D[\mathbf x_i(t)]\mathbf1_{B_b}[\mathbf v_i(t)]\,dt.
$$

对稳恒源，Q_i 用 s⁻¹，τ_i,D,b 用 s，N_D,b 为平均粒子计数；V_D 用 m³、Δ³v_b 用 m³ s⁻³，空间和速度 bin 平均 PSD f̄_D,b 的单位为 s³ m⁻⁶：

$$
N_{D,b}=\sum_i Q_i\tau_{i,D,b},\qquad
\boxed{\bar f_{D,b}=\frac{\sum_i Q_i\tau_{i,D,b}}{V_D\Delta^3v_b}}.
$$

不再除以总积分时间，也不再除以粒子数，因为 Q_i 已含抽样归一化。此估计器不要求源一定是 Maxwellian，只要求率权重有物理含义。有限飞行年龄会截断驻留积分，因此 500 s 的结果并不自动证明稳态收敛。

### 3.2 与公共库完全一致的统计路径

- [forward_psd](../../../src/tracing/detector_psd_forward.jl)：从内存轨迹计算。
- [forward_psd_saved、foreach_saved_trajectory](../../../src/tracing/trajectory_io.jl)：从分批文件读取。
- 两者共用 [ForwardPSDAccumulator、accumulate_forward_psd!、finish_forward_psd](../../../src/tracing/forward_psd_accumulator.jl)。本例 [analyze_saved_probes.jl](analyze_saved_probes.jl) 一次扫描所有轨迹，同时更新三个探头。

每个相邻保存状态之间，位置和同步速度均线性插值。令 α、α_a、α_b 为无量纲段内比例，x_0、x_1 用 m、v_0、v_1 用 m/s、Δt 和 δt 用 s：

$$
\mathbf x(\alpha)=\mathbf x_0+\alpha(\mathbf x_1-\mathbf x_0),\quad
\mathbf v(\alpha)=\mathbf v_0+\alpha(\mathbf v_1-\mathbf v_0),\quad
\delta t=\Delta t(\alpha_b-\alpha_a).
$$

先将线段裁剪到立方体，再按速度 bin 边界拆分，累加各段的率乘时间。两端都在探头外、但中途穿过的段也会统计。探头空间采用下闭上开边界；速度末端全局上边界包含在最后一格。不会在探头交点重新做 Boris 更新，也无需加载 MHD。

三个速度轴均为 −500 到 500 km/s，bin 宽度 5 km/s，即 200³ 个 bin。三维 PSD 用稀疏结构保存非零值，省略项为本次抽样的零估计；不是跳过第三维直接做二维直方图。

```julia
using MarsTP
Rm = 3_390_000.0 # m
settings = (; detector_m=[0.,0.,2.] * Rm, side_m=0.2Rm,
             vlim=(-500.,500.), vgrid=200, velocity_unit=:km_s,
             species="O2+", coordinate_system="native MHD Cartesian axes")

# 已在内存中的轨迹：solutions 的 t/u 必须为同步的 s、m、m/s。
r3 = forward_psd(solutions; settings..., rate_weights_s=Q, option="3D")
# 约 50 GB 的保存目录：逐条读取，不将整个集合装入内存。
r3 = forward_psd_saved("outputs/my_run"; settings..., option="3D", storage=:sparse)
# 也可直接取得完整速度积分的二维结果。
rxy = forward_psd_saved("outputs/my_run"; settings..., option="Vx-Vy")
```

`vgrid=200` 是每轴 bin 数，范围 ±500 km/s 对应宽度 5 km/s。`storage=:dense` 保留原有数组接口；`:sparse` 返回一基 bin 元组到 PSD 值的 Dict。运行生成的 NPZ 转为零基 COO 索引。重复调用 `forward_psd_saved` 会重复读取磁盘，因此本例多探头使用下述一次扫描入口。

输出速度边界和中心始终为 m/s；`velocity_unit=:km_s` 仅指定输入 vlim 的单位。手动提供 t/u 时需保证物种和速度同步，不能对已经同步的速度再次修正半步。失败返回码、非有限数据或不递增时间会报错。`Terminated` 本身不能区分撞击或逃逸，不外推保存时间之外的轨迹。

### 3.3 完整积分得到二维图

令 f_ijk 为三维 PSD（s³ m⁻⁶），Δv_y=Δv_z=5000 m/s；二维 F_xy、F_xz 的单位为 s² m⁻⁵：

$$
F_{xy}(v_{x,i},v_{y,j})=\sum_k f_{ijk}\Delta v_{z,k},\qquad
F_{xz}(v_{x,i},v_{z,k})=\sum_j f_{ijk}\Delta v_{y,j}.
$$

求和覆盖被省略速度轴的完整 ±500 km/s 范围，**不是取 vz=0 或 vy=0 切片**。图的横纵轴仅显示 ±300 km/s，坐标换成 km/s，色标仍保持 SI 单位，不附加速度单位换算因子。笛卡尔速度网格不乘球坐标的 v² Jacobian。

令 n_D 为所选速度范围内的探头平均密度（m⁻³），f 为 s³ m⁻⁶，F 为 s² m⁻⁵，Δv 为 m/s，Q 为 s⁻¹、τ 为 s、V_D 为 m³：

$$
n_D=\sum_{ijk}f_{ijk}\Delta v^3
=\sum_{ij}F_{xy,ij}\Delta v^2
=\sum_{ik}F_{xz,ik}\Delta v^2
=\frac{\sum_iQ_i\tau_{i,D}}{V_D}.
$$

最后一个等号要求速度范围包含全部驻留贡献；库另外报告范围外密度。这三个结果的范围外密度均为零。重复进出会增加驻留时间，但不会成为新的独立源样本。

### 3.4 探头结果

中心 (1,0,2) Rm：

![Probe at 1 0 2](probe_psd_projections.png)

中心 (0,0,2) Rm：

![Probe at 0 0 2](probe_0_0_2/probe_psd_projections.png)

中心 (−1.5,0,1) Rm：

![Probe at minus1.5 0 1](probe_m1p5_0_1/probe_psd_projections.png)

所有图使用 turbo、Arial、对数色标；灰色为无抽样贡献，无平滑，不加底部说明文字。各图色标独立归一化，数值与样本统计列于下方复现记录。

### 3.5 权重噪声与穿面通量

令 a_i=Q_iτ_i,D 为粒子计数贡献，Q_i 用 s⁻¹、τ_i,D 用 s，有效样本数 N_eff 无量纲：

$$
a_i=Q_i\tau_{i,D},\qquad N_{\rm eff}=\frac{(\sum_i a_i)^2}{\sum_i a_i^2}.
$$

特别是 (0,0,2) 探头的 N_eff 约 9.7，细结构的抽样噪声仍很大。

穿面事件统计是另一个量。令 A_face 为 m²，Q_i 为 s⁻¹，指定面和方向的穿面事件集合为 C，通量 Γ 为 m⁻² s⁻¹：

$$
\Gamma=\frac{\sum_{i\in C}Q_i}{A_{\rm face}}.
$$

不要用穿面事件计数替代驻留时间密度。以匀速垂直穿过立方体的束流为例，n 用 m⁻³、u 用 m/s、L 用 m、Q 用 s⁻¹、τ 用 s：

$$
Q=nuL^2,\qquad \tau=L/u,\qquad \frac{Q\tau}{L^3}=n.
$$

这说明速度较大的粒子虽然穿面率更大，但驻留时间更短，两者共同恢复正确密度。

## 4. 分批保存轨迹与复现

[write_trajectory_batch](../../../src/tracing/trajectory_io.jl) 是可复用的库函数。每批只保留当前批次的轨迹，逐批写入后释放内存；总输出可以很大，内存不随全部历史状态线性增长。

```julia
using MarsTP
write_trajectory_batch("outputs/my_run/trajectories_00001.jld2", trajectory_batch;
    particle_ids=ids, rate_weights_s=Q,
    source_density_weights_m3=density_weights, cell_ids=cell_ids,
    termination_codes=statuses, species="O2+",
    coordinate_system="native MHD Cartesian axes", compress=false)
```

必需参数是 `particle_ids` 和物理粒子率 `rate_weights_s`，其余诊断可省略。函数拒绝覆盖现有路径。文件保存 Float64 的 `p<ID>/state`，Julia 为 7×N，h5py 为 N×7，顺序 t,x,y,z,vx,vy,vz，单位 s、m、m/s。每个粒子同时保存率权重；文件还保存格式版本、单位、种类、坐标说明、终止类型及批次完成标记。`compress=true` 可无损压缩，默认关闭。中途失败的文件不带完成标记，读取器会拒绝使用。

完整的源采样、追踪和保存入口是 [ShellMonteCarlo.run_monte_carlo](monte_carlo_shell.jl)：

```julia
include("examples/forward_tracing/monte_carlo_forward_tracing/monte_carlo_shell.jl")
ShellMonteCarlo.run_monte_carlo("outputs/new_rate100_run",
    ShellMonteCarlo.Config(per_cell=100, dt=0.1, tmax=500., batch_size=1024,
                          flux_model="reservoir_maxwellian_rate", compress_trajectories=false))
```

该函数实际调用 `write_trajectory_batch`，不是另一套内嵌写盘实现。原始运行写出了约 48.3 GB 文件；若仍保留这些文件，读取接口兼容原来的 legacy p<ID> 格式，无需重写。旧文件未保存的逐粒子终止码返回 `unavailable`，不伪造成功状态；原 `particles.csv` 仍保留实际终止原因。

### 文件与复现

重新生成教学图（仅 PNG，不运行轨迹模拟）：

```powershell
& C:\Users\Win\.conda\envs\mars\python.exe examples/forward_tracing/monte_carlo_forward_tracing/plot_monte_carlo_sampling.py
```

在仓库根目录运行，使用本仓库 Project.toml / Manifest.toml。Python 仅用于 PNG 和可移植 NPZ 输出，需要 NumPy、Matplotlib、h5py。轨迹示意图另用 `py_space_zc`。

首次模拟需要 `data/mars_fields_spherical_from_dat.vts`，SHA256 为 `fa92bd82fe16975ad0d50f4e40ace344e9d41389d6024976d423e6126bc7a5c8`。对已保存轨迹计算 PSD 不需要重新读取 MHD。原始模拟为 Julia 1.12.6、MarsTP 提交 `58e078fa4d853a0307df0e9c87b9b44528835dc7`。

```powershell
$example = 'examples/forward_tracing/monte_carlo_forward_tracing'
$run = 'outputs/mc500_rate100_full_20260906a'
$reprobe = 'outputs/new_library_psd'
$py = 'C:\Users\Win\.conda\envs\mars\python.exe'
# 一次读取全部批次，并计算三个探头。输出目录必须不存在。
julia --startup-file=no --compiled-modules=existing --project=. "$example/analyze_saved_probes.jl" $run $reprobe
# 只画库函数已经计算好的 PSD，不再在 Python 中分箱。
& $py "$example/plot_library_psd.py" $reprobe
```

需要从头重新模拟时，先设置新的 `$run`，再执行：

```powershell
$env:MC_GIT_COMMIT = git rev-parse HEAD
$env:MC_GIT_STATUS = (git status --short) -join "`n"
julia --startup-file=no --compiled-modules=existing --threads=12 --project=. "$example/monte_carlo_shell.jl" $run 1 500 0.1 100 reservoir_maxwellian_rate false
```

只要保留完整轨迹，重新分 bin 或增加探头就不需要重跑原始积分；若已删除轨迹，则需要重新模拟。可先用 cell stride=36、tmax=0.2、dt=0.1、per_cell=2 做小样本写盘检查。

| 文件 | 用途 |
| --- | --- |
| [analyze_saved_probes.jl](analyze_saved_probes.jl) | 一次读取全部轨迹，三个公共累积器计算 PSD 和探头记录 |
| [plot_library_psd.py](plot_library_psd.py) | 读取库计算的 JLD2，核验积分、输出 PNG/NPZ/JSON |
| [monte_carlo_shell.jl](monte_carlo_shell.jl) | 源采样、Boris、公共分批保存函数 |
| [plot_trajectories.py](plot_trajectories.py) | 5000 条轨迹示意图 |
| [test_monte_carlo.jl](test_monte_carlo.jl) | 源模型、几何与积分测试 |
| [trajectory_io.jl](../../../test/trajectory_io.jl) | 磁盘与内存结果一致、压缩、legacy、异常输入测试 |

`reprobe_saved.py` 兼容入口转发到 Julia；`analyze_monte_carlo.py`、`analyze_probe.py` 的命令行转发到库结果绘图。其中旧 Python 分箱函数只保留为独立回归参考，不用于当前发布结果。`synchronize_reprobe.jl` 已停用，会明确提示使用新的公共入口，防止误用旧 Boris 交点处理。

每个探头输出 `library_psd.jld2`、`library_summary.toml`、`probe_residence.csv`，绘图步骤另存 PNG、`probe_psd_sparse.npz`、`analysis_summary.json`、逐粒子权重/驻留时间 CSV、穿面通量/速度 CSV。完整扫描成功后才写入 `analysis_complete.toml`。运行生成的本地 NPZ 包含三维非零值、零基索引、完整边界和二维投影：

```python
import numpy as np
with np.load('probe_psd_sparse.npz') as a:
    ix, iy, iz = a['indices_xyz'].T
    fxy = np.zeros(tuple(a['shape_xyz'][:2]))
    np.add.at(fxy, (ix, iy), a['f3d_s3_m6'] * float(a['dv_ms']))
    np.testing.assert_allclose(fxy, a['fxy_s2_m5'])
```

### 已保存轨迹上的三个探头

| 探头位置 (Rm) | 密度 (cm⁻³) | 独立命中粒子 | 驻留权重有效样本数 | 驻留段 |
| --- | ---: | ---: | ---: | ---: |
| (1,0,2) | 0.0309669 | 416 | 119.29 | 11,151 |
| (0,0,2) | 0.001002168 | 115 | 9.71 | 3,328 |
| (−1.5,0,1) | 0.1177217 | 2,011 | 353.71 | 124,374 |

当前图来自 `mc500_library_psd_20260906`，公共接口完整读取原运行的 992 个批次和 1,015,300 个粒子组。三个探头均满足 `sum(f3d)*dv³ = sum(fxy)*dv² = sum(fxz)*dv² = sum(Q*tau)/V`，速度范围外密度为零。色标独立归一化，不宜只凭颜色跨图比较。特别是 (0,0,2) 的有效样本数约 9.7，5 km/s 细结构仍受抽样噪声影响。

相较以前的 Boris 交点重算版本，密度与命中数不变，六个二维投影的相对 L1 差异最大 2.48×10⁻⁶，约 0.00025%。这一比较为本次方法统一时的验证结果。

## 5. 验证与限制

以下模拟与库测试为已有运行记录；本次文档合并仅重新验证教学采样及链接，未重跑百万粒子模拟。

稳态公式要求稳恒源、静态场和收敛的驻留积分。持续注入永久束缚且无损失的粒子时，有限稳态可能不存在；时变源或场需要明确释放与观测时刻。若加入有依据的损失或电荷态变化，应使用对应存活权重及物种轨迹段。收敛应分别检查样本数、探头尺寸、速度 bin、时间步与最大飞行年龄，误差按独立源样本或随机种子批次估计。

- 包内全部 312 项测试通过，其中 25 项覆盖分批读写及内存/磁盘 PSD 等价；示例 226 项测试通过。
- 使用新 `run_monte_carlo` 和 `write_trajectory_batch` 完成 566 个样本、0.2 s 的实际 MHD 小规模写盘运行；未重写原始 48.3 GB。
- 原运行正率 506,259 个、零率 509,041 个；inner 315,992、outer 108,153、time_limit 82,114，无数值失败。12 线程积分写盘约 589.5 s。
- 最大飞行年龄 500 s，达到时限的源率占 18.13%，未证明稳态或完整 PSD 收敛。以前 231 条子集时间步细化改变探头密度约 −0.00107%、−0.00218%，不代表完整集合收敛。
- 原轨迹整体最大做功闭合残差约 0.173 eV；320 条相对残差超过 1%（分母以 1 eV 为下限）。是否仍可重算其他探头，取决于本地是否保留完整轨迹。

```powershell
julia --startup-file=no --compiled-modules=existing --project=. test/runtests.jl
julia --startup-file=no --compiled-modules=existing --project=. "$example/test_monte_carlo.jl"
```

环境记录：仅将 Manifest 中已有的标准库 TOML 声明为直接依赖，保留所有锁定版本。现有本地 Registry 缺少锁定的 UnsafeAtomics 0.3.2，`Pkg.resolve()` 未完成；未因此升级包。已安装环境可完成上述测试和运行。
