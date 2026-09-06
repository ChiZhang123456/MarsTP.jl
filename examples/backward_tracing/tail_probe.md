# 尾区探头：(-2, 0, 0) Rm 的 5000 个 O2+

本例将探头移到模型笛卡尔坐标 (-2,0,0) Rm，Rm=3390 km。按用户确认的分布，速率在 0 至 100 km/s **均匀采样**，方向在三维球面各向同性。不是三个分量独立均匀采样，也不是在速度球体内按体积均匀采样。

用固定种子 20260905 的 Julia Xoshiro，采样 speed~Uniform(0,100 km/s)、mu=cos(theta)~Uniform(-1,1)、phi~Uniform(0,2pi)，然后取 v=speed*(sqrt(1-mu²)cos(phi),sqrt(1-mu²)sin(phi),mu)。初始 Vy 可以非零，实际积分为三维运动，绘图仅为 XZ 投影。

## 图与能量符号

![XZ 轨迹](images/tail_probe_trajectories_xz_5000.png)

轨迹图颜色为初始速率，使用 turbo。星号是探头，火星及参考边界使用 py_space_zc.maven.plot_mars、bs_mpb，所有非数学字体为 Arial。

![局部功率与累计做功](images/tail_probe_path_power_xz_5000.png)

两行三列分别为对流、Hall、总电场：第一行是每个 0.5 s 显示段的平均功率 qE·v（eV/s）；第二行是从回溯终点沿正时间累积到当前位置的做功（eV）。正值红色代表增能，负值蓝色代表失能。两行使用 coolwarm 和 SymLogNorm，零附近分别在 ±1 eV/s、±1 eV 内线性，三个场项在同一行共享色标。没有空间分箱；重叠的投影线不是空间平均。

回溯使用负 dt，但电荷和物理速度不翻转。每步正时间做功为 -qE·v dt，累计量在回溯终点为零，在探头处等于整条已计算路径的净做功。总做功与 K_probe-K_endpoint 检查闭合。详细公式见 [局部与累计说明](path_power.md)。这些图不按源项 PSD 加权，也不把回溯终点认定为粒子出生位置。

## 物理输入与终止

使用仓库 data/mars_fields_spherical_from_dat.vts 的静态总场进行积分，对流、Hall 做功在同一条总场轨迹上分解。保留 O2+ 正电荷、固定 200 km 吸收边界、4 Rm 外边界、500 s 回溯上限。Boris 步长 -0.05 s。图中源层不终止轨迹，此路径示例不计算 400 km 面源或体积源的 VDF 权重。

**达到 500 s 上限的轨迹仅表示有限历史的累计做功，不能称为完整边界到探头的能量变化。** 边界状态与时间上限状态分别保存在数据中。

## 代码与复现

复用 [probe_path_power.jl](probe_path_power.jl) 和 [plot_path_power.py](plot_path_power.py)，新增 tail_isotropic 模式。旧模式仍是默认，保留原探头 (0,0,2) Rm 的设置。所有命令从仓库根目录运行，使用已安装的 Julia 项目依赖和 Python 的 NumPy、Matplotlib、py_space_zc。输入场不自动下载。

```sh
julia --startup-file=no --compiled-modules=existing --threads=1 --project=. examples/backward_tracing/probe_path_power.jl outputs/tail_probe_new_run 5000 tail_isotropic
python examples/backward_tracing/plot_path_power.py outputs/tail_probe_new_run
```

小样本：将 5000 改为 5；末尾追加 -0.025 可检查半步长。每次使用新输出目录。case.json 保存探头、采样模式、随机种子、粒子数与步长。segments.csv 保存局部和累计绘图所需段数据，particles.csv 保存各粒子初始速度、终止状态、净做功、动能差及残差。大型分段文件保留本地，可按上述命令重算；GitHub 保存 PNG、逐粒子结果、检查数据和代码。

此前 (0,0,2) Rm 的结果和代码仍见 [原算例说明](path_power.md)。两个算例的图文件名不同，不覆盖此前结果。

## 本次结果与检查

5000 条中，935 条到达 200 km 内边界，3605 条到达 4 Rm 外边界，460 条达到 500 s 上限。没有非有限状态。抽样速率范围 0.0171 至 99.9785 km/s，均值 50.4727 km/s。三方向单位向量平方的均值分别为 0.3392、0.3347、0.3261，与各向同性预期 1/3 相近。

全样本最大能量闭合残差 1.1545 eV。以 max(abs(work),abs(delta K),1 eV) 为分母，99 百分位相对残差为 0.0104%。其中两个低净能量变化粒子的相对残差超过 1%，已单独标记并减小步长复算，不能只根据接近零的净做功符号作出精确结论。半步长小样本的五条轨迹保持终止状态一致，最大总/对流/Hall 做功变化分别为 0.04481、0.04682、0.00201 eV。检查仅覆盖选定样本，不代表全体轨迹和显示间隔已经收敛。

- [每粒子净能量结果](data/tail_probe_particles.csv)
- [采样配置](data/tail_probe_case.json)
- [采样及半步长检查](data/tail_probe_sampling_step_qa.json)
- [五粒子半步长原始结果](data/tail_probe_fine_pilot.csv)
- [需细化的速度点](data/tail_probe_refine_points.csv)

逐粒子数据状态码 1=时间上限，2=内边界，4=外边界。末尾 residual_eV 为 total_eV-deltaK_eV。正净做功为增能，负值为失能。

[两个速度点的细化结果](data/tail_probe_energy_refinement.csv)：dt=0.00625 s 时，最大绝对闭合残差 0.00443091 eV，最大相对残差 0.0462687%。原始粗步长数据保留，未静默替换。可运行 `julia --startup-file=no --compiled-modules=existing --project=. examples/backward_tracing/refine_energy_gain.jl outputs/tail_probe_new_run tail`，输入目录需先包含 refine_points.csv。

[绘图和累计端点检查](data/tail_probe_path_power_qa.json)：2330378 个显示段，所有粒子的累计端点与净做功最大差为 3.1e-11 eV。仅发布 PNG。
