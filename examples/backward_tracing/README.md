# Backward tracing

本目录包含固定探头反向轨迹、薄层面源推导及对 Vy 积分的二维 VDF。示例使用火星半径 3390 km；所有命令从仓库根目录运行。

## 1. 随机速度的反向轨迹

![5000 个粒子的 XZ 轨迹](images/probe_backtrace_xz_5000.png)

探头位于 `(0,0,2 Rm)`，5000 个初始速率在 10 至 200 km/s 均匀采样，初始 Vy=0，方向在 XZ 平面均匀采样。这张图仅展示路径，不代表按物理源项加权的 VDF。

```sh
python examples/backward_tracing/plot_probe_backtrace.py
```

入口调用 [random_probe_backtrace.jl](random_probe_backtrace.jl)。已有轨迹可以用 `--data outputs/<run_id>/trajectories.jsonl` 重绘。详细条件和依赖见 [轨迹说明](probe_backtrace.md)。

## 2. 面源和相空间密度推导

[完整推导与单位](backtracing_derivation.md) 对应 400 km 可穿越薄层面源与固定 200 km 内边界。每次穿层累计 `F*g/abs(v dot er)`，穿层后继续积分体积源。核心实现位于仓库 `src/tracing/detector_psd_backward.jl`，MHD 读取位于 `src/data/mhd.jl`。

## 3. Vy 积分的二维速度分布，尚未收敛

![初步二维 VDF](images/vdf_xz_preliminary.png)

这是 `f_xz=integral f dVy`，单位 s² m⁻⁵，与 Vy=0 切片不同。首轮网格共 35301 条轨迹，**速度网格和时间步长检查未收敛，不能用于定量科学结论**。详细设置与检查见 [初步 VDF 说明](vdf_preliminary.md)。

重绘随示例保留的数值数据：

```sh
python examples/backward_tracing/plot_probe_vdf.py
```

`data/` 保存首轮 `vdf.csv`、五个检查点、步长/速度间隔检查结果和绘图元数据。原始 JLD2 和运行日志仍保留于本地 outputs。

重新计算（约 35301 条轨迹），请为每次运行使用新的输出目录：

```sh
julia --startup-file=no --threads=4 --project=. examples/backward_tracing/probe_vdf_integrated.jl outputs/probe_vdf_new_run
python examples/backward_tracing/plot_probe_vdf.py outputs/probe_vdf_new_run
```

检查脚本 [check_probe_vdf_grid.jl](check_probe_vdf_grid.jl) 读取指定输出目录的 `qa_points.csv`，计算高值点在不同 Vy 和时间步长下的结果；随示例的数据目录提供本轮检查点。该脚本依赖同样的 MHD、MAT 体积源及 GITM 输入。

## 4. 对流和 Hall 电场净增能

![全路径电场净增能](images/energy_gain_xz_preliminary.png)

[符号、单位和计算说明](energy_gain.md)。保存每个三维速度网格点的对流、Hall、总电场做功，并按 PSD 对 Vy 加权绘制二维均值。正值为过去到探头的净增能，负值为净失能。数据包含能量闭合检查与异常轨迹步长细化结果；速度网格仍未收敛。

## 5. 沿路径的局部功率与累计做功

![沿轨迹的对流、Hall 和总电场功率](images/probe_path_power_xz_5000.png)

[计算定义和运行说明](path_power.md)。复用原 5000 条轨迹初始条件，XZ 投影第一行表示局部功率（eV/s），第二行表示从回溯终点到当前位置的累计做功（eV）。红色增能，蓝色失能，采用 coolwarm 和对称对数色标。

## 6. 尾区探头 (-2,0,0) Rm

[新算例说明](tail_probe.md)：5000 个 O2+，速率在 0–100 km/s 均匀采样，三维方向各向同性。提供 XZ 轨迹图，以及对流、Hall、总电场的局部功率与累计做功两行图。

## 依赖

Julia 使用项目环境。随机轨迹绘图需 Chi Zhang 自定义 py_space_zc（bs_mpb、plot_mars 及火星贴图）、NumPy、Matplotlib、Arial；VDF 绘图只需 NumPy、Matplotlib、Arial。物理运行需仓库 data 中的相应输入，不自动下载或安装。

返回 [示例总目录](../README.md)。
