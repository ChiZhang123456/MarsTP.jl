# O₂⁺ 轨迹的 XZ 与 YZ 投影

![XZ and YZ particle trajectories](images/trajectories_xz_yz_800km.png)

## 图的布局

3 行、2 列，第一列为 XZ 投影，第二列为 YZ 投影。三行依次为全部粒子、起始位置在向阳面（`X0 > 0`）的粒子、起始位置在背阳面（`X0 <= 0`）的粒子。

每行两列显示同一批粒子的完整三维轨迹，仅投影方向不同。粒子随后跨过日夜分界面不会改变所属分组。所有面板采用相同尺度和 Arial 字体，标题分别为 `All O2+`、`Dayside O2+`、`Nightside O2+`；不显示主标题、底部图例或说明文字。

- 蓝色：到达外边界；橙色：返回内边界。若后续参数产生达到时间上限的粒子，则使用紫色。
- 深色小点：释放位置。细灰圆：释放球面的投影轮廓。
- XZ 列用 `py_space_zc.maven.bs_mpb` 绘制 BS（虚线）和 MPB（点线）参考曲线，用 `plot_mars` 绘制火星。
- YZ 列只用 `plot_mars(texture=False)` 绘制几何圆盘，不将 XZ 边界曲线误用到 YZ 平面，也不指定未经核实的表面贴图方向。
- 曲线穿过火星圆盘的二维投影不代表实际粒子进入火星。XZ 投影包含所有 Y，YZ 投影包含所有 X。

## 运行代码

绘图入口：[plot_trajectories_xz_yz.py](plot_trajectories_xz_yz.py)。它直接调用本目录已上传的 Julia 积分代码 [sphere_trajectories.jl](sphere_trajectories.jl)。

在仓库根目录，用已配置 `py_space_zc` 的 Python 执行：

```sh
python examples/plot_trajectories_xz_yz.py
```

默认输出为 `examples/images/trajectories_xz_yz_800km.png`。要保留已上传的图片，可在运行前指定新路径，例如 PowerShell：

```powershell
$env:TRAJECTORY_PREVIEW = "$PWD/xz_yz_preview.png"
python examples/plot_trajectories_xz_yz.py
```

脚本使用 4 个 Julia 线程。Julia 在 PATH 中可用，Python 需要 NumPy、Matplotlib 和用户的自定义 `py_space_zc` 库。完整依赖说明见 [examples README](README.md#环境与输入数据)。运行需要本地的 `data/mars_fields_spherical_from_dat.vts`，该大文件不包含在 GitHub 中。

轨迹数据仅通过标准输出传入 Python 内存，不写入数据文件；只输出图片，因此每次运行都会重新积分。

## 参数与核对结果

沿用 800 km 球面释放示例：1000 个 O₂⁺，三个初始速度分量均为零，随机数为 `Xoshiro(20260905)`，按固体角均匀采样。积分步长 0.1 s，火星半径 3390 km，吸收内边界为高度 200 km，外边界为火心距离 `4 Rm`。20,000 s 为时间保护上限。电磁场和物理假设见 [参数说明](README.md#初始条件和边界)。

| 起始位置分组 | 粒子数 | 到达外边界 | 返回内边界 |
|---|---:|---:|---:|
| All | 1000 | 736 | 264 |
| Dayside | 500 | 453 | 47 |
| Nightside | 500 | 283 | 217 |

本次全部粒子均触边，最长飞行时间约 4148.53 s。脚本检查输入分组数量、有限值和最终三维边界半径。图片是轨迹可视化，不代表源项加权后的粒子通量或完整数值收敛验证。
