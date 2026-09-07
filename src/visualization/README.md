# 火星边界与轨迹绘图

Python 绘图模块，依赖 NumPy、Matplotlib、h5py 和 Python 3.11+，无需 `py_space_zc`。在仓库根目录执行以下命令。

## 轨迹图

`plot_trajectory.py` 读取轨迹文件，调用 `mars.py` 绘制火星、BS 和 MPB，输出 PNG。

```powershell
# Forward：分批保存的 JLD2 目录，选择 O2+，画三个二维投影。
python src/visualization/plot_trajectory.py outputs/my_run --species O2+ --count 5000 --output outputs/forward.png

# Backward：points 用 Rm、times_s 用 s 的 JSONL 文件。
python src/visualization/plot_trajectory.py outputs/my_backtrace/trajectories.jsonl --position-unit Rm --species O2+ --planes XZ XY YZ 3D --output outputs/backward.png

# 指定物种和计算粒子 ID。
python src/visualization/plot_trajectory.py outputs/my_run --species O2+ --particle-ids 12 45 91 --planes 3D --output outputs/selected.png
```

`--species O2+` 筛选离子种类；`--particle-ids` 选择保存时的粒子编号。物种信息按逐粒子字段、文件字段、旁边的 `metadata.toml` 或 `metadata.json` 读取。混合物种文件可以给每条轨迹单独标记物种。缺少物种标签的旧数据，使用 `--assume-species O2+` 显式声明其物种；此参数不会改写已有标签，也不会把其他物种转换成 O₂⁺。

`--count` 控制最多显示的轨迹数，默认 5000，使用固定随机种子的等概率蓄水池抽样；`--max-points` 控制每条曲线最多显示的点数，默认 300，并保留首尾点。文件名与粒子 ID 共同标识轨迹，传入单个文件可进一步限定选择范围。

### 文件格式

- `write_trajectory_batch` 输出的 JLD2：`p<ID>/state` 经 h5py 读取为 N×7，列顺序为 t,x,y,z,vx,vy,vz，单位 s、m、m/s。支持相同结构的 HDF5 文件。读取时保留时间方向，递增为 forward，递减为 backward。
- JSONL：每行一个对象，至少包含 `id`、`points`（N×3）和 `times_s`（N 个值），可包含 `species`、`coordinate_system`、`position_unit`。位置单位为 `m`、`km` 或 `Rm`，可通过 `--position-unit` 指定。
- 只有 PSD 或终点的文件不包含完整轨迹，需要在追踪时另存路径。

默认 Rm=3,390,000 m，可用 `--rm-m` 修改。最终绘图坐标统一为 Rm。JLD2 只读取被抽中的轨迹及显示点，不将全部轨迹状态装入内存。

### Python 调用

```python
import sys
sys.path.insert(0, "src")
from visualization import plot_trajectory

fig, axes, records = plot_trajectory(
    "outputs/my_run", species="O2+", count=100,
    planes=("XZ", "3D"), output="outputs/trajectory.png")
```

## 火星、BS 和 MPB

`plot_mars` 在二维坐标轴默认使用随仓库附带的 [火星图片](mars_globe_true_color.png)，无需安装 `py_space_zc` 或指定图片路径。图片来自用户提供的 `py_space_zc.maven` 本地资源，文件原样复制。使用 `texture=False` 可画纯色圆盘，`texture_path=...` 可指定自己的图片。

三维坐标轴默认绘制球体。这张图片是火星圆盘照片，不是经纬度展开的全球纹理，因此不直接贴到三维球面。

`bs_mpb` 使用与 `py_space_zc.maven.bs_mpb` 相同的圆锥曲线参数。x、ρ、r、x₀、L 均以 Rm 表示，θ 为 rad，偏心率 ε 无量纲：

$$
r=\frac{L}{1+\epsilon\cos\theta},\qquad
x=x_0+r\cos\theta,\qquad \rho=r\sin\theta.
$$

| 边界 | x₀ (Rm) | L (Rm) | ε | 区域 |
| --- | ---: | ---: | ---: | --- |
| BS | 0.600 | 2.081 | 1.026 | 全部 |
| MPB 日侧 | 0.640 | 1.080 | 0.770 | x≥0 |
| MPB 夜侧 | 1.600 | 0.528 | 1.009 | x<0 |

只使用正半径分支。三维图绕 X 轴旋转生成边界曲面，XY、XZ 绘制对称截线。YZ 可用 `--x-slice 0` 显示 X=0 Rm 的边界截面，默认不画 YZ 边界。

模型假设 +X 向阳。若数据坐标不符合该约定，使用 `--no-boundaries`，或先转换轨迹坐标。

```python
import matplotlib.pyplot as plt
from visualization import plot_mars_context, bs_mpb, plot_mars

fig, ax = plt.subplots()
plot_mars_context(ax, plane="XZ", xmin=-4.2)
ax.set(xlim=(-4.2, 4.2), ylim=(-4.2, 4.2), aspect="equal")

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
bs_mpb(ax, plane="3D")
plot_mars(ax)
ax.set_box_aspect((1, 1, 1))
```
