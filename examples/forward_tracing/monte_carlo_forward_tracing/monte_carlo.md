# Monte Carlo forward tracing：立方体探头的三维速度分布与相空间密度

## 1. 核心结果

对于稳恒源，每条 forward trajectory 携带粒子率权重 $Q_i$，即 `rate_weights_s[i]`，单位为 particles s⁻¹。立方体探头内第 $b$ 个三维速度 bin 的分布函数为

$$
\boxed{\bar f_{D,b}=\frac{\sum_i Q_i\tau_{i,D,b}}{V_D\Delta^3v_b}}
$$

其中 $\tau_{i,D,b}$ 是轨迹 $i$ **同时位于探头空间内、且速度位于 bin $b$ 内**的累计驻留时间。先用“粒子率 × 驻留时间”得到稳态粒子数，再除以空间体积和速度空间体积。这里不再除以总积分时间，也不再除以粒子总数，因为源采样权重已经含有 Monte Carlo 归一化。

此式适用于有确定物理粒子率权重的不同源模型，不要求一定采用 Maxwellian。源权重决定注入什么，驻留时间估计器决定如何统计探头内的分布。

## 2. 所求物理量与单位

本文使用速度空间分布函数 $f(\mathbf x,\mathbf v)$，定义为

$$
dN=f(\mathbf x,\mathbf v)\,d^3x\,d^3v,
\qquad n(\mathbf x)=\int f(\mathbf x,\mathbf v)\,d^3v.
$$

立方体边长 $L$ 用 m 表示，体积 $V_D=L^3$，单位 m³；三维笛卡尔速度 bin 体积为

$$
\Delta^3v_b=\Delta v_{x,b}\Delta v_{y,b}\Delta v_{z,b}.
$$

采用 m/s 时，其单位为 m³ s⁻³，因此

$$
[\bar f_{D,b}]=\frac{\mathrm{s^{-1}}\,\mathrm{s}}
{\mathrm{m^3}\,\mathrm{m^3\,s^{-3}}}=\mathrm{s^3\,m^{-6}}.
$$

严格来说结果是立方体和速度 bin 内的平均值：

$$
\bar f_{D,b}=\frac{1}{V_D\Delta^3v_b}
\int_D d^3x\int_{B_b}d^3v\,f(\mathbf x,\mathbf v).
$$

只有当探头足够小、分布在其内部变化不大时，才可近似称为探头中心“某个点”的 PSD。探头越小，统计通常越差，需要增加采样。

| 统计量 | 公式 | SI 单位 |
| --- | --- | --- |
| bin 内平均粒子数 | $N_{D,b}=\sum_iQ_i\tau_{i,D,b}$ | particles |
| bin 内密度贡献 | $n_{D,b}=N_{D,b}/V_D$ | m⁻³ |
| 三维 PSD | $f_{D,b}=n_{D,b}/\Delta^3v_b$ | s³ m⁻⁶ |
| 探头总密度 | $n_D=\sum_bf_{D,b}\Delta^3v_b$ | m⁻³ |

速度 bin 如果用 km/s，直接相除得到的是 m⁻³ (km/s)⁻³，其数值等于 SI PSD 数值乘以 $10^9$。建议计算时全部用 m/s，仅绘图坐标转 km/s。若密度也改用 cm⁻³，则 cm⁻³ (km/s)⁻³ 的数值为 SI PSD 数值乘以 $10^3$。对 Cartesian velocity bins 不需要额外的 $v^2$ 因子；球坐标速度网格才需要相应 Jacobian。

## 3. 使用接口与数值实现

`forward_psd` 接收 `trace_forward` 返回的结果，也支持 TestParticle ensemble、轨迹数组或单条含 `t` 和 `u` 的轨迹。示例：

```julia
using MarsTP

# solutions = trace_forward(source.initial_states; config)
# 每条轨迹必须与源采样时保存的粒子率权重一一对应。
settings = (
    detector_m = [-1.5Rm, 0.0, 1.0Rm], # 与轨迹相同的笛卡尔坐标系，m
    side_m = 100e3,                     # 立方体边长，m
    vlim = (-100.0, 100.0),             # 三个速度方向共用的边界
    vgrid = 80,                        # 每个速度方向的 bin 数，不是步长
    velocity_unit = :km_s,              # 仅指定输入 vlim 的单位
    species = "O2+",
    rate_weights_s = source.rate_weights_s, # particles/s，不能使用 density_weights_m3
    coordinate_system = "MSO",         # 仅在实际轨迹采用 MSO 时这样设置
)
result3d = forward_psd(solutions; settings..., option="3D")
result_xy = forward_psd(solutions; settings..., option="Vx-Vy")
result_yz = forward_psd(solutions; settings..., option="Vy-Vz")
result_xz = forward_psd(solutions; settings..., option="Vx-Vz")

f_xyz = result3d.psd # 数组维度顺序为 (Vx, Vy, Vz)，单位 s^3 m^-6
f_xz = result_xz.psd # 数组维度顺序为 (Vx, Vz)，单位 s^2 m^-5
vx, vz = result_xz.velocity_centers_m_s # 始终返回 m/s
n_covered = result3d.density_in_range_m3
n_outside = result3d.density_outside_vlim_m3
n_total = result3d.density_total_m3     # 包括速度范围外的贡献
```

`vlim=100.0` 表示对称边界 `(-100,100)`，也可以为三个轴分别设置边界，例如 `vlim=((-100,100),(-50,50),(-80,80))`。`vgrid=(80,40,64)` 设置各轴 bin 数。默认 `velocity_unit=:m_s`，输出始终为 SI。`option` 也接受 `:xyz`、`:xy`、`:yz`、`:xz`。

二维结果在指定的第三轴速度范围内积分，与对应三维结果乘以该轴 bin 宽度后求和一致；若要恢复完整二维分布，须确保三个轴的速度范围足够宽。二维选项直接累计二维数组，不分配完整三维数组。

粒子种类用于检查已有 TestParticle 轨迹的质量和电荷参数，并记录到结果中，不会重新计算轨迹。每次调用只处理同一物种。手动构造的 `t/u` 数据需由调用者保证物种和速度同步。当前 TestParticle 0.23.3 的 Boris 和 AdaptiveBoris 保存输出已经同步速度，不要再次做半步修正。

结果还提供 `residence_s`、`outside_vlim_residence_s` 和 `retcodes`，分别记录每条轨迹的探头内总驻留时间、速度范围外驻留时间和原始求解器返回码。`Terminated` 仅保留其原始含义，不能据此区分撞击或逃逸；达到积分末时刻也不等于物理逃逸。失败返回码、非有限数据和不递增时间会报错。函数只统计保存时间段，不外推终止后的轨迹，也不保存文件。

实现步骤：

1. 设置探头中心、边长和三个速度轴的 bin edges，记录实际场文件的坐标系。所有单位先统一为 SI。
2. 保留每条轨迹对应的 $Q_i$。位置和速度需要同步到同一时刻，Boris 半步速度不能未经同步直接用于速度 bin。
3. 对每个积分步求轨迹与立方体的步内交集。即使两个端点都在探头外，轨迹也可能穿过探头，不能仅检查保存点是否在内部。
4. 对探头内的时间段，进一步按速度 bin 边界切分。粒子在探头内被加速时，不能把整段驻留时间全部记入入口速度 bin。
5. 对同时属于空间探头和速度 bin 的每个子段累加 `hist[b] += Q_i * dt_sub`。
6. 最后计算 `f[b] = hist[b] / (L^3 * dvx[b_x] * dvy[b_y] * dvz[b_z])`。

以下伪代码概括统计过程，`phase_space_segments` 代表几何求交及速度分箱步骤，实际实现位于 `forward_psd` 内部：

```julia
occupancy = zeros(length(vx_edges)-1, length(vy_edges)-1, length(vz_edges)-1)
for i in eachindex(trajectories)
    Q = rate_weights_s[i]
    for step in consecutive_steps(trajectories[i])
        for (ix, iy, iz, dt_inside) in phase_space_segments(step, cube, velocity_edges)
            occupancy[ix, iy, iz] += Q * dt_inside
        end
    end
end
for iz in axes(occupancy,3), iy in axes(occupancy,2), ix in axes(occupancy,1)
    dv3 = diff(vx_edges)[ix] * diff(vy_edges)[iy] * diff(vz_edges)[iz]
    f[ix,iy,iz] = occupancy[ix,iy,iz] / (L^3 * dv3)
end
```

局部线性近似下，可先用 slab intersection 求空间进入和离开的参数 $\alpha\in[0,1]$，再求 $v_k(\alpha)$ 与各速度 bin edge 的交点。排序这些交点后，用每段中点判定 bin，用 $\Delta t(\alpha_2-\alpha_1)$ 累加驻留时间。该方法只对采用的线段与线性速度模型准确，真实曲线误差仍需步长收敛。也可使用积分器 dense output 配合事件求根或充分细的子步。

采用半开 bin 区间等一致规则，避免边界重复计数。零长度接触不贡献驻留时间。对速度范围外的驻留时间单独记录，不能默默丢弃后声称恢复了总密度。线程并行时使用线程独立 histogram 后归并，避免数据竞争。

## 4. 稳态、截断和误差检查

- 第 1 节假设稳恒注入、静态场和独立测试粒子，且 relevant residence integrals 收敛。所有飞行年龄从 0 开始是正确的，不要求随机设置释放时刻。
- 最大飞行年龄是积分截断参数，不是观测时间。增加它直到探头 PSD 收敛；达到时限应记录 `time_limit`，不能视为已逃逸或已损失。
- 持续注入且存在永久束缚、无损失的粒子群时，密度可能持续增长，有限稳态未必存在。此时要定义源开启时间或有依据的损失过程。
- 对时变源或时变场，需要释放时刻及观测时刻。静态场中从时刻 0 开启恒定源，在时刻 T 的结果可将飞行年龄积分截到 T；任意时变场一般必须显式采样释放时刻，不能套用一条稳态轨迹。
- 若有物理损失概率、碰撞分支或电荷态变化，使用对应的存活权重和物种轨迹段。不能无依据添加这些过程。
- 收敛检查分别改变粒子数、探头尺寸、速度 bin 宽度、积分步长和最大飞行年龄。相同轨迹的多个驻留子段不是独立样本，误差估计应按独立源样本或独立随机种子批次进行。

两个必要的一致性检验：

$$
\sum_b f_{D,b}\Delta^3v_b
=\frac{1}{V_D}\sum_iQ_i\tau_{i,D}
$$

前提是速度 bins 覆盖所有探头内速度。另一个解析案例是速度 u 的均匀束流沿立方体法向穿越：通过面积 $L^2$ 的粒子率 $Q=n u L^2$，驻留时间 $L/u$，因此 $Q(L/u)/L^3=n$。只数进入事件而不计时间，会把快粒子过度计权，得到通量性质的分布。

## 5. 与通量和二维图的区别

穿过探头某一面的事件求和除以该面面积，得到方向选择下的通量，单位 m⁻² s⁻¹，不是密度。稳态局部分布满足 $\Gamma_+=\int_{v_n>0}v_n f\,d^3v$。不要把立方体六个面的所有穿越当作独立粒子来恢复密度。密度应统计全部速度方向，包括回流。

若绘制积分掉 $v_y$ 的二维分布，使用

$$
f_{xz}(v_x,v_z)=\sum_{b_y} f_{xyz}(v_x,v_{y,b_y},v_z)\Delta v_{y,b_y},
$$

单位为 s² m⁻⁵。固定 $v_y$ 的二维切片仍是三维 PSD 的切片，单位 s³ m⁻⁶，两者不能混称。

## 6. 依据与代码位置

上述粒子率驻留时间公式由稳态粒子数守恒推导。标准轨迹长度估计器采用权重乘体内轨迹长度并除以体积，见 [Geant4 官方 cell flux 说明](https://geant4.web.cern.ch/documentation/pipelines/master/bfad_html/ForApplicationDevelopers/Detector/hit.html)。本文统计的是数密度，沿路径使用 $dt=dl/|v|$，并加入速度 bin 指示函数；不能把 Geant4 的轨迹长度通量直接当作 PSD。

- [源权重实现](../../../src/tracing/monte_carlo_weight.jl)
- [探头 PSD 实现](../../../src/tracing/detector_psd_forward.jl)
- [源权重和单位说明](../monte_carlo_weights.md)
- [MHD 源参数与 forward tracing 接口示例](../maxwellian_source.jl)
