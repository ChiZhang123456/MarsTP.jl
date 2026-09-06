# O₂⁺ backtracing：体积源与电离层薄层面源的推导和单位

本文对应 `src/tracing/detector_psd_backward.jl` 的 `thin_shell_source_v1` 模型。它取代此前“到达电离层高度即施加边界分布并终止”的模型。电离层现在是一层可穿越的面源，**默认离地 400 km；吸收内边界始终位于离地 200 km**。

## 1. 模型、位置与粒子定义

使用静态 MHD 电磁场，积分非相对论 O₂⁺ 的洛伦兹运动。质量和正电荷取自项目实际 TestParticle `SpeciesDict["O2+"]`，不手写近似质量。位置、速度均为笛卡尔分量，MHD 球网格使用火心距离、余纬和方位角。

设火星半径 $R_M=3.390\times10^6\ \mathrm{m}$，火心距离 $r=|\mathbf x|$，离地高度 $h=r-R_M$。积分域为

$$
r_{\rm in}=R_M+200\times10^3\ \mathrm{m},\qquad r_{\rm out}=4R_M.
$$

电离层薄层半径为

$$
r_s=R_M+h_s,\qquad h_s=400\times10^3\ \mathrm{m}\quad\text{（默认）}.
$$

`ionosphere_altitude_km` 仍允许 200 至 800 km；改变它不改变内边界。体积源在整个有效积分域内存在，包括 200 至 400 km；面源仅在 $r=r_s$ 处存在。300 km 处的局地产生项只有体积源，但是该处的累计 PSD 仍可能包含轨迹曾穿越 400 km 时产生的粒子。

当前模型不包含碰撞、损失、重力、化学反馈、粒子自洽电磁场或时间变化源。没有将各高度的出流率再次映射到 400 km。

## 2. 相空间密度的定义与单位

这里使用速度空间分布 $f$，不是动量空间分布，也不是能量谱：

$$
dN=f(\mathbf x,\mathbf v,t)\,d^3x\,d^3v,\qquad
n(\mathbf x,t)=\int f(\mathbf x,\mathbf v,t)\,d^3v.
$$

因此

$$
[f]=\mathrm{m^{-3}(m\,s^{-1})^{-3}}=\boxed{\mathrm{s^3\,m^{-6}}}.
$$

粒子个数在量纲分析中视为无量纲。若文献用动量分布 $f_p$，非相对论条件下 $\mathbf p=m_i\mathbf v$，则 $f_v=m_i^3f_p$；二者不能直接比较数值。

| 量 | 定义 | 内部 SI 单位 |
|---|---|---|
| $\mathbf x,r,R_M$ | 位置或火心距离 | m |
| $t,\Delta t$ | 物理时间与积分步长 | s |
| $\mathbf v$ | 被追踪粒子的瞬时速度 | m s⁻¹ |
| $\mathbf U_i$ | MHD O₂⁺ 三维体速度 | m s⁻¹ |
| $\mathbf E,\mathbf B$ | 电场与磁场 | V m⁻¹、T |
| $n_i$ | MHD O₂⁺ 数密度 | m⁻³ |
| $T_i,T_n$ | 离子温度与中性温度 | K |
| $m_i,q_i$ | 粒子质量与电荷 | kg、C |
| $Q_V$ | 单位体积、单位时间产生的粒子数 | m⁻³ s⁻¹ |
| $F_s$ | 单位面积、单位时间产生的粒子数 | m⁻² s⁻¹ |
| $g_V,g_M$ | 归一化三维速度概率密度 | s³ m⁻³ |
| $\delta(r-r_s)$ | 径向 Dirac delta | m⁻¹ |
| $S_V,S_s$ | 相空间产生项，即 $df/dt$ | s² m⁻⁶ |
| $f,f_V,f_s$ | 三维速度空间 PSD | s³ m⁻⁶ |
| $f_{xz}$ | 对 $v_y$ 积分后的二维 VDF | s² m⁻⁵ |

`velocity_axes` 将 km/s 乘 $10^3$ 转成 m/s，再进行轨迹和源项计算。图中的 cm⁻² s⁻¹ 通量由 SI 通量除以 $10^4$ 得到，不直接输入积分器。

## 3. 从输运方程到轨迹积分

有源无损的相空间连续方程为

$$
\frac{\partial f}{\partial t}
+\nabla_{\mathbf x}\cdot(\mathbf v f)
+\nabla_{\mathbf v}\cdot(\mathbf a f)=S_V+S_s,
\qquad
\mathbf a=\frac{q_i}{m_i}(\mathbf E+\mathbf v\times\mathbf B).
$$

洛伦兹运动满足 $\nabla_{\mathbf x}\cdot\mathbf v=0$ 与 $\nabla_{\mathbf v}\cdot\mathbf a=0$，因此沿同一条特征轨迹

$$
\dot{\mathbf x}=\mathbf v,\qquad
\dot{\mathbf v}=\mathbf a,\qquad
\frac{df}{dt}=S_V[\mathbf x(t),\mathbf v(t)]+S_s[\mathbf x(t),\mathbf v(t)].
$$

取探头观测时刻 $t_0=0$，从指定 $(\mathbf x_0,\mathbf v_0)$ 反向积分到 $t_*<0$：

$$
f(\mathbf x_0,\mathbf v_0,0)
=f[\mathbf x(t_*),\mathbf v(t_*),t_*]
+\int_{t_*}^{0}(S_V+S_s)\,dt.
$$

当前代码设置过去端点的外加背景分布为零，所以返回的是本次回溯窗口内的源项贡献。到达 200 km 或 $4R_M$ 时终止；到达回溯时限时也终止，但该结果是截断时间窗内的累计值，不自动等于稳态解。现有代码不提供外边界注入的非零背景 VDF。

实现使用负时间步长，质量、电荷及物理速度符号均不变。定义回溯时长 $\tau=-t\ge0$，源积分也可写为 $\int_0^{\tau_*}S[\mathbf x(-\tau),\mathbf v(-\tau)]d\tau$。这解释了代码累计源项时使用 $|\Delta t|$，无需把正源项取负。

## 4. 体积产生项

现有 MAT 数据提供

$$
Q_V(\mathbf x)\quad[\mathrm{m^{-3}s^{-1}}].
$$

该量本身不是 PSD。新生粒子的归一化速度分布为

$$
g_V(\mathbf v;\mathbf x)
=\frac{1}{\pi^{3/2}w_n^3}
\exp\left(-\frac{|\mathbf v|^2}{w_n^2}\right),\qquad
w_n=\sqrt{\frac{2k_BT_n(\mathbf x)}{m_i}},\qquad
\int g_Vd^3v=1.
$$

这里保留原模型：$T_n$ 来自 GITM，体积产生粒子的漂移速度为零。这与薄层使用 MHD 的 $T_i,\mathbf U_i$ 是两个不同的出生分布，不能混淆。

$$
S_V=Q_Vg_V,
\qquad [S_V]=(\mathrm{m^{-3}s^{-1}})(\mathrm{s^3m^{-3}})=\mathrm{s^2m^{-6}},
$$

$$
\boxed{f_V=\int_{t_*}^{0}Q_V[\mathbf x(t)]\,g_V[\mathbf v(t);\mathbf x(t)]\,dt}.
$$

$[f_V]=\mathrm{s^3m^{-6}}$。速度参数是每个轨迹点的瞬时速度，不是始终使用探头的初始速度。积分变量是时间，不是距离，不再乘速度或网格体积。

## 5. MHD 薄层面产生率与速度分布

在 $r_s$ 球面上插值读取 MHD O₂⁺：

- `n_O^2^p [m^-3]`，得到 $n_i$；
- `T_O^2^p [K]`，得到 $T_i$；
- `U_O^2^p [m/s]`，得到笛卡尔向量 $\mathbf U_i$。

按本项目约定，面产生率为

$$
\boxed{F_s(\Omega)=n_i(r_s,\Omega)|\mathbf U_i(r_s,\Omega)|},\qquad
[F_s]=\mathrm{m^{-2}s^{-1}}.
$$

这是一项模型处方，不是有符号的球面法向流量 $n_i\mathbf U_i\cdot\hat{\mathbf r}$，也不是麦氏分布积分得到的单向热通量。它把 MHD 局地矩转换成独立面源强度，不能据此声称 MHD 本身给出了独立的粒子出生率。

归一化漂移麦氏分布为

$$
g_M(\mathbf v;\Omega)
=\frac{1}{\pi^{3/2}w_i^3}
\exp\left[-\frac{|\mathbf v-\mathbf U_i|^2}{w_i^2}\right],\qquad
w_i=\sqrt{\frac{2k_BT_i}{m_i}},\qquad \int g_Md^3v=1.
$$

此处 $g_M$ **不含密度因子**；密度已经进入 $F_s$。`ionosphere_distribution` 同时返回 `g=g_M` 与 `f=n_i*g_M`，面源计算只使用 `flux*g`，不使用 `flux*f`，从而避免重复乘密度。

当前使用完整三维麦氏分布，不裁剪速度半空间。因此向内、向外两类穿越均可贡献；若以后要求只发射向外粒子，需要重新定义并归一化速度分布，不能简单删除一半积分。

## 6. 为什么薄层贡献要除以粒子的径向速度

把面源写成三维空间中的体源：

$$
S_s(\mathbf x,\mathbf v)=F_s(\Omega)g_M(\mathbf v;\Omega)\delta(r-r_s).
$$

单位为

$$
[S_s]=(\mathrm{m^{-2}s^{-1}})(\mathrm{s^3m^{-3}})(\mathrm{m^{-1}})
=\mathrm{s^2m^{-6}}.
$$

对于球面有 $|\nabla(r-r_s)|=1$，且

$$
\int F_s\delta(r-r_s)d^3x=\int F_s r_s^2d\Omega=\int F_s dA.
$$

因此源项中不需要再放一个 $r_s^2$ 或面元面积；几何面积已经由空间体积元与 delta 函数体现。

设轨迹在 $t_j$ 横穿薄层，即 $r(t_j)=r_s$ 且 $\dot r(t_j)\ne0$。利用一维 delta 的变量代换：

$$
\delta[r(t)-r_s]
=\sum_j\frac{\delta(t-t_j)}{|\dot r(t_j)|},\qquad
\dot r(t_j)=\mathbf v(t_j)\cdot\hat{\mathbf r}_j.
$$

于是

$$
\boxed{
f_s=\sum_j\frac{F_s(\Omega_j)g_M[\mathbf v(t_j);\Omega_j]}
{|\mathbf v(t_j)\cdot\hat{\mathbf r}_j|}
}.
$$

每次穿越贡献的单位为

$$
[\Delta f_s]
=\frac{(\mathrm{m^{-2}s^{-1}})(\mathrm{s^3m^{-3}})}{\mathrm{m\,s^{-1}}}
=\boxed{\mathrm{s^3m^{-6}}}.
$$

分母使用的是**测试粒子的瞬时径向速度**，不是 MHD 的径向体速度。这一因子来自粒子穿过薄层的停留时间几何关系，并不改变 $F_s=n_i|\mathbf U_i|$ 的定义。

如果同一条轨迹再次穿越该层，每次横穿都独立累计，不覆盖之前的结果。跨步共享端点仅算一次。源层不会改变粒子速度或电荷，也不终止轨迹。

## 7. 总 PSD 与代码实际输出

零外加背景条件下，探头处每个速度点的三维 PSD 为

$$
\boxed{
f(\mathbf x_0,\mathbf v_0)
=\int_{t_*}^{0}Q_V[\mathbf x(t)]g_V[\mathbf v(t);\mathbf x(t)]dt
+\sum_j\frac{n_i(\Omega_j)|\mathbf U_i(\Omega_j)|g_M[\mathbf v(t_j);\Omega_j]}
{|\mathbf v(t_j)\cdot\hat{\mathbf r}_j|}
}.
$$

随后代码对探头速度网格的 $v_y$ 做矩形求和：

$$
f_{xz}(v_x,v_z)=\int f(v_x,v_y,v_z)dv_y
\simeq\sum_{\ell}f(v_x,v_{y,\ell},v_z)\Delta v_y.
$$

$$
[f_{xz}]=(\mathrm{s^3m^{-6}})(\mathrm{m\,s^{-1}})
=\boxed{\mathrm{s^2m^{-5}}}.
$$

| 输出字段 | 内容 | 单位 |
|---|---|---|
| `f2d_volume` | $\sum_\ell f_V\Delta v_y$ | s² m⁻⁵ |
| `f2d_ionosphere` | $\sum_\ell f_s\Delta v_y$ | s² m⁻⁵ |
| `f2d_xz` | 上述两项之和 | s² m⁻⁵ |
| `ionosphere_crossings` | 每个 $(v_x,v_z)$ 格点所有已计算 Vy 轨迹的穿层次数之和 | 无量纲计数，不乘 $\Delta v_y$ |
| `vx_km,vy_km,vz_km` | 探头速度网格 | km/s |
| `model` | `thin_shell_source_v1` | 模型标识 |
| `source_units` | 内部各类源项和 PSD 的单位字典 | 字符串 |

即使只设置一个 Vy 网格点，代码也会乘 `dvy_kms*1000`，此时输出表示该速度宽度的矩形近似，**不是**三维 PSD 的切片。恢复数密度还需要对 $v_x,v_z$ 做积分，并验证速度网格覆盖范围和分辨率。输出不是每单位 eV 的微分通量。

## 8. 离散实现与边界处理

1. 在 SI 单位下从探头状态使用负时间步长运行 Boris 或 AdaptiveBoris。
2. 逐步检查是否穿过 200 km 内边界或 $4R_M$ 外边界，截取该步的有效部分。
3. 对这段位置线段求与薄层球面的交点，按时间顺序处理所有交点。
4. 在每个交点同步 Boris 的半步速度，再评估 $F_sg_M/|v_r|$。
5. 在交点及有效段末端分段计算体积源梯形积分：

$$
\Delta f_V\simeq\frac{S_{V,k}+S_{V,k+1}}{2}|t_{k+1}-t_k|.
$$

6. 穿层后继续至 200 km、外边界或回溯时限。400 km 以下不重复添加薄层项，但仍累计体积源。
7. 分别累计三维 PSD 的两类贡献，再对 Vy 做积分并保存。

求交使用位置线段近似，交点速度由同一 Boris 步的交错速度同步得到，仍存在有限步长误差。一步内的真实曲线可能多次穿层，而线段近似无法完全解析，因此必须做步长收敛检查。代码不在场域外求值。边界附近 1e-8 m 的偏移只用于避免浮点舍入出界。

状态计数索引依次为：1 保留未使用，2 达到时限，3 到达 200 km，4 非有限轨迹，5 到达外边界。穿过 400 km 不属于终止状态。非有限轨迹的已累计部分可能仍在输出中，必须检查 `status_counts[4]==0` 后再用于完整科学分析。无效源项或 MHD 矩会报错，不替换为零。

## 9. 零厚度模型的限制

当 $|v_r|\to0$ 时，面源 PSD 可能变大；严格切向事件不满足上述 delta 代换的简单根条件。代码不会使用任意最小速度去截断分母。数值上，若 $|v_r|\le\sqrt{\epsilon_{64}}\max(|\mathbf v|,1\ \mathrm{m\,s^{-1}})$，或线段求交判为不可分辨切向事件，将报错。该阈值是拒绝无法解析事件的数值判据，不是物理正则化。

探头恰好位于零厚度源面上也会被拒绝，避免初始时刻源项的单侧归属歧义。源高为 200 km 时，按源面从域内趋近吸收边界的单侧极限约定，在终止交点计入一次完整横穿权重。

若必须研究切向轨迹，需另行指定有限层厚 $\Delta r$ 与归一化径向剖面 $W(r)$，其中 $\int Wdr=1$、$[W]=\mathrm{m^{-1}}$，再计算 $\int F_sg_MW[r(t)]dt$。当前实现未加入这个额外物理参数。

300 km 的局部面源严格为零，但从 300 km 开始反向向外运动的轨迹可能穿过 400 km，从而取得面源贡献。这是沿特征线累计源项的结果，不是把 400 km 面源涂抹到所有高度。

MHD 的某些极点包含近零但有限的矩值；分布图的南极遮罩只是绘图处理，不会自动应用到输运计算。模型还假定所选面源与 MAT 体积源可相加，未自动去除潜在的物理来源重叠。研究总体产生率时需要独立核实这项假设。

## 10. 运行与验证

```julia
using MarsTP, StaticArrays
cfg = BacktraceConfig(
    detector_Rm=SA[0.0,0.0,2.0],
    ionosphere_altitude_km=400.0,
    include_ionosphere=true,
    dt=-0.05,
    tspan=(0.0,-500.0),
    vx_min_kms=-20.0, vx_max_kms=20.0,
    vy_min_kms=-20.0, vy_max_kms=20.0,
    vz_min_kms=-20.0, vz_max_kms=20.0,
    dv_kms=10.0, dvy_kms=10.0,
)
result = run_backtrace_vdf(cfg)
```

上述速度范围仅为小规模示例，不能视为已经收敛的科学 VDF。运行前准备项目 Julia 环境及 MHD、源率和 GITM 输入。软件不会自动安装依赖。输出路径可以在配置中另行指定，以保留不同高度和步长的结果。

```sh
julia --project=. test/runtests.jl
julia --project=. scripts/smoke_ionosphere.jl
```

自动测试覆盖密度与速度 PDF 的区分、通量单位、200/400/800 km 源面与固定 200 km 内边界、400 km 以下的体积源、反向向外穿层、斜向两次穿层、共享端点去重、切向拒绝、Boris 速度同步和步长减半。真实 MHD 小规模测试另行检查穿层后继续前进及非零面源贡献。理论推导不替代时间步长、速度网格和回溯时限的收敛检验。


### 本次执行记录

Julia 1.12.6、TestParticle 0.23.3 下，项目测试共 **122 项通过**。实际 MHD 单粒子测试从 `(1.12 Rm, 0, 0)`、初始速度 `(10, 0, 0) km/s` 出发，源层为 400 km，步长 −0.05 s，回溯 2 s。轨迹穿层后继续运行至时限，未在 400 km 停止；面源的 Vy 积分贡献为约 `1.17466e-5 s^2 m^-5`，此处仅一个 Vy 点，速度宽度为 1 km/s。该数值是端到端运行检查，不是已收敛的观测预测。
