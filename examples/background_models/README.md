# 背景模型：MHD、GITM 与 AMPS

本目录整理模拟所用电磁场、离子矩及中性大气的背景图和绘图代码。原始输入统一保留在仓库 `data/`；本目录不复制大型场文件。

## MHD O₂⁺ 通量与温度

![MHD 通量和温度](images/flux_Ti_200_400km.png)

200、400 km 的 MSO 经纬度分布。左列为 `n*norm(Ui)`，右列为离子温度 Ti；源自 MHD 的密度、三维体速度和温度。南极异常采样点只在图中遮罩，不修复原始矩数据。详细方法及依赖见 [MHD 图说明](ionosphere_maps.md)。

```sh
python examples/background_models/plot_ionosphere_maps.py
```

绘图入口调用同目录 [sample_ionosphere_maps.jl](sample_ionosphere_maps.jl)。也可提供已有采样目录重绘：

```sh
python examples/background_models/plot_ionosphere_maps.py outputs/ionosphere_maps_<timestamp>
```

## GITM 与 AMPS

![GITM 与 AMPS](images/atmosphere_gitm200km_amps500km.png)

包含 GITM 中性密度、温度剖面及 200 km 的分布图，以及 AMPS 热氧剖面和 500 km 分布图。AMPS 不提供温度；GITM 高度延伸的假设见 [大气说明](atmosphere.md)。

```sh
python examples/background_models/plot_atmosphere.py
```

图保存在 `images/`。本地生成的 `atmosphere_source_data.npz` 和 `metadata.json` 保留在本目录，不作为原始输入上传。PNG 为已有输出，不因本次目录整理而重算物理数据。

## 模型输入与运行环境

- MHD：`data/mars_fields_spherical_from_dat.vts`，提供场和离子矩；大文件未随仓库分发。
- GITM：`data/gitm_sph.mat`，提供中性大气参数。
- AMPS：`data/amps_sph.mat`，提供热氧密度。
- Python：NumPy、Matplotlib、SciPy（大气图）及 Arial 字体。
- Julia：使用仓库 Project.toml/Manifest.toml，并使 julia 在 PATH 中可用。

所有命令从仓库根目录运行。各入口不自动安装依赖。返回 [示例总目录](../README.md)。
