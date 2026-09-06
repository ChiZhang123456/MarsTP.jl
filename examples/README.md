# MarsTP 示例目录

| 目录 | 内容 |
|---|---|
| [background_models](background_models/README.md) | MHD、GITM、AMPS 背景参数，通量、温度及大气分布图 |
| [forward_tracing](forward_tracing/README.md) | 800 km 球面随机采样，前向轨迹和对流/Hall/总电场累计做功 |
| [backward_tracing](backward_tracing/README.md) | 固定探头反向轨迹、薄层面源 PSD 推导及初步二维 VDF |

每个子目录包含说明、绘图/积分代码和 images；反向追踪目录还包含初步 VDF 的小型数据及检查结果。原始物理输入仍位于仓库 data，完整运行结果保留于 outputs。所有运行命令均从仓库根目录执行，具体依赖见各目录 README。

最新二维 VDF 尚未通过速度网格与时间步长收敛检查，保留为明确标注的初步示例。
