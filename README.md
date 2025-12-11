# CLSM2LGD



# **1\. 项目简介**

本项目基于Python开发，旨在利用 **GLDAS Catchment Land Surface Model (CLSM)** 的水文状态变量，正演计算GRACE/GRACE-FO卫星在轨道高度处的视线向重力差（Line-of-Sight Gravity Difference, LGD）。

**核心特性：**

* **多变量支持**：支持灵活组合CLSM模型中的地下水、土壤水、积雪等多个变量，实现分量重力信号的解耦或总水储量（TWS）的重构。  
* **CLSM适配**：针对GLDAS v2.2 CLSM的数据结构进行了深度适配，支持处理含GRACE同化（DA）的高精度水文数据。  
* **高精度物理模型**：内置PREM地球弹性负荷响应模型（Load Love Numbers），确保反演结果符合地球物理真实情况。  
* **沿轨分析**：直接在卫星轨道坐标系下计算重力响应，避免了传统球谐系数截断带来的信号损失。

本项目是在项目[GLDAS2LGD](https://github.com/SingyuTang/GLDAS2LGD)基础上进行修改开发，这里只对项目如何安装运行进行介绍，更多内容见[博客-基于GLDAS_CLSM_D数据计算LGD](https://singyutang.github.io/2025/12/11/%E5%9F%BA%E4%BA%8EGLDAS-NOAH-025-3H%E5%9C%9F%E5%A3%A4%E6%B9%BF%E5%BA%A6%E6%95%B0%E6%8D%AE%E5%92%8CGLDAS-CLSM-D%E6%95%B0%E6%8D%AE%E8%AE%A1%E7%AE%97LGD%E4%BA%8C/)。

# **2\. 环境依赖与安装**

本项目的环境配置与基础LGD工具包一致（[GLDAS2LGD](https://github.com/SingyuTang/GLDAS2LGD)）。核心依赖如下：

```Bash
pip install numpy scipy netCDF4 matplotlib
```

* **NumPy**: 高性能矩阵运算。  

* **SciPy**: 勒让德多项式计算 (scipy.special.lpn)。  

* **NetCDF4**: 读取GLDAS.nc4 数据文件。  

* **Matplotlib**: 结果可视化。

* 自定义依赖 (需确保主要脚本同级目录下存在以下模块): 

  S02compute_grace_lgd: 用于加载卫星轨道 (OrbitLoader)。 

  S05plot_lgd_ra_cwt_filter: 用于筛选特定区域的轨道 (filter_complete_tracks_passing_region)。 -- read_love_numbers: 用于读取负荷勒夫数。

依赖S02compute_grace_lgd和S05plot_lgd_ra_cwt_filter可到[仓库](https://github.com/SingyuTang/grace_intersatellite_ranging_lgd_processor/tree/master/Scripts_grace_lri1b_lgd_20251111)进行下载。

请确保项目目录下包含以下核心脚本：

* lgd\_processor.py: 主程序，负责批处理流程控制。  
* clsm2lgd.py: 核心算法类，包含GLDAS读取与重力积分逻辑。  
* lgd\_plot.py: 绘图工具。  
* read\_love\_numbers.py: 负荷勒夫数读取工具。  
* data/: 存放勒夫数数据文件（如 PREM-LLNs-truncated.dat）。

# **3\. GLDAS CLSM 数据准备**

本项目专为 **GLDAS Catchment Land Surface Model L4 daily 0.25 x 0.25 degree GRACE-DA1 V2.2 (GLDAS_CLSM025_DA1_D)** 数据设计（0.25°分辨率）。

1. **数据下载**：  
   * 访问 [NASA GES DISC 网站](https://disc.gsfc.nasa.gov/datasets/GLDAS_CLSM025_DA1_D_2.2/summary?keywords=CLSM)。  
   * 文件格式示例: GLDAS_CLSM025_DA1_D.A20200504.022.nc4。  
2. 关键变量说明：  
   与Noah模型不同，CLSM包含独特的地下水变量。常用变量如下：  
   * TWS\_tavg: 总陆地水储量（已包含地下水、土壤水、雪等）。  
   * GWS\_tavg: 地下水储量（Ground Water Storage）。  
   * SoilMoist\_RZ\_tavg: 根区土壤水。  
   * SoilMoist\_S\_tavg: 表层土壤水。  
   * SoilMoist_P_tavg：整个土柱剖面水。代表代表整个剖面（从地表到基岩，可能深达数米）的土壤水。它已经包含了 Root Zone 和 Surface 的水分。
   * SWE\_tavg: 雪水当量。  
   * CanopInt\_tavg: 冠层截留水。



# 4.GRACE-FO轨道数据

程序需要 GRACE-FO 的精密轨道数据（GNV1B）来确定卫星在特定时间的位置。代码中通过 OrbitLoader 类加载 GROOPS 工作区数据，你需要根据实际情况修改 `load_orbit` 函数以适配你的轨道文件格式（如 GNV1B等）。

如果按照博客[详细介绍利用GRACE1B多日数据计算LGD工作流程二_基于LRI1B多日数据](https://singyutang.github.io/2025/11/10/详细介绍利用GRACE1B多日数据计算LGD工作流程二-基于LRI1B多日数据/)已经创建了工作根目录`workdir`并进行了相关处理，只需要将本项目文件全部拷贝到该目录下即可运行，因为本项目中的 `load_orbit` 函数默认读取`workdir/gracefo_dataset`路径下的GNV1B轨道文件。

# **5\. 配置与运行指南**

## **5.1 配置文件修改 (lgd\_processor.py)**

打开 l`gd\_processor.py`，找到 `CONFIG` 字典。这是实现 **数据更换** 与 **多变量反演** 的关键配置区。

```python
CONFIG = {
        # 路径设置
        'gldas_dir': r"I:\LGD\GLDAS_CLSM025_D_2.2_2020",
        'groops_workspace': r'G:\GROOPS\PNAS2020Workspace',
        'output_dir': r"./results_lgd/TWS_tavg",

        # 【核心配置】定义需要读取并累加的变量列表
        # 示例 1: 仅地下水 -> ['GWS_tavg']，Ground Water Storage，    # （推荐仅地下水）
        # 示例 2: 仅地表水 -> ['SoilMoist_S_tavg']，Surface Soil moisture,仅代表地表表层（通常是 0-2 cm）的土壤水。
        # 示例 3: 仅根区水 -> ['SoilMoist_RZ_tavg']，Root Zone Soil moisture，代表根区（通常是 0-100 cm）的土壤水。   （推荐仅土壤水）
        # 示例 4: 整个剖面水 -> ['SoilMoist_P_tavg']，Profile Soil moisture，代表代表整个剖面（从地表到基岩，可能深达数米）的土壤水。它已经包含了 Root Zone 和 Surface 的水分。
        # 示例 5: 地下水 + 土壤水 -> ['GWS_tavg', 'SoilMoist_RZ_tavg']
        # 示例 6: 陆地水储量 -> ['TWS_tavg']，Terrestrial water storage

        # * 建议：1. 分离使用： 单独使用土壤水使用RZ；单独使用地下水使用GWS
        #        2. 联合使用：CLSM 的 Profile 变量包含了整个土柱的水。在某些情况下，它可以用作总地下储量（地表以下，不是GWS）的近似，但它不如 RZ + GWS 的组合精细。
        'variable_names': ['TWS_tavg'],

        # 背景场设置 (用于计算距平)
        'bgd_start': '2020-05-01',
        'bgd_end': '2020-05-31',

        # 区域设置
        'region_lon': (88, 92),
        'region_lat': (22, 26),
        'track_lat_limit': (-80.0, 80.0),

        # 轨道与计算参数
        'orbit_direction': 'asc',  # 'asc' or 'desc'
        'orbit_interval': 5,  # 秒
        'cutoff_deg': 20.0,  # 积分截断距离
        'n_max': 200,  # 阶数

    }
```

修改待处理日期列表`date_list`

```python
date_list = [
    '2020-06-04', '2020-06-10', '2020-06-15', '2020-06-21', '2020-06-26',
    '2020-07-02', '2020-07-07', '2020-07-13', '2020-07-18', '2020-07-24', '2020-07-29',
    '2020-08-04', '2020-08-09', '2020-08-15', '2020-08-20'
]
```



## **5.2 运行主程序**

配置完成后，直接运行处理脚本：

```bash
python lgd_processor.py
```

**程序运行逻辑：**

1. **背景场预计算**：程序首先读取 bgd\_start 到 bgd\_end 期间的所有GLDAS文件，提取 variable\_names 中定义的所有变量并求和，计算其时间平均值。  
2. **每日反演**：  
   * 读取当日GLDAS文件，累加目标变量。  
   * 减去背景场均值，得到**质量异常 (Mass Anomaly)**。  
   * 结合当日卫星轨道位置，利用积分核函数计算LGD。  
3. **结果存储**：结果将以压缩的 .npz 格式保存在 output\_dir 中，包含 lat, lon, lgd, alt 等数组。

## **5.3 结果可视化**

使用 lgd\_plot.py 绘制沿轨LGD剖面堆叠图：

```Bash
python lgd_plot.py
```

*注意：请确保 lgd\_plot.py 中的 results\_dir 指向了正确的输出目录。*

# 6.参考文章



Ghobadi-Far, Khosro, Shin-Chan Han, Christopher M. McCullough, David N. Wiese, Richard D. Ray, Jeanne Sauber, Linus Shihora, and Henryk Dobslaw. 2022. “Along-Orbit Analysis of GRACE Follow-On Inter-Satellite Laser Ranging Measurements for Sub-Monthly Surface Mass Variations.” Journal of Geophysical Research: Solid Earth 127(2):e2021JB022983. doi:10.1029/2021JB022983.

[博客-基于GLDAS_CLSM_D数据计算LGD](https://singyutang.github.io/2025/12/11/%E5%9F%BA%E4%BA%8EGLDAS-NOAH-025-3H%E5%9C%9F%E5%A3%A4%E6%B9%BF%E5%BA%A6%E6%95%B0%E6%8D%AE%E5%92%8CGLDAS-CLSM-D%E6%95%B0%E6%8D%AE%E8%AE%A1%E7%AE%97LGD%E4%BA%8C/)
