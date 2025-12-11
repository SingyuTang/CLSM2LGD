import time
from S05plot_lgd_ra_cwt_filter import filter_complete_tracks_passing_region
from clsm2lgd import *
import numpy as np
from datetime import datetime
import os

def get_gldas_filename(root_dir, date_str):
    """根据日期生成GLDAS文件名"""
    # 格式: GLDAS_CLSM025_DA1_D.A20200501.022.nc4
    dt = datetime.strptime(date_str, "%Y-%m-%d")
    date_part = dt.strftime("A%Y%m%d")
    filename = f"GLDAS_CLSM025_DA1_D.{date_part}.022.nc4"
    return os.path.join(root_dir, filename)


def save_results(output_dir, date_str, lat, lon, lgd, alt):
    """保存计算结果到文件"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    save_path = os.path.join(output_dir, f"LGD_Result_{date_str}.npz")
    np.savez(save_path, lat=lat, lon=lon, lgd=lgd, alt=alt)
    print(f"结果已保存至: {save_path}")


def process_single_date(calc, date_str, bg_mean_grid, config):
    """
    处理单个日期的核心逻辑
    :param calc: TrackLGDCalculator 实例
    :param date_str: 处理日期 'YYYY-MM-DD'
    :param bg_mean_grid: 预计算好的背景场累加均值网格 (2D array, 单位 mm)
    :param config: 配置字典
    :return: (lats, lgd_values) 或 (None, None) 如果失败
    """
    print(f"\n[{date_str}] 开始处理...")

    # 1. 读取当日 GLDAS 数据
    gldas_file = get_gldas_filename(config['gldas_dir'], date_str)
    if not os.path.exists(gldas_file):
        print(f"[{date_str}] 错误: GLDAS文件不存在 -> {gldas_file}")
        return None, None

    # 获取配置中的变量列表
    var_list = config.get('variable_names', ['GWS_tavg'])

    try:
        # 【关键步骤】读取单天数据，函数内部会自动将列表中所有变量累加
        ds = calc.read_clsm_var(gldas_file, variable_names=var_list)
    except Exception as e:
        print(f"[{date_str}] 读取GLDAS出错: {e}")
        return None, None

    # 准备网格数据
    gldas_lat_range, gldas_lon_range = np.array(ds['lat']), np.array(ds['lon'])
    lon_grid, lat_grid = np.meshgrid(gldas_lon_range, gldas_lat_range)

    # 【关键步骤】构建键名提取累加后的数据
    # clsm2lgd.py 中的逻辑是 key = "Var1+Var2+..."
    daily_data_key = "+".join(var_list)

    if daily_data_key not in ds:
        print(f"[{date_str}] 错误: 返回的数据字典中找不到键 '{daily_data_key}'")
        print(f"可用键: {list(ds.keys())}")
        return None, None

    # 获取当日累加网格 (例如: GWS + SoilMoi)
    curr_daily_sum_grid = np.squeeze(ds[daily_data_key])

    # 确保维度一致
    if curr_daily_sum_grid.shape != bg_mean_grid.shape:
        print(f"[{date_str}] 维度不匹配: 当日数据 {curr_daily_sum_grid.shape} vs 背景场 {bg_mean_grid.shape}")
        return None, None

    # 【关键步骤】计算距平 (当日累加总和 - 背景累加总和均值)
    # 单位转换: mm -> m
    h_grid = (curr_daily_sum_grid - bg_mean_grid) * 1e-3

    # 展平用于积分计算
    g_h = h_grid.flatten()
    g_lat = lat_grid.flatten()
    g_lon = lon_grid.flatten()

    # 2. 加载轨道数据
    try:
        # 注意：这里加载的是当天的轨道
        orbitc, orbitd = load_orbit(date_str=date_str,
                                    groops_workspace=config['groops_workspace'],
                                    coord_type='geodetic')
    except Exception as e:
        print(f"[{date_str}] 加载轨道数据出错: {e}")
        return None, None

    # 采样与提取
    interval = config['orbit_interval']
    posc_geo = np.array([obj.position for obj in orbitc])[::interval]
    posd_geo = np.array([obj.position for obj in orbitd])[::interval]

    # 3. 筛选经过研究区的轨迹
    lonlat = np.column_stack([posc_geo[:, 0], posd_geo[:, 1]])

    tracks, indices = filter_complete_tracks_passing_region(
        lonlat,
        config['region_lon'],
        config['region_lat'],
        lat_limit=config['track_lat_limit'],
        separate=False,
        direction=config['orbit_direction']
    )

    if indices is None or len(indices) == 0:
        print(f"[{date_str}] 未找到经过目标区域的有效轨迹。")
        return None, None

        # 提取筛选后的轨迹点
    track_s1 = posc_geo[np.squeeze(indices)]
    track_s2 = posd_geo[np.squeeze(indices)]

    # 更新 calc 对象中的轨迹属性
    calc.track_s1_lon, calc.track_s1_lat, calc.track_s1_h = track_s1[:, 0], track_s1[:, 1], track_s1[:, 2]
    calc.track_s2_lon, calc.track_s2_lat, calc.track_s2_h = track_s2[:, 0], track_s2[:, 1], track_s2[:, 2]

    # 4. 执行 LGD 计算
    t_start = time.time()
    lgd_series = calc.compute_track_lgd(
        calc.track_s1_lat,
        calc.track_s1_lon,
        g_h, g_lat, g_lon,
        cutoff_deg=config['cutoff_deg']
    )
    print(f"[{date_str}] 计算完成，耗时 {time.time() - t_start:.2f} 秒")

    # 保存结果
    save_results(config['output_dir'], date_str, calc.track_s1_lat, calc.track_s1_lon, lgd_series, calc.track_s1_h)

    return calc.track_s1_lat, lgd_series


def main():
    # ================= 配置区域 =================
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

    # ================= 待处理日期列表 =================
    date_list = [
        '2020-06-04', '2020-06-10', '2020-06-15', '2020-06-21', '2020-06-26',
        '2020-07-02', '2020-07-07', '2020-07-13', '2020-07-18', '2020-07-24', '2020-07-29',
        '2020-08-04', '2020-08-09', '2020-08-15', '2020-08-20'
    ]

    print(f"计划处理日期: {date_list}")

    # ================= 初始化 =================
    calc = TrackLGDCalculator(n_max=CONFIG['n_max'])

    # 1. 预计算背景场 (只做一次)
    print(f"\n>>> 正在计算背景场 ({CONFIG['bgd_start']} 到 {CONFIG['bgd_end']})...")

    # 从配置中获取变量列表
    var_list = CONFIG.get('variable_names', ['GWS_tavg'])
    bg_mean_grid = None

    print(f"计划处理变量: {var_list}")

    try:
        # 【关键步骤】传入变量列表，计算这些变量累加后的时间平均场
        bgd_data = calc.calculate_vars_mean_by_date_range(
            CONFIG['gldas_dir'],
            CONFIG['bgd_start'],
            CONFIG['bgd_end'],
            variable_names=var_list
        )

        # 【关键步骤】构建键名提取背景数据
        # clsm2lgd.py 中的逻辑是 key = "mean_" + "Var1+Var2+..."
        bg_data_key = "mean_" + "+".join(var_list)

        if bg_data_key not in bgd_data:
            raise KeyError(f"在背景场结果中找不到键 '{bg_data_key}'")

        bg_mean_grid = np.squeeze(bgd_data[bg_data_key])
        print(f"背景场计算完成。")
        print(f"使用了以下变量进行累加: {var_list}")
        print(f"背景场网格形状: {bg_mean_grid.shape}, 均值: {np.nanmean(bg_mean_grid):.2f} mm")

    except Exception as e:
        print(f"背景场计算失败，程序终止: {e}")
        return

    # ================= 批量循环 =================
    success_count = 0
    results_summary = {}

    for date_str in date_list:
        # 传入背景场网格，process_single_date 内部会使用同样的变量列表读取单天数据
        lats, lgds = process_single_date(calc, date_str, bg_mean_grid, CONFIG)
        if lats is not None:
            success_count += 1
            results_summary[date_str] = (lats, lgds)

    print(f"\n==========================================")
    print(f"批量处理结束。成功: {success_count}/{len(date_list)}")
    print(f"结果已保存至: {CONFIG['output_dir']}")


if __name__ == "__main__":
    main()