import matplotlib.pyplot as plt
from skyfield.api import EarthSatellite, load
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import timedelta  # 导入 timedelta
import numpy as np
# 加载时间尺度
ts = load.timescale()

# ################################################################
# 设置计算卫星位置的时刻
# 由于TLE文件有效期，注意跟TLE文件时间不要相差太多
t = ts.utc(2024, 7, 28, 2, 50)  # 替换为所需的时刻
# 读取TLE文件并解析卫星数据，TLE文件应包括所有想绘制的卫星
tle_file = "20240825_GRID.TLE"  # 替换为您的TLE文件路径
# ################################################################

# 打印调试信息
debug = 0

'''
# 是否显示昼夜区域
bSun = 0
'''

# 定义地球的标准引力参数（单位：km^3/s^2）
GM_earth = 398600.4418

satellites = []
satellite_names = []
with open(tle_file, 'r') as f:
    lines = f.readlines()
    for i in range(0, len(lines), 3):
        name = lines[i].strip()
        line1 = lines[i+1].strip()
        line2 = lines[i+2].strip()
        satellites.append(EarthSatellite(line1, line2, name, ts))
        satellite_names.append(name)

# 定义投影列表
projections = [ccrs.Robinson(), ccrs.Mollweide(), ccrs.EqualEarth(), ccrs.PlateCarree()]
proj_names = ["Robinson", "Mollweide", "Equal Earth", "Plate Carree"]
current_proj_idx = 0  # 初始化当前投影索引

'''
# 计算昼夜区域和晨昏线
def calculate_day_night_boundary_simple(date, lat_step=1, lon_step=1):
    # 计算格林尼治恒星时
    JD = 367 * date.year - (7 * (date.year + ((date.month + 9) // 12)) // 4) + (275 * date.month // 9) + date.day + 1721013.5 + (date.hour + date.minute / 60 + date.second / 3600) / 24
    T = (JD - 2451545.0) / 36525
    gst = (280.46061837 + 360.98564736629 * (JD - 2451545) + T**2 * (0.000387933 - T / 38710000)) % 360
    day_of_year = date.timetuple().tm_yday
    # 计算太阳赤纬角
    declination_angle = 23.44 * np.sin(np.radians((360 / 365.25) * (day_of_year - 81)))
    
    lon = np.arange(-180, 180, lon_step)
    lat = np.arange(-90, 90, lat_step)
    
    sin_declination = np.sin(np.radians(declination_angle))
    cos_declination = np.cos(np.radians(declination_angle))
    
    sun_altitude_estimate = np.zeros((len(lat), len(lon)))
    
    for i in range(len(lat)):
        sin_lat = np.sin(np.radians(lat[i]))
        cos_lat = np.cos(np.radians(lat[i]))
        
        for j in range(len(lon)):
            local_sidereal_time = gst + lon[j]
            hour_angle = (local_sidereal_time - 180) % 360  # 地方时角
            if hour_angle > 180:
                hour_angle -= 360
            
            sun_altitude_estimate[i, j] = (sin_lat * sin_declination) + \
                                          (cos_lat * cos_declination * np.cos(np.radians(hour_angle)))
    
    # 将结果从 [-1, 1] 映射至 [-90, 90] 角度，便于判断
    sun_altitude_degrees = np.degrees(np.arcsin(sun_altitude_estimate))
    
    return lon, lat, sun_altitude_degrees
'''

def update_plot(t):
    global bSun
    plt.clf()  # 清除整个图像
    ax = plt.axes(projection=projections[current_proj_idx])  # 重新创建 `ax` 并设置投影
    ax.set_global()  # 保证新的投影在全球范围内展现
    # 更淡的颜色配置（接近白色）
    ax.add_feature(cfeature.LAND.with_scale('110m'), edgecolor='none', facecolor='#F0F0F0')  # 非常浅的陆地颜色
    ax.add_feature(cfeature.OCEAN.with_scale('110m'), facecolor='#FAFAFA')  # 非常非常浅的海洋颜色
    ax.add_feature(cfeature.COASTLINE, edgecolor='#C0C0C0')
    # 不同卫星使用不同颜色
    colors = plt.cm.get_cmap('tab10', len(satellites))

    '''
    if bSun:
        # 计算白天和黑夜边界（优化方法）
        lon, lat, sun_altitude = calculate_day_night_boundary_simple(t.utc_datetime())
        # 添加阴影以表示黑夜区域
        ax.contourf(lon, lat, sun_altitude, levels=[-90, 0], colors='gray', alpha=0.1, transform=ccrs.PlateCarree())
    '''

    for idx, satellite in enumerate(satellites):
        # 获取卫星当前位置
        geocentric = satellite.at(t)
        subpoint = geocentric.subpoint()
        lat = subpoint.latitude.degrees
        lon = subpoint.longitude.degrees

        # 画出卫星当前的位置（加大五角星）并标注名称（加大字体）
        ax.plot(lon, lat, marker='*', color=colors(idx), markersize=24, transform=ccrs.Geodetic())
        if '&' in satellite_names[idx]:
            ax.text(lon + 6, lat - 6, satellite_names[idx].replace('&', '\n'), transform=ccrs.Geodetic(), fontsize=20, color="black", weight='bold')
        else:
            ax.text(lon + 6, lat, satellite_names[idx], transform=ccrs.Geodetic(), fontsize=20, color="black", weight='bold')

        # 根据开普勒第三定律计算轨道周期
        r = geocentric.distance().km  # km
        period_minutes = (2 * np.pi * np.sqrt((r**3) / GM_earth)) / 60
        # 粗略取90分钟，以免TLE文件过期计算出很奇怪的值
        # period_minutes = 90
        if np.isnan(period_minutes): # 检查，避免是NaN的值
            continue
        
        # 计算并绘制此后一轨的轨迹（大致绕地球一周）
        trajectory_lats = []
        trajectory_lons = []
        for minutes in range(0, int(period_minutes*1.08), int(period_minutes/45)+1): #2):  # 根据计算的轨道周期，动态调整分钟数
            offset_time = t + minutes / (24 * 60)  # Skyfield 中使用天文单位，1 天 = 24 * 60 分钟
            geocentric = satellite.at(offset_time)  # 使用偏移后的时间
            subpoint = geocentric.subpoint()
            if debug == 1:
                print(offset_time.utc_datetime(), geocentric.subpoint().latitude.degrees, geocentric.subpoint().longitude.degrees)
            trajectory_lats.append(subpoint.latitude.degrees)
            trajectory_lons.append(subpoint.longitude.degrees)
        ax.plot(trajectory_lons, trajectory_lats, '-', color=colors(idx), transform=ccrs.Geodetic())
    # 动态生成标题中的时间
    if debug == 1:
        title_str = f'UTC {t.utc_iso().replace("T", " ")[:19]} + 1 orbit \nProjection: {proj_names[current_proj_idx]}'
    else:
        title_str = f'UTC {t.utc_iso().replace("T", " ")[:19]} + 1 orbit'
    ax.set_title(title_str, fontsize=24, pad=22, weight='bold')
    fig.canvas.draw()

def on_key(event):
    global t, current_proj_idx, bSun
    if event.key == 'left':
        t = t - 5 / (24 * 60)  # -5分钟
    elif event.key == 'right':
        t = t + 5 / (24 * 60)  # +5分钟
    elif event.key == 'm':
        current_proj_idx = (current_proj_idx + 1) % len(projections)  # 切换投影
    '''
    elif event.key == 'd':
        if bSun:
            bSun = 0
        else:
            bSun = 1
    '''
    update_plot(t)

# 创建地图
fig = plt.figure(figsize=(16, 9))
# ax = plt.axes(projection=ccrs.PlateCarree())  # 使用PlateCarree投影
# 初始化图像
update_plot(t)
# 绑定事件处理函数
fig.canvas.mpl_connect('key_press_event', on_key)
# 展示图像
plt.show()