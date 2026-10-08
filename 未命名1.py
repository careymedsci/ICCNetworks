import numpy as np
import pandas as pd
import cv2
import matplotlib.pyplot as plt
from tslearn.clustering import KShape
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
from scipy.ndimage import gaussian_filter1d
import os

# === 参数 ===
seconds_per_frame = 2.0  # 请根据实际设置帧率（单位：秒/帧）
num_clusters = 10         # 聚类簇数量

# === 1. 读取波形数据（跳过第一列） ===
raw_df = pd.read_csv('dfof_mean_values.txt', sep='\t', header=None)
all_traces = raw_df.iloc[1:, 1:].astype(float).to_numpy()  # shape: (frames, cells)

frames, rois = all_traces.shape
frames_per_min = 60 / seconds_per_frame

# === 2. ΔF校正：逐帧减去当前帧所有细胞的平均值 ===
frame_means = np.mean(all_traces, axis=1, keepdims=True)  # shape: (frames, 1)
all_traces -= frame_means

# === 3. 高斯滤波平滑每个ROI ===
sigma = 7 / seconds_per_frame
for i in range(rois):
    all_traces[:, i] = gaussian_filter1d(all_traces[:, i], sigma)

# === 4. KShape聚类前预处理 ===
data = all_traces.T  # shape: (num_cells, time)
scaler = TimeSeriesScalerMeanVariance()
data_scaled = scaler.fit_transform(data[:, :, np.newaxis])

# === 5. 执行KShape聚类 ===
kshape = KShape(n_clusters=num_clusters, verbose=True, random_state=0)
kshape.fit(data_scaled)
labels = kshape.labels_

# === 6. 读取坐标文件 ===
coord_df = pd.read_csv('dfof_coordintes.txt', sep='\t', header=None)
coords = coord_df.iloc[1, 1:].astype(float).to_numpy().reshape(-1, 2)  # shape: (cells, 2)

# === 7. 读取/生成 max-z 图像 ===
if os.path.exists("max_z_project.tif"):
    max_z = cv2.imread("max_z_project.tif", cv2.IMREAD_GRAYSCALE)
else:
    max_z = np.random.randint(0, 256, (1024, 1280), dtype=np.uint8)
    cv2.imwrite("max_z_project.tif", max_z)

max_z_rgb = cv2.cvtColor(max_z, cv2.COLOR_GRAY2BGR)

# === 8. 给每个簇分配颜色并绘图 ===
colors = plt.cm.get_cmap("tab10", num_clusters)

for i, (x, y) in enumerate(coords):
    cluster = labels[i]
    rgb = np.array(colors(cluster)[:3]) * 255
    color = tuple(int(c) for c in rgb)
    cv2.circle(max_z_rgb, (int(x), int(y)), 5, color, -1)

# === 9. 显示聚类图 ===
plt.figure(figsize=(12, 8))
plt.imshow(max_z_rgb[..., ::-1])  # BGR转RGB
plt.title("Clustered Cells Overlayed on Max-Z Image")
plt.axis("off")
plt.show()

# === 10. 每个簇的平均波形 ===
plt.figure(figsize=(10, 6))
for cl in range(num_clusters):
    members = data[labels == cl]
    mean_trace = members.mean(axis=0)
    plt.plot(mean_trace.flatten(), label=f"Cluster {cl}")
plt.title("Mean Trace for Each Cluster")
plt.xlabel("Frame")
plt.ylabel("dF (a.u.)")
plt.legend()
plt.tight_layout()
plt.show()

# === 11. 保存聚类信息 ===
cluster_info = pd.DataFrame({
    'Cell_ID': np.arange(1, rois + 1),
    'Cluster_Label': labels
})
cluster_info.to_csv("cell_clusters.csv", index=False)
print("聚类结果已保存为 cell_clusters.csv")


# 以下进行类聚类可视化
from sklearn.cluster import DBSCAN

# === 设置颜色映射（每类不同色） ===
cluster_colors = plt.cm.get_cmap("tab10", num_clusters)

# 拷贝 max_z 图像用于绘制所有聚类簇
overlay_all = max_z_rgb.copy()

for target_cluster in range(num_clusters):
    # 获取该簇的细胞坐标
    target_coords = coords[labels == target_cluster]
    if len(target_coords) == 0:
        continue

    # 使用 DBSCAN 做空间聚类（可调参数）
    db = DBSCAN(eps=30, min_samples=3).fit(target_coords)
    spatial_labels = db.labels_

    unique_groups = np.unique(spatial_labels)

    # 当前聚类簇的颜色（用于所有子团轮廓）
    color_rgb = np.array(cluster_colors(target_cluster)[:3]) * 255
    color_bgr = tuple(int(c) for c in color_rgb[::-1])  # RGB → BGR for OpenCV

    for grp in unique_groups:
        if grp == -1:
            continue  # 忽略孤立点

        points = target_coords[spatial_labels == grp].astype(np.int32)

        if len(points) >= 3:
            hull = cv2.convexHull(points)
            cv2.polylines(overlay_all, [hull], isClosed=True, color=color_bgr, thickness=2)
        else:
            for (x, y) in points:
                cv2.circle(overlay_all, (int(x), int(y)), 3, color_bgr, -1)

# === 显示并保存结果 ===
plt.figure(figsize=(14, 10))
plt.imshow(overlay_all[..., ::-1])
plt.title("Spatial Subclusters for All Clusters")
plt.axis("off")
plt.show()

cv2.imwrite("all_clusters_spatial_subgroups.png", overlay_all)
print("已保存为 all_clusters_spatial_subgroups.png")


# 计算小世界特性。
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# 计算网络的平均聚类系数和平均最短路径长度
def calculate_network_properties(adj_matrix):
    G = nx.from_numpy_array(adj_matrix)
    clustering_coeff = np.mean(list(nx.clustering(G).values()))
    try:
        avg_path_length = nx.average_shortest_path_length(G)
    except nx.NetworkXError:  # 如果图是非连通的
        avg_path_length = np.nan
    return clustering_coeff, avg_path_length

# 创建小世界网络并计算属性
def create_small_world_network(num_nodes, p=0.1):
    # 创建一个环形网络
    G = nx.watts_strogatz_graph(num_nodes, k=4, p=p)
    adj_matrix = nx.to_numpy_array(G)
    return adj_matrix

# 计算并绘制对比图
def plot_comparison(real_adj_matrix, num_nodes):
    # 1. 计算原始网络属性
    real_clustering, real_avg_path_length = calculate_network_properties(real_adj_matrix)
    
    # 2. 使用Watts-Strogatz生成小世界网络
    small_world_adj_matrix = create_small_world_network(num_nodes, p=0.1)  # p值控制小世界特性
    small_world_clustering, small_world_avg_path_length = calculate_network_properties(small_world_adj_matrix)
    
    # 3. 绘制对比图
    labels = ['Clustering Coefficient', 'Average Path Length']
    real_values = [real_clustering, real_avg_path_length]
    small_world_values = [small_world_clustering, small_world_avg_path_length]
    
    x = np.arange(len(labels))  # x轴位置
    width = 0.35  # 每个条形图的宽度
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(x - width/2, real_values, width, label='Real Network')
    ax.bar(x + width/2, small_world_values, width, label='Small-World Network')
    
    ax.set_ylabel('Values')
    ax.set_title('Comparison of Real Network vs Small-World Network')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    
    plt.show()

# 假设原始的邻接矩阵（real_adj_matrix）是一个代表真实细胞网络的矩阵
# 这里使用随机矩阵作为示例
num_cells = 100  # 假设有100个细胞
real_adj_matrix = np.random.random((num_cells, num_cells))
real_adj_matrix = (real_adj_matrix + real_adj_matrix.T) / 2  # 保证对称性

# 绘制比较图
plot_comparison(real_adj_matrix, num_cells)




Created on Fri Apr 18 22:01:18 2025

@author: xiaoliu
"""

