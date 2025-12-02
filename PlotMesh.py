import numpy as np
import matplotlib.pyplot as plt

def plotMesh(reg):
    #点图
    plotMeshNode(reg)

def plotMeshNode(reg):
    # 生成三维坐标数据
    x = reg.mesh.nodeCentroids[:, 0]  # x坐标
    y = reg.mesh.nodeCentroids[:, 1]  # y坐标
    z = reg.mesh.nodeCentroids[:, 2]  # z坐标

    # 创建3D图形
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # 绘制三维散点图
    scatter = ax.scatter(x, y, z, s=10, alpha=0.8, edgecolors='w', linewidth=0.5)

    # 设置坐标轴标签
    ax.set_xlabel('X-coordinate', fontsize=12)
    ax.set_ylabel('Y-coordinate', fontsize=12)
    ax.set_zlabel('Z-coordinate', fontsize=12)

    # 设置标题
    ax.set_title('Mesh', fontsize=15, pad=20)

    ax.set_box_aspect([1, 1, 1])

    # 调整视角
    ax.view_init(elev=30, azim=45)  # elevation(仰角), azimuth(方位角)

    # 显示图形
    plt.tight_layout()
    plt.show()

def plotMeshLine(reg):
    # 生成三维坐标数据
    x = reg.mesh.nodeCentroids[:, 0]  # x坐标
    y = reg.mesh.nodeCentroids[:, 1]  # y坐标
    z = reg.mesh.nodeCentroids[:, 2]  # z坐标

    # 创建3D图形
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # 绘制三维散点图
    scatter = ax.scatter(x, y, z, s=10, alpha=0.8, edgecolors='w', linewidth=0.5)

    # 设置坐标轴标签
    ax.set_xlabel('X-coordinate', fontsize=12)
    ax.set_ylabel('Y-coordinate', fontsize=12)
    ax.set_zlabel('Z-coordinate', fontsize=12)

    # 设置标题
    ax.set_title('Mesh', fontsize=15, pad=20)

    ax.set_box_aspect([1, 1, 1])

    # 调整视角
    ax.view_init(elev=30, azim=45)  # elevation(仰角), azimuth(方位角)

    # 显示图形
    plt.tight_layout()
    plt.show()