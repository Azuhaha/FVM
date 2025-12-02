import matplotlib.pyplot as plt
import numpy as np

def plotResiduals(StepOfIterations, reg):
    global residuals_data
    if StepOfIterations==1:
        residuals_data=[[] for _ in range(4)]
        plt.figure(figsize=(10, 6))
    Ux_value=reg.model.U.residuals.rmsResidual[0]
    Uy_value = reg.model.U.residuals.rmsResidual[1]
    Uz_value = reg.model.U.residuals.rmsResidual[2]
    p_value = reg.model.p.residuals.rmsResidual
    residuals_data[0].append(Ux_value)
    residuals_data[1].append(Uy_value)
    residuals_data[2].append(Uz_value)
    residuals_data[3].append(p_value)

    # 准备绘图数据
    steps = list(range(1,StepOfIterations+1))

    # 清除当前图表并重新绘制
    plt.clf()
    # 绘制第一条曲线（Ux）
    plt.plot(steps, residuals_data[0], 'o-', color='red',markersize=4, linewidth=1, label='Ux')
    # 绘制第二条曲线（Uy）
    plt.plot(steps, residuals_data[1], 's-', color='green',markersize=4, linewidth=1, label='Uy')

    # 绘制第三条曲线（Uz）
    plt.plot(steps, residuals_data[2], '^-', color='blue', markersize=4, linewidth=1, label='Uz')

    # 绘制第四条曲线（p）
    plt.plot(steps, residuals_data[3], 'd-', color='black',markersize=4, linewidth=1, label='p')

    plt.xlabel('Global Iterations')
    plt.ylabel('Scaled RMS Residuals')
    plt.title('Residuals Plot')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim(0, StepOfIterations+0.5)  # 设置x轴范围，使显示更美观
    plt.xticks(range((1+StepOfIterations//10)*10))  # 确保x轴显示0到10的所有整数
    plt.yscale('log')
    plt.legend(fontsize=10, loc='best')



    # 实时更新图表
    plt.draw()
    plt.pause(0.1)  # 暂停一小段时间，让图表有时间更新

