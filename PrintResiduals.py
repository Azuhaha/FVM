
def printResiduals(numOfIterations,reg):
    #Print
    print('|======================================================================================|')
    print(f'                                    Global Iter {numOfIterations}')
    print('|--------------------------------------------------------------------------------------|')
    print('|--------------------------------------------------------------------------------------|')
    print('|     Equation     |      RMS      |      MAX      | initialResidual |  finalResidual  |')
    print('|--------------------------------------------------------------------------------------|')

    rmsResidual_U=reg.model.U.residuals.rmsResidual
    maxResidual_U = reg.model.U.residuals.maxResidual
    initResidual_U = reg.model.U.residuals.initResidual
    finalResidual_U = reg.model.U.residuals.finalResidual

    print(f'|       U-x        |  {rmsResidual_U[0]:.3e}    |   {maxResidual_U[0]:.3e}   |    {initResidual_U[0]:.3e}    |    {finalResidual_U[0]:.3e}    |')
    print(f'|       U-y        |  {rmsResidual_U[1]:.3e}    |   {maxResidual_U[1]:.3e}   |    {initResidual_U[1]:.3e}    |    {finalResidual_U[1]:.3e}    |')
    print(f'|       U-z        |  {rmsResidual_U[2]:.3e}    |   {maxResidual_U[2]:.3e}   |    {initResidual_U[2]:.3e}    |    {finalResidual_U[2]:.3e}    |')

    rmsResidual_p=float(reg.model.p.residuals.rmsResidual)
    maxResidual_p = float(reg.model.p.residuals.maxResidual)
    initResidual_p = float(reg.model.p.residuals.initResidual)
    finalResidual_p = float(reg.model.p.residuals.finalResidual)
    print(f'|        p         |  {rmsResidual_p:.3e}    |   {maxResidual_p:.3e}   |    {initResidual_p:.3e}    |    {finalResidual_p:.3e}    |')
    print('|======================================================================================|')
    '''
    #Plot
    x = np.linspace(0, numOfIterations,numOfIterations+1)# 从0到10生成100个均匀分布的点

    y_Ux = np.append(RMSResidual_lst[0],rmsResidual_U[0])
    y_Uy = np.append(RMSResidual_lst[1],rmsResidual_U[1])
    y_Uz = np.append(RMSResidual_lst[2],rmsResidual_U[2])
    y_p  = np.append(RMSResidual_lst[3],rmsResidual_p)

    # 创建画布和子图
    plt.figure(figsize=(10, 6))  # 设置画布大小

    # 绘制曲线

    plt.plot(x, y_Ux, label='Ux', color='red', linestyle='-', linewidth=2, marker='*', markersize=3)
    plt.plot(x, y_Uy, label='Uy', color='green', linestyle='-', linewidth=2, marker='o', markersize=3)
    plt.plot(x, y_Uz, label='Uz', color='blue', linestyle='-', linewidth=2, marker='^', markersize=3)
    plt.plot(x, y_Uz, label='p', color='black', linestyle='-', linewidth=2, marker='^', markersize=3)

    # 设置x轴刻度为10的间隔
    plt.xticks(np.arange(0, numOfIterations + 1, 10))  # 从0开始，到numOfIterations，间隔10

    plt.yscale('log')
    # 添加标题和坐标轴标签
    plt.title('Residuals Plot', fontsize=15)
    plt.xlabel('Global Iterations', fontsize=12)
    plt.ylabel('Scaled RMS Residuals', fontsize=12)

    # 添加网格线
    plt.grid(True, linestyle=':', alpha=0.7)

    # 添加图例
    plt.legend(fontsize=10)

    # 设置坐标轴范围
    #plt.xlim()
    #plt.ylim()

    # 显示图形
    plt.tight_layout()  # 调整布局，避免标签被截断
    plt.show()

    return [y_Ux,y_Uy,y_Uz,y_p]
    '''