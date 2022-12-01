function [x_result,final_iter,final_eps,runningTime,R] = SORIterate(B_SRM,F_SRM,maxIter,eps)
    % 该函数实现超松弛迭代法求解线性方程组X(k+1)=B_SRM*X(k)+F_SRM
    % 输入参数分别为SOR迭代公式中的矩阵B和F(松弛因子已经选定),最大迭代次数(产生离散误差)以及误差限
    % 输出参数为该线性方程组的解,最终迭代次数,最后一次迭代的误差,运算时间(ms)和渐近收敛速度
    % 用相邻两次的求解结果的二范数误差小于误差限或者迭代次数达到上限作为终止标准
    % 因为selectRelaxationFactor脚本的对w的搜寻范围为0到2,故在该脚本中无需检验迭代收敛的必要条件
    tic;
    x_num = size(B_SRM,2);                      % 未知数个数
    x_result = zeros(x_num,1);                  % 迭代初始值
    Rou_B_SRM = vrho(B_SRM);                    % 求出B的谱半径,也可以调用max(abs(eig(BG_S)))
    if Rou_B_SRM >= 1                           % 收敛的充要条件
        error('矩阵B的谱半径大于等于1,迭代不收敛')
    end
    final_iter = 0;
    while final_iter < maxIter
        x_result1 = B_SRM*x_result+F_SRM;       % 迭代求解(产生离散误差)
        final_eps = norm(x_result1-x_result);
        final_iter = final_iter + 1;
        if final_eps <= eps
            x_result = x_result1;
            break
        end
        x_result = x_result1;
    end
    runningTime = toc*1e3;                      % 程序运行时间(单位为ms)
    R = -log(Rou_B_SRM);                        % 渐进收敛速度
end