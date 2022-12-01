function [x_result,final_iter,final_eps,runningTime,R] = GaussSeidelIterate(AugmentationMatrix,maxIter,eps)
    % 该函数实现用Gauss-Seidel迭代法求解线性方程组X(k+1)=BG_S*X(k)+FG_S
    % 输入参数分别为线性方程组的增广矩阵和最大迭代次数(产生离散误差)以及误差限
    % 输出参数为该线性方程组的解,最终迭代次数,最后一次迭代的误差,运算时间(ms)和渐近收敛速度
    % 用相邻两次的求解结果的二范数误差小于误差限或者迭代次数达到上限作为终止标准
    tic;
    if size(AugmentationMatrix,1) + 1 ~= size(AugmentationMatrix,2)
        error('输入矩阵维度应为n*(n+1)')
    end
    x_num = size(AugmentationMatrix,1);  % 未知数个数
    x_result = zeros(x_num,1);           % 迭代初始值
    A = AugmentationMatrix(:,1:x_num);   % 系数矩阵
    b = AugmentationMatrix(:,x_num+1);
    A_diag = diag(A);                    % 主对角元素对应的向量
    if double(ismember(0,A_diag)) == 1   % 迭代适用的必要条件
        error('系数矩阵的主对角元素不能出现0');
    end
    D = diag(A_diag);                    % 取出主对角元素矩阵
    L = -tril(A)+D;                      % 取出负下三角矩阵(不包含主对角)
    U = -triu(A)+D;                      % 取出负上三角矩阵(不包含主对角)
    BG_S = (D-L)\U;                      % 求出Gauss-Seidel迭代公式中的矩阵B
    FG_S = (D-L)\b;                      % 求出Gauss-Seidel迭代公式中的矩阵F
    Rou_BG_S = vrho(BG_S);               % 求出B的谱半径,也可以调用max(abs(eig(BG_S)))
    if Rou_BG_S >= 1                     % 收敛的充要条件
        error('矩阵B的谱半径大于等于1,迭代不收敛')
    end
    final_iter = 0;
    while final_iter < maxIter
        x_result1 = BG_S*x_result+FG_S;   % 迭代求解(产生离散误差)
        final_eps = norm(x_result1-x_result);
        final_iter = final_iter + 1;
        if final_eps <= eps
            x_result = x_result1;
            break
        end
        x_result = x_result1;
    end
    runningTime = toc*1e3;               % 程序运行时间(单位为ms)
    R = -log(Rou_BG_S);                  % 渐进收敛速度
end