function [optimal_w,B_SRM,F_SRM] = selectRelaxationFactor(AugmentationMatrix,step)
    % 用于寻找较优的松弛因子,以使得收敛速度较快
    % 输入参数为超松弛迭代法的迭代矩阵B,以及搜寻的步长(范围为0到2)
    % 输出较优的松弛因子
    if size(AugmentationMatrix,1) + 1 ~= size(AugmentationMatrix,2)
        error('输入矩阵维度应为n*(n+1)')
    end
    A = AugmentationMatrix(:,1:end-1);      % 系数矩阵
    b = AugmentationMatrix(:,end);
    A_diag = diag(A);                       % 主对角元素对应的向量
    D = diag(A_diag);                       % 取出主对角元素矩阵
    L = -tril(A)+D;                         % 取出负下三角矩阵(不包含主对角)
    U = -triu(A)+D;                         % 取出负上三角矩阵(不包含主对角)
    idx = 1;
    Rou_B = zeros(1,2/step+1);
    for w = 0:step:2
        B_SRM = (D-w*L)\((1-w)*D+w*U);      % 求出SOR迭代公式中的矩阵B
        Rou_B(idx) = vrho(B_SRM);           % 求出B的谱半径
        idx = idx + 1;
    end
    optimal_idx = find(Rou_B==min(Rou_B));
    optimal_w = (optimal_idx-1)*step;
    B_SRM = (D-optimal_w*L)\((1-optimal_w)*D+optimal_w*U);
    F_SRM = optimal_w*((D-optimal_w*L)\b);  % 别漏了这里的括号
end