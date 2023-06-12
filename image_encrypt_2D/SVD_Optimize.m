function Phi = SVD_Optimize(A)
    %-- 使用奇异值分解SVD来优化矩阵
    [M,N] = size(A);
    t = 5;
    
    %-- 测量矩阵A进行SVD分解, 计算对角矩阵S的平均值mean，得到对角矩阵S大于等于mean的总数f
    [~,S,~] = svd(A);
    mean = 0;
    for i=1:M
        mean = mean + S(i,i);
    end
    
    f = 0;
    for i=1:M
        if S(i,i)>=mean
            f = f+1;
        end
    end
    
    %-- 构造一个M*N的矩阵J，矩阵J的前f列乘以t
    J = ones(M,N);
    for j=1:f
        J(:,j) = J(:,j)*t;
    end
    
    %-- A与J点乘，
    A1 = A.*J;
    [U,S1,V] = svd(A1);
    for i=1:M
        S1(i,i) = 1;
    end
    Phi = U*S1*V';
              
end

