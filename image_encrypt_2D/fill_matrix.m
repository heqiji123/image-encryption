function Mat = fill_matrix(Matrix,block_fill)
% 填充矩阵
[M,N] = size(Matrix);
if block_fill == 0
    Mat = matrix;
else 
    Mat = zeros(M+block_fill,N+block_fill);
    for i = 1:M
        for j = 1:N
            Mat(i,j) = Matrix(i,j);
        end
    end
end
end

