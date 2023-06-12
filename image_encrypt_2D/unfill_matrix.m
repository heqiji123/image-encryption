function Mat = unfill_matrix(Matrix,block_fill)
%取出填充的元素
[~,N] = size(Matrix);
if block_fill == 0
    Mat = Matrix;
else
    Mat = zeros(N - block_fill, N - block_fill);
    for i = 1:N - block_fill
        for j = 1:N - block_fill
            Mat(i,j) = Matrix(i,j);
        end
    end
end
end

