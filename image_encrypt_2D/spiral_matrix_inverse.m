function S = spiral_matrix_inverse(PS,M,N)
%PS -- 1*MN
%S -- M*N

    S = zeros(M,N);
    startx = 1; starty = 1;
    loop = M/2;
    count = 1;
    offset = 1;     %用来控制右边界 和 下边界
    
    while loop>0
        i = startx;
        j = starty;
        
        while j<=N-offset     %-- 最上行从左到右（左闭右开）
            S(i,j) = PS(count);
            j = j+1;  count = count+1;
        end
        
        while i<=M-offset     %-- 最右列从上到下（上闭下开）
            S(i,j) = PS(count);
            i = i+1; count = count+1;
        end
        
        while j>starty        %-- 最下行从右到左（右闭左开）
            S(i,j) = PS(count);
            j = j-1; count = count+1;
        end
        
        while i>startx        %-- 最左列从下到上（下闭上开）
            S(i,j) = PS(count);
            i = i-1; count = count+1;
        end
        
        %-- 进行下一轮
        startx = startx+1;
        starty = starty+1;
        offset = offset+1;
        loop = loop-1;
    end
   
end

