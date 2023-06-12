function PS = spiral_matrix(S)
    %-- S -M*3N
    %-- PS -1*3*M*N
    [M,N] = size(S);
    startx = 1; starty = 1;
    loop = M/2;
    count = 1;
    offset = 1;     %用来控制右边界 和 下边界

    PS = zeros(1,M*N);
    while loop>0
        i = startx;
        j = starty;
        
        while j<=N-offset     %-- 最上行从左到右（左闭右开）
            PS(count) = S(i,j);
            j = j+1;  count = count+1;
        end
        
        while i<=M-offset     %-- 最右列从上到下（上闭下开）
            PS(count) = S(i,j);
            i = i+1; count = count+1;
        end
        
        while j>starty        %-- 最下行从右到左（右闭左开）
            PS(count) = S(i,j);
            j = j-1; count = count+1;
        end
        
        while i>startx        %-- 最左列从下到上（下闭上开）
            PS(count) = S(i,j);
            i = i-1; count = count+1;
        end
        
        %-- 进行下一轮
        startx = startx+1;
        starty = starty+1;
        offset = offset+1;
        loop = loop-1;
    end
    
    %-- 如果M为奇数  那么还剩下一行
    if mod(M,2)==1
        i = startx; j = starty;
        while j<=N-offset+1
            PS(count) = S(i,j);
            j = j+1; count = count+1;
        end

end

