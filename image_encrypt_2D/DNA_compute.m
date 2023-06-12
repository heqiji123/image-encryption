%% 对已经DNA编码的矩阵进行扩散操作
function result = DNA_compute(DNA_R,num)
    [m,n] = size(DNA_R);
    result = strings(m,n);
    for i=1:m
        for j=1:n
            switch DNA_R(i,j)
                case 'A'
                    switch num
                        case 1
                            result(i,j) = 'A';
                        case 2
                            result(i,j) = 'C';
                        case 3
                            result(i,j) = 'G';
                        case 4
                            result(i,j) = 'T';
                    end
                case 'C'
                    switch num
                        case 1
                            result(i,j) = 'C';
                        case 2
                            result(i,j) = 'G';
                        case 3
                            result(i,j) = 'T';
                        case 4
                            result(i,j) = 'A';
                    end
                case 'T'
                    switch num
                        case 1
                            result(i,j) = 'T';
                        case 2
                            result(i,j) = 'A';
                        case 3
                            result(i,j) = 'C';
                        case 4
                            result(i,j) = 'G';
                    end
                case 'G'
                    switch num
                        case 1
                            result(i,j) = 'G';
                        case 2
                            result(i,j) = 'T';
                        case 3
                            result(i,j) = 'A';
                        case 4
                            result(i,j) = 'C';
                    end
            end
        end
    end
end