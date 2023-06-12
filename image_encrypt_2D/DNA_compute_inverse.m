%% DNA编码矩阵的逆运算
function result = DNA_compute_inverse(N_DNA,num)
     [m,n] = size(N_DNA);
     result = strings(m,n);
     for i=1:m
        for j=1:n
            switch N_DNA(i,j)
                case 'A'
                    switch num
                        case 1
                            result(i,j) = 'A';
                        case 2
                            result(i,j) = 'T';
                        case 3
                            result(i,j) = 'G';
                        case 4
                            result(i,j) = 'C';
                    end
                case 'C'
                    switch num
                        case 1
                            result(i,j) = 'C';
                        case 2
                            result(i,j) = 'A';
                        case 3
                            result(i,j) = 'T';
                        case 4
                            result(i,j) = 'G';
                    end
                case 'T'
                    switch num
                        case 1
                            result(i,j) = 'T';
                        case 2
                            result(i,j) = 'G';
                        case 3
                            result(i,j) = 'C';
                        case 4
                            result(i,j) = 'A';
                    end
                case 'G'
                    switch num
                        case 1
                            result(i,j) = 'G';
                        case 2
                            result(i,j) = 'C';
                        case 3
                            result(i,j) = 'A';
                        case 4
                            result(i,j) = 'T';
                    end
            end
        end
    end

end