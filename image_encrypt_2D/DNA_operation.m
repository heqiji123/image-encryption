function result = DNA_operation(A,B,num)
%DNA_OPERATION DNA操作函数 + - xor
% A B 是DNA数组 
%-- num是操作的类型
[m,n] = size(A);
result = strings(m,n);

for i=1:m
    for j=1:n
        switch num
            case 1
                switch A(i,j)
                    case "A"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "C";
                            case "T"
                                result(i,j) = "G";
                            case "C"
                                result(i,j) = "A";
                            case "G"
                                result(i,j) = "T";
                        end
                    case "T"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "G";
                            case "T"
                                result(i,j) = "C";
                            case "C"
                                result(i,j) = "T";
                            case "G"
                                result(i,j) = "A";
                        end
                    case "C"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "A";
                            case "T"
                                result(i,j) = "T";
                            case "C"
                                result(i,j) = "C";
                            case "G"
                                result(i,j) = "G";
                        end
                    case "G"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "T";
                            case "T"
                                result(i,j) = "A";
                            case "C"
                                result(i,j) = "G";
                            case "G"
                                result(i,j) = "C";
                        end
                end
            case 2
                switch A(i,j)
                    case "A"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "C";
                            case "T"
                                result(i,j) = "A";
                            case "C"
                                result(i,j) = "G";
                            case "G"
                                result(i,j) = "T";
                        end
                    case "T"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "G";
                            case "T"
                                result(i,j) = "C";
                            case "C"
                                result(i,j) = "T";
                            case "G"
                                result(i,j) = "A";
                        end
                    case "C"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "A";
                            case "T"
                                result(i,j) = "T";
                            case "C"
                                result(i,j) = "C";
                            case "G"
                                result(i,j) = "G";
                        end
                    case "G"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "T";
                            case "T"
                                result(i,j) = "G";
                            case "C"
                                result(i,j) = "A";
                            case "G"
                                result(i,j) = "C";
                        end
                end
            case 3
                switch A(i,j)
                    case "A"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "A";
                            case "T"
                                result(i,j) = "T";
                            case "C"
                                result(i,j) = "C";
                            case "G"
                                result(i,j) = "G";
                        end
                    case "T"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "T";
                            case "T"
                                result(i,j) = "A";
                            case "C"
                                result(i,j) = "G";
                            case "G"
                                result(i,j) = "C";
                        end
                    case "C"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "C";
                            case "T"
                                result(i,j) = "G";
                            case "C"
                                result(i,j) = "A";
                            case "G"
                                result(i,j) = "T";
                        end
                    case "G"
                        switch  B(i,j)
                            case "A"
                                result(i,j) = "G";
                            case "T"
                                result(i,j) = "C";
                            case "C"
                                result(i,j) = "T";
                            case "G"
                                result(i,j) = "A";
                        end
                end
        end
    end
end
end

