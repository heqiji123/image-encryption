function result = DNA_decode_(A,num)
[m,n] = size(A);
result = zeros(m,n/4);

for i=1:m
    for j=1:n/4
        switch num
            case 1
                switch A(i,((j-1)*4+1))
                    case "A"
                        a1 = 0;
                    case "T"
                        a1 = 3;
                    case "C"
                        a1 = 1;
                    case "G"
                        a1 = 2;
                end
                switch A(i,((j-1)*4+2))
                    case "A"
                        a2 = 0;
                    case "T"
                        a2 = 3;
                    case "C"
                        a2 = 1;
                    case "G"
                        a2 = 2;
                end
                switch A(i,((j-1)*4+3))
                    case "A"
                        a3 = 0;
                    case "T"
                        a3 = 3;
                    case "C"
                        a3 = 1;
                    case "G"
                        a3 = 2;
                end
                switch A(i,((j-1)*4+4))
                    case "A"
                        a4 = 0;
                    case "T"
                        a4 = 3;
                    case "C"
                        a4 = 1;
                    case "G"
                        a4 = 2;
                end
            case 2
                switch A(i,((j-1)*4+1))
                    case "A"
                        a1 = 0;
                    case "T"
                        a1 = 3;
                    case "C"
                        a1 = 2;
                    case "G"
                        a1 = 1;
                end
                switch A(i,((j-1)*4+2))
                    case "A"
                        a2 = 0;
                    case "T"
                        a2 = 3;
                    case "C"
                        a2 = 2;
                    case "G"
                        a2 = 1;
                end
                switch A(i,((j-1)*4+3))
                    case "A"
                        a3 = 0;
                    case "T"
                        a3 = 3;
                    case "C"
                        a3 = 2;
                    case "G"
                        a3 = 1;
                end
                switch A(i,((j-1)*4+4))
                    case "A"
                        a4 = 0;
                    case "T"
                        a4 = 3;
                    case "C"
                        a4 = 2;
                    case "G"
                        a4 = 1;
                end
            case 3
                switch A(i,((j-1)*4+1))
                    case "A"
                        a1 = 3;
                    case "T"
                        a1 = 0;
                    case "C"
                        a1 = 1;
                    case "G"
                        a1 = 2;
                end
                switch A(i,((j-1)*4+2))
                    case "A"
                        a2 = 3;
                    case "T"
                        a2 = 0;
                    case "C"
                        a2 = 1;
                    case "G"
                        a2 = 2;
                end
                switch A(i,((j-1)*4+3))
                    case "A"
                        a3 = 3;
                    case "T"
                        a3 = 0;
                    case "C"
                        a3 = 1;
                    case "G"
                        a3 = 2;
                end
                switch A(i,((j-1)*4+4))
                    case "A"
                        a4 = 3;
                    case "T"
                        a4 = 0;
                    case "C"
                        a4 = 1;
                    case "G"
                        a4 = 2;
                end
            case 4
                switch A(i,((j-1)*4+1))
                    case "A"
                        a1 = 2;
                    case "T"
                        a1 = 1;
                    case "C"
                        a1 = 0;
                    case "G"
                        a1 = 3;
                end
                switch A(i,((j-1)*4+2))
                    case "A"
                        a2 = 2;
                    case "T"
                        a2 = 1;
                    case "C"
                        a2 = 0;
                    case "G"
                        a2 = 3;
                end
                switch A(i,((j-1)*4+3))
                    case "A"
                        a3 = 2;
                    case "T"
                        a3 = 1;
                    case "C"
                        a3 = 0;
                    case "G"
                        a3 = 3;
                end
                switch A(i,((j-1)*4+4))
                    case "A"
                        a4 = 2;
                    case "T"
                        a4 = 1;
                    case "C"
                        a4 = 0;
                    case "G"
                        a4 = 3;
                end
            case 5
                switch A(i,((j-1)*4+1))
                    case "A"
                        a1 = 1;
                    case "T"
                        a1 = 2;
                    case "C"
                        a1 = 3;
                    case "G"
                        a1 = 0;
                end
                switch A(i,((j-1)*4+2))
                    case "A"
                        a2 = 1;
                    case "T"
                        a2 = 2;
                    case "C"
                        a2 = 3;
                    case "G"
                        a2 = 0;
                end
                switch A(i,((j-1)*4+3))
                    case "A"
                        a3 = 1;
                    case "T"
                        a3 = 2;
                    case "C"
                        a3 = 3;
                    case "G"
                        a3 = 0;
                end
                switch A(i,((j-1)*4+4))
                    case "A"
                        a4 = 1;
                    case "T"
                        a4 = 2;
                    case "C"
                        a4 = 3;
                    case "G"
                        a4 = 0;
                end
            case 6
                switch A(i,((j-1)*4+1))
                    case "A"
                        a1 = 2;
                    case "T"
                        a1 = 1;
                    case "C"
                        a1 = 3;
                    case "G"
                        a1 = 0;
                end
                switch A(i,((j-1)*4+2))
                    case "A"
                        a2 = 2;
                    case "T"
                        a2 = 1;
                    case "C"
                        a2 = 3;
                    case "G"
                        a2 = 0;
                end
                switch A(i,((j-1)*4+3))
                    case "A"
                        a3 = 2;
                    case "T"
                        a3 = 1;
                    case "C"
                        a3 = 3;
                    case "G"
                        a3 = 0;
                end
                switch A(i,((j-1)*4+4))
                    case "A"
                        a4 = 2;
                    case "T"
                        a4 = 1;
                    case "C"
                        a4 = 3;
                    case "G"
                        a4 = 0;
                end
            case 7
                switch A(i,((j-1)*4+1))
                    case "A"
                        a1 = 1;
                    case "T"
                        a1 = 2;
                    case "C"
                        a1 = 0;
                    case "G"
                        a1 = 3;
                end
                switch A(i,((j-1)*4+2))
                    case "A"
                        a2 = 1;
                    case "T"
                        a2 = 2;
                    case "C"
                        a2 = 0;
                    case "G"
                        a2 = 3;
                end
                switch A(i,((j-1)*4+3))
                    case "A"
                        a3 = 1;
                    case "T"
                        a3 = 2;
                    case "C"
                        a3 = 0;
                    case "G"
                        a3 = 3;
                end
                switch A(i,((j-1)*4+4))
                    case "A"
                        a4 = 1;
                    case "T"
                        a4 = 2;
                    case "C"
                        a4 = 0;
                    case "G"
                        a4 = 3;
                end
            case 8
                switch A(i,((j-1)*4+1))
                    case "A"
                        a1 = 3;
                    case "T"
                        a1 = 0;
                    case "C"
                        a1 = 2;
                    case "G"
                        a1 = 1;
                end
                switch A(i,((j-1)*4+2))
                    case "A"
                        a2 = 3;
                    case "T"
                        a2 = 0;
                    case "C"
                        a2 = 2;
                    case "G"
                        a2 = 1;
                end
                switch A(i,((j-1)*4+3))
                    case "A"
                        a3 = 3;
                    case "T"
                        a3 = 0;
                    case "C"
                        a3 = 2;
                    case "G"
                        a3 = 1;
                end
                switch A(i,((j-1)*4+4))
                    case "A"
                        a4 = 3;
                    case "T"
                        a4 = 0;
                    case "C"
                        a4 = 2;
                    case "G"
                        a4 = 1;
                end
        end
       result(i,j) = a1*64+a2*16+a3*4+a4; 
    end
end

end

