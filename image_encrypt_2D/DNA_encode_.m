%% 子函数 DNA编码
function sub_block_DNA=DNA_encode_(sub_block,num)
[M,N] = size(sub_block);
sub_block_DNA=strings(M,4*N);

for i=1:M
    for j=1:N
        a1=bitand(sub_block(i,j),192)/64;  %取第1、2两位
        a2=bitand(sub_block(i,j),48)/16;   %取第3，4两位
        a3=bitand(sub_block(i,j),12)/4;    %取第5、6两位
        a4=bitand(sub_block(i,j),3);       %取第7、8两位
        switch num
            case 1
                if a1==0                   %-- 编码第 1 2位
                    sub_block_DNA(i,((j-1)*4+1)) = 'A';
                elseif a1==1
                    sub_block_DNA(i,((j-1)*4+1)) = 'C';
                elseif a1==2
                    sub_block_DNA(i,((j-1)*4+1)) = 'G';
                else
                    sub_block_DNA(i,((j-1)*4+1)) = 'T';
                end
                if a2==0                   %-- 编码第 3 4位
                    sub_block_DNA(i,((j-1)*4+2)) = 'A';
                elseif a2==1
                    sub_block_DNA(i,((j-1)*4+2)) = 'C';
                elseif a2==2
                    sub_block_DNA(i,((j-1)*4+2)) = 'G';
                else
                    sub_block_DNA(i,((j-1)*4+2)) = 'T';
                end
                if a3==0                   %-- 编码第 5 6位
                    sub_block_DNA(i,((j-1)*4+3)) = 'A';
                elseif a3==1
                    sub_block_DNA(i,((j-1)*4+3)) = 'C';
                elseif a3==2
                    sub_block_DNA(i,((j-1)*4+3)) = 'G';
                else
                    sub_block_DNA(i,((j-1)*4+3)) = 'T';
                end
                if a4==0
                    sub_block_DNA(i,((j-1)*4+4)) = 'A';
                elseif a4==1
                    sub_block_DNA(i,((j-1)*4+4)) = 'C';
                elseif a4==2
                    sub_block_DNA(i,((j-1)*4+4)) = 'G';
                else
                    sub_block_DNA(i,((j-1)*4+4)) = 'T';
                end
            case 2
                if a1==0                   %-- 编码第 1 2位
                    sub_block_DNA(i,((j-1)*4+1)) = 'A';
                elseif a1==1
                    sub_block_DNA(i,((j-1)*4+1)) = 'G';
                elseif a1==2
                    sub_block_DNA(i,((j-1)*4+1)) = 'C';
                else
                    sub_block_DNA(i,((j-1)*4+1)) = 'T';
                end
                if a2==0                   %-- 编码第 3 4位
                    sub_block_DNA(i,((j-1)*4+2)) = 'A';
                elseif a2==1
                    sub_block_DNA(i,((j-1)*4+2)) = 'G';
                elseif a2==2
                    sub_block_DNA(i,((j-1)*4+2)) = 'C';
                else
                    sub_block_DNA(i,((j-1)*4+2)) = 'T';
                end
                if a3==0                   %-- 编码第 5 6位
                    sub_block_DNA(i,((j-1)*4+3)) = 'A';
                elseif a3==1
                    sub_block_DNA(i,((j-1)*4+3)) = 'G';
                elseif a3==2
                    sub_block_DNA(i,((j-1)*4+3)) = 'C';
                else
                    sub_block_DNA(i,((j-1)*4+3)) = 'T';
                end
                if a4==0
                    sub_block_DNA(i,((j-1)*4+4)) = 'A';
                elseif a4==1
                    sub_block_DNA(i,((j-1)*4+4)) = 'G';
                elseif a4==2
                    sub_block_DNA(i,((j-1)*4+4)) = 'C';
                else
                    sub_block_DNA(i,((j-1)*4+4)) = 'T';
                end
            case 3
                if a1==0                   %-- 编码第 1 2位
                    sub_block_DNA(i,((j-1)*4+1)) = 'T';
                elseif a1==1
                    sub_block_DNA(i,((j-1)*4+1)) = 'C';
                elseif a1==2
                    sub_block_DNA(i,((j-1)*4+1)) = 'G';
                else
                    sub_block_DNA(i,((j-1)*4+1)) = 'A';
                end
                if a2==0                   %-- 编码第 3 4位
                    sub_block_DNA(i,((j-1)*4+2)) = 'T';
                elseif a2==1
                    sub_block_DNA(i,((j-1)*4+2)) = 'C';
                elseif a2==2
                    sub_block_DNA(i,((j-1)*4+2)) = 'G';
                else
                    sub_block_DNA(i,((j-1)*4+2)) = 'A';
                end
                if a3==0                   %-- 编码第 5 6位
                    sub_block_DNA(i,((j-1)*4+3)) = 'T';
                elseif a3==1
                    sub_block_DNA(i,((j-1)*4+3)) = 'C';
                elseif a3==2
                    sub_block_DNA(i,((j-1)*4+3)) = 'G';
                else
                    sub_block_DNA(i,((j-1)*4+3)) = 'A';
                end
                if a4==0
                    sub_block_DNA(i,((j-1)*4+4)) = 'T';
                elseif a4==1
                    sub_block_DNA(i,((j-1)*4+4)) = 'C';
                elseif a4==2
                    sub_block_DNA(i,((j-1)*4+4)) = 'G';
                else
                    sub_block_DNA(i,((j-1)*4+4)) = 'A';
                end
            case 4
                if a1==0                   %-- 编码第 1 2位
                    sub_block_DNA(i,((j-1)*4+1)) = 'C';
                elseif a1==1
                    sub_block_DNA(i,((j-1)*4+1)) = 'T';
                elseif a1==2
                    sub_block_DNA(i,((j-1)*4+1)) = 'A';
                else
                    sub_block_DNA(i,((j-1)*4+1)) = 'G';
                end
                if a2==0                   %-- 编码第 3 4位
                    sub_block_DNA(i,((j-1)*4+2)) = 'C';
                elseif a2==1
                    sub_block_DNA(i,((j-1)*4+2)) = 'T';
                elseif a2==2
                    sub_block_DNA(i,((j-1)*4+2)) = 'A';
                else
                    sub_block_DNA(i,((j-1)*4+2)) = 'G';
                end
                if a3==0                   %-- 编码第 5 6位
                    sub_block_DNA(i,((j-1)*4+3)) = 'C';
                elseif a3==1
                    sub_block_DNA(i,((j-1)*4+3)) = 'T';
                elseif a3==2
                    sub_block_DNA(i,((j-1)*4+3)) = 'A';
                else
                    sub_block_DNA(i,((j-1)*4+3)) = 'G';
                end
                if a4==0
                    sub_block_DNA(i,((j-1)*4+4)) = 'C';
                elseif a4==1
                    sub_block_DNA(i,((j-1)*4+4)) = 'T';
                elseif a4==2
                    sub_block_DNA(i,((j-1)*4+4)) = 'A';
                else
                    sub_block_DNA(i,((j-1)*4+4)) = 'G';
                end
            case 5
                if a1==0                   %-- 编码第 1 2位
                    sub_block_DNA(i,((j-1)*4+1)) = 'G';
                elseif a1==1
                    sub_block_DNA(i,((j-1)*4+1)) = 'A';
                elseif a1==2
                    sub_block_DNA(i,((j-1)*4+1)) = 'T';
                else
                    sub_block_DNA(i,((j-1)*4+1)) = 'C';
                end
                if a2==0                   %-- 编码第 3 4位
                    sub_block_DNA(i,((j-1)*4+2)) = 'G';
                elseif a2==1
                    sub_block_DNA(i,((j-1)*4+2)) = 'A';
                elseif a2==2
                    sub_block_DNA(i,((j-1)*4+2)) = 'T';
                else
                    sub_block_DNA(i,((j-1)*4+2)) = 'C';
                end
                if a3==0                   %-- 编码第 5 6位
                    sub_block_DNA(i,((j-1)*4+3)) = 'G';
                elseif a3==1
                    sub_block_DNA(i,((j-1)*4+3)) = 'A';
                elseif a3==2
                    sub_block_DNA(i,((j-1)*4+3)) = 'T';
                else
                    sub_block_DNA(i,((j-1)*4+3)) = 'C';
                end
                if a4==0
                    sub_block_DNA(i,((j-1)*4+4)) = 'G';
                elseif a4==1
                    sub_block_DNA(i,((j-1)*4+4)) = 'A';
                elseif a4==2
                    sub_block_DNA(i,((j-1)*4+4)) = 'T';
                else
                    sub_block_DNA(i,((j-1)*4+4)) = 'C';
                end
            case 6
                if a1==0                   %-- 编码第 1 2位
                    sub_block_DNA(i,((j-1)*4+1)) = 'G';
                elseif a1==1
                    sub_block_DNA(i,((j-1)*4+1)) = 'T';
                elseif a1==2
                    sub_block_DNA(i,((j-1)*4+1)) = 'A';
                else
                    sub_block_DNA(i,((j-1)*4+1)) = 'C';
                end
                if a2==0                   %-- 编码第 3 4位
                    sub_block_DNA(i,((j-1)*4+2)) = 'G';
                elseif a2==1
                    sub_block_DNA(i,((j-1)*4+2)) = 'T';
                elseif a2==2
                    sub_block_DNA(i,((j-1)*4+2)) = 'A';
                else
                    sub_block_DNA(i,((j-1)*4+2)) = 'C';
                end
                if a3==0                   %-- 编码第 5 6位
                    sub_block_DNA(i,((j-1)*4+3)) = 'G';
                elseif a3==1
                    sub_block_DNA(i,((j-1)*4+3)) = 'T';
                elseif a3==2
                    sub_block_DNA(i,((j-1)*4+3)) = 'A';
                else
                    sub_block_DNA(i,((j-1)*4+3)) = 'C';
                end
                if a4==0
                    sub_block_DNA(i,((j-1)*4+4)) = 'G';
                elseif a4==1
                    sub_block_DNA(i,((j-1)*4+4)) = 'T';
                elseif a4==2
                    sub_block_DNA(i,((j-1)*4+4)) = 'A';
                else
                    sub_block_DNA(i,((j-1)*4+4)) = 'C';
                end
            case 7
                if a1==0                   %-- 编码第 1 2位
                    sub_block_DNA(i,((j-1)*4+1)) = 'C';
                elseif a1==1
                    sub_block_DNA(i,((j-1)*4+1)) = 'A';
                elseif a1==2
                    sub_block_DNA(i,((j-1)*4+1)) = 'T';
                else
                    sub_block_DNA(i,((j-1)*4+1)) = 'G';
                end
                if a2==0                   %-- 编码第 3 4位
                    sub_block_DNA(i,((j-1)*4+2)) = 'C';
                elseif a2==1
                    sub_block_DNA(i,((j-1)*4+2)) = 'A';
                elseif a2==2
                    sub_block_DNA(i,((j-1)*4+2)) = 'T';
                else
                    sub_block_DNA(i,((j-1)*4+2)) = 'G';
                end
                if a3==0                   %-- 编码第 5 6位
                    sub_block_DNA(i,((j-1)*4+3)) = 'C';
                elseif a3==1
                    sub_block_DNA(i,((j-1)*4+3)) = 'A';
                elseif a3==2
                    sub_block_DNA(i,((j-1)*4+3)) = 'T';
                else
                    sub_block_DNA(i,((j-1)*4+3)) = 'G';
                end
                if a4==0
                    sub_block_DNA(i,((j-1)*4+4)) = 'C';
                elseif a4==1
                    sub_block_DNA(i,((j-1)*4+4)) = 'A';
                elseif a4==2
                    sub_block_DNA(i,((j-1)*4+4)) = 'T';
                else
                    sub_block_DNA(i,((j-1)*4+4)) = 'G';
                end
            case 8
                if a1==0                   %-- 编码第 1 2位
                    sub_block_DNA(i,((j-1)*4+1)) = 'T';
                elseif a1==1
                    sub_block_DNA(i,((j-1)*4+1)) = 'G';
                elseif a1==2
                    sub_block_DNA(i,((j-1)*4+1)) = 'C';
                else
                    sub_block_DNA(i,((j-1)*4+1)) = 'A';
                end
                if a2==0                   %-- 编码第 3 4位
                    sub_block_DNA(i,((j-1)*4+2)) = 'T';
                elseif a2==1
                    sub_block_DNA(i,((j-1)*4+2)) = 'G';
                elseif a2==2
                    sub_block_DNA(i,((j-1)*4+2)) = 'C';
                else
                    sub_block_DNA(i,((j-1)*4+2)) = 'A';
                end
                if a3==0                   %-- 编码第 5 6位
                    sub_block_DNA(i,((j-1)*4+3)) = 'T';
                elseif a3==1
                    sub_block_DNA(i,((j-1)*4+3)) = 'G';
                elseif a3==2
                    sub_block_DNA(i,((j-1)*4+3)) = 'C';
                else
                    sub_block_DNA(i,((j-1)*4+3)) = 'A';
                end
                if a4==0
                    sub_block_DNA(i,((j-1)*4+4)) = 'T';
                elseif a4==1
                    sub_block_DNA(i,((j-1)*4+4)) = 'G';
                elseif a4==2
                    sub_block_DNA(i,((j-1)*4+4)) = 'C';
                else
                    sub_block_DNA(i,((j-1)*4+4)) = 'A';
                end
        end
    end
end
    