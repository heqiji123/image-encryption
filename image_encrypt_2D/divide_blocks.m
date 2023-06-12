%% 子函数 分块函数，a为每块的边长，I为要分块的图像，num为返回第几大块
function fv=divide_blocks(block_size,image,num)
[~,N]=size(image);
N=N/block_size;
x=floor(num/N)+1;      %第几大行
y=mod(num,N);          %第几大列
if y==0
    x=x-1;
    y=N;
end
fv=image(block_size*(x-1)+1:block_size*x,block_size*(y-1)+1:block_size*y);

