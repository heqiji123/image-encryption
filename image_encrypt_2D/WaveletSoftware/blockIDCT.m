%  此函数用来对图像进行块IDCT变换
%  时间：2016年9月12日
%  编程人：张波
%  单位：重庆通信学院

function reconstructed_image = blockIDCT(DCT_coeff, size_images, block_size)

num_of_block =  size_images / block_size;    %  分块个数（行列数，由于这里是仿真，行列数相等）

DCT_coeff_cell = mat2cell(DCT_coeff, ones(size_images / block_size,1) * block_size, ones(size_images / block_size,1) * block_size);    %  矩阵分块
%   mat2cell可以将矩阵分快，并将个快保存成cell类型，它的参数是列举出各块的行列数

%  对各块进行IDCT变换
for i = 1:num_of_block 
    
    for j = 1:num_of_block 
        
        block_image {i,j} = idct(DCT_coeff_cell{i,j});
        
    end
    
end

reconstructed_image  = cell2mat(block_image);   %  将cell转换成矩阵形式
end