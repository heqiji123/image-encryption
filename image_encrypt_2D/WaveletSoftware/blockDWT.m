%  此函数用来对图像进行块DCT变换
%  时间：2016年9月12日
%  编程人：张波
%  单位：重庆通信学院


function DCT_Coeff = blockDWT(original_image, size_images, block_size, num_levels)

num_of_block =  size_images / block_size;    %  分块个数（行列数，由于这里是仿真，行列数相等）

block_original_image = mat2cell(original_image, ones(size_images / block_size,1) * block_size, ones(size_images / block_size,1) * block_size);       %  矩阵分块
%   mat2cell可以将矩阵分快，并将个快保存成cell类型，它的参数是列举出各块的行列数

%  对各块进行DCT变换
for i = 1:num_of_block 
    
    for j = 1:num_of_block 
        
        block_DCT_coeff_cell{i,j} = waveletcdf97(block_original_image{i,j}, num_levels);
        
    end
    
end

DCT_Coeff = cell2mat(block_DCT_coeff_cell);   %  将cell转换成矩阵形式
end