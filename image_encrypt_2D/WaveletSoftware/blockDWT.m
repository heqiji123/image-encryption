%  �˺���������ͼ����п�DCT�任
%  ʱ�䣺2016��9��12��
%  ����ˣ��Ų�
%  ��λ������ͨ��ѧԺ


function DCT_Coeff = blockDWT(original_image, size_images, block_size, num_levels)

num_of_block =  size_images / block_size;    %  �ֿ�����������������������Ƿ��棬��������ȣ�

block_original_image = mat2cell(original_image, ones(size_images / block_size,1) * block_size, ones(size_images / block_size,1) * block_size);       %  ����ֿ�
%   mat2cell���Խ�����ֿ죬�������챣���cell���ͣ����Ĳ������оٳ������������

%  �Ը������DCT�任
for i = 1:num_of_block 
    
    for j = 1:num_of_block 
        
        block_DCT_coeff_cell{i,j} = waveletcdf97(block_original_image{i,j}, num_levels);
        
    end
    
end

DCT_Coeff = cell2mat(block_DCT_coeff_cell);   %  ��cellת���ɾ�����ʽ
end