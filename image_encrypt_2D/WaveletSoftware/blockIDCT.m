%  �˺���������ͼ����п�IDCT�任
%  ʱ�䣺2016��9��12��
%  ����ˣ��Ų�
%  ��λ������ͨ��ѧԺ

function reconstructed_image = blockIDCT(DCT_coeff, size_images, block_size)

num_of_block =  size_images / block_size;    %  �ֿ�����������������������Ƿ��棬��������ȣ�

DCT_coeff_cell = mat2cell(DCT_coeff, ones(size_images / block_size,1) * block_size, ones(size_images / block_size,1) * block_size);    %  ����ֿ�
%   mat2cell���Խ�����ֿ죬�������챣���cell���ͣ����Ĳ������оٳ������������

%  �Ը������IDCT�任
for i = 1:num_of_block 
    
    for j = 1:num_of_block 
        
        block_image {i,j} = idct(DCT_coeff_cell{i,j});
        
    end
    
end

reconstructed_image  = cell2mat(block_image);   %  ��cellת���ɾ�����ʽ
end