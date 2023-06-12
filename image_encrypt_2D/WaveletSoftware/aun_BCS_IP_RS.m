%  ���ڼ�Ȩ���Ե�ͼ�����ѹ���㷨

%  ���㷨�Ļ���˼����ͨ����Ȩ���ԣ�����ͼ����м��ܣ�ͬʱ���ѹ������

%  ����ˣ���  ��

%  ʱ�䣺  2018��9��22��

%  ��λ��  ����ͨ��ѧԺ��Ϣ����ϵ


clc;

clear;

%  ����Ҫ�����ݼ���·��

addpath('C:\Users\zb\Desktop\code\images');    


%%  ͼ���ļ�ѡ��
 
filename = 'lenna';                                       %   ͼ���ļ������������ĺô�����ʵ��ʱ���ķ��㣩

% filename = 'peppers';                               %   ͼ���ļ������������ĺô�����ʵ��ʱ���ķ��㣩

% filename = 'barbara';                               %   ͼ���ļ������������ĺô�����ʵ��ʱ���ķ��㣩

% filename = 'goldhill';                                %   ͼ���ļ������������ĺô�����ʵ��ʱ���ķ��㣩

% filename = 'mandrill';                               %   ͼ���ļ������������ĺô�����ʵ��ʱ���ķ��㣩

%%  ������������

size_images =512;                       %   �����������(����ֻ���Ƿ���)

block_size = 32;                           %   �ӿ�Ĵ�С

subrate = 0.3;                              %   Ƿ������

num_trial = 1;                              %   �ظ��������

num_levels = 3;                           %   С���ֽ�ˮƽ����Daubechies 9/7 wavelet transform��

allow = 30;                                   %   ��Ȩ���Ե���ֵ�����ڸ���ֵ���м�Ȩ��

length_of_blocks = block_size^2;                           %    ��ĳ��ȣ�Ԫ�ظ�����

num_of_blocks = (size_images/block_size)^2;      %    �ֿܷ����


%%  ��ȡͼ�񣬲���ͼ�������ָ����С

original_filename = [ filename '.pgm'];                            %    �洢ͼ�����

original_image = double(imread(original_filename));    %    ��ͼ��

original_image = double(imresize(original_image, [size_images,size_images]));     %    ����ͼ���С�������ͼ��X��СΪ�� size_images * size_images

% original_mean = mean(original_image(:));                    %    ��ͼ��ľ�ֵ

% X = original_image - original_mean;                             %    ԭʼͼ���ȥ��ֵ

%%  С���任���ɹ�ѡ�������С���任��Daubechies 9/7 С��������С���任��
%    С���任��Ӧ�����еĲ���1

%  ����С���任

% ww = DWT(size_images);                                                              %    ��������С���任����(N*N)
% 
% wavelet_coeff = ww * original_image * ww';                               %    С��ϵ��������ΪX�������ǶԾ�ֵ������ͼ��任��

% Daubechies 9/7 С���任

wavelet_coeff = waveletcdf97(original_image, num_levels);          %    С���ֽ�
 

%%  �����û�ǰ������ϡ��Ⱥ�ƽ����ϡ���

block_sparsity_wavelet_coeff = blocksparsity (wavelet_coeff, block_size, num_of_blocks, length_of_blocks, allow);     %  ��¼����ϡ���

ave_block_sparsity = sum(block_sparsity_wavelet_coeff)/num_of_blocks;            %    ����ƽ����ϡ���

max_block_sparsity = max( block_sparsity_wavelet_coeff);                                    %     ����ϡ���


%%  �������������ɴ�����PSNRֵ

reconstructed_image_transform = waveletcdf97(wavelet_coeff, -num_levels);              %    С����任

PSNRMAX = psnr(uint8(reconstructed_image_transform),uint8(original_image));        %    �ܹ��ﵽ�����PSNRֵ

 %%  zigzag�û�Ԥ����
 
%  wavelet_coeff = zigzagpermutation(wavelet_coeff, size_images);   


%%  ��Ȩ��������Ȩ������Ϊһ�ּ��ܲ��ԣ�

%  ����1����ͳ��Ȩ����

  wavelet_coeff = reweight(wavelet_coeff);                            %     ��ͳ��Ȩ��������߲���Ч�ʣ�

%  ����2�������������ļ�Ȩ
  
 reweight_matrix = reweight_matrix_fun(size_images, wavelet_coeff, allow);         %    ��Ȩ���󣨼����в��õľ���
%  
 reweight_matrix1 = reweight_matrix_fun(size_images, wavelet_coeff, allow);      %    ��Ȩ���󣨽����м�Ȩ����ƥ��Ľ���Ч����
%   
%   wavelet_coeff =  random_reweight( reweight_matrix, wavelet_coeff);                   %     �����Ȩ


%%   �����û�����

%  �û���������ɣ����ַ�ʽ����֯�û�������û���

%%  ����1������û�����
% 
% permutation1=randperm(size_images);                %   �������������
% 
% permutation2=randperm(size_images);                %   �������������
% 
% P_row = zeros(size_images, size_images);             %   ��ʼ������û�����
% 
% P_column = zeros(size_images, size_images);       %   ��ʼ������û�����
% 
% for i = 1:size_images
%     
%     P_row(i,:) = eyevec(size_images, permutation1(i));
%     
%     P_column(i,:) = eyevec(size_images, permutation2(i));
%     
% end    

%  ����2����֯�û�

IP = ipm(size_images, block_size);                       %  ��֯���������

P_row = IP;

P_column  = IP';

%% ��֯�û�����

wavelet_coeff_ip =P_row * wavelet_coeff * P_column;                          %    �����û�

%%  ѹ����֪����

col_wavelet_coeff_ip = im2col (wavelet_coeff_ip, [block_size block_size], 'distinct');                %   ���ӿ��ų���

for i = 1:num_trial
    
Phi = GenerateProjection(block_size, subrate, filename);           %   ���ɲ�������orth�����ľ���

y = Phi * col_wavelet_coeff_ip ;                                                     %    ����ֵ����

%% ����Ӳ��ֵ�㷨�ع��������ع���

rec_col_wavelet_coeff_ip = zeros(length_of_blocks, num_of_blocks);              %  col_wavelet_coeff_mp�Ĺ���ֵ

for j=1:num_of_blocks
    
    j
    
    if all( y ( :, j) )
        
    rec = separatedecoder( y ( : , j ) , Phi);                                                        %   ��ʱ������IHT�㷨����¼�ع��������У�
    
%     rec = nomp(y(:,j), Phi, length_of_blocks, ave_block_sparsity );             %   ��ʱ������OMP�ع�����¼�ع��������У�
    
    end 
    
         rec_col_wavelet_coeff_ip (: , j)=rec;
         
end

% ����Ӳ��ֵ�㷨��OMP�ع��������ع���  END

%%   ��������

% �е�ͼ��

rec_wavelet_coeff_ip = col2im(rec_col_wavelet_coeff_ip, [block_size block_size], [size_images, size_images], 'distinct');          %  �������ųɾ���ÿ���ų�m*m�ľ����ųɺ�ľ���ά��Ϊ[m,N*N/m]��


%% ���任�ָ�ԭͼ�񣬼���PSNRֵ

rec_wavelet_coeff = P_row' * rec_wavelet_coeff_ip * P_column';                     %    �������û�

rec_wavelet_coeff = unreweight(rec_wavelet_coeff);                                   %   ��ͨ��Ȩ����任

% rec_wavelet_coeff  = random_unreweight(reweight_matrix, rec_wavelet_coeff);                   %   �����Ȩ����任
 
% reconstructed_image = ww' * rec_wavelet_coeff * ww;                           %   С�����任(��������С���任)

reconstructed_image = waveletcdf97(rec_wavelet_coeff, -num_levels);      %  С�����任��Daubechies 9/7 С����

% reconstructed_image = reconstructed_image + original_mean;             %  ԭʼͼ���ȥ��ֵ

reconstructed_image_wiener = wiener2(reconstructed_image, [3 3]);         %  ά���˲����ͼ����󣨵�ѹ���ʱȽ�Сʱ������ά���˲��ɸ������ܣ�

PSNR (i)= psnr(uint8(reconstructed_image), uint8(original_image));                                   %   û���˲��ķ�ֵ�����

PSNR_wiener (i)= psnr(uint8(reconstructed_image_wiener), uint8(original_image));          %   �˲���ķ�ֵ�����
end
 
 PSNR_AV =  sum(PSNR)/num_trial
 
PSNR_wiener_AV =  sum(PSNR_wiener)/num_trial

%    
% % ��ȡͼ��ֲ�
% 
% for i=1:1:256
%     
%      for j=1:1:256
%          
%        X_detail (i, j) =  original_image(i+256, j);
%        
%      end
%      
% end


   
%  ��ʾͼ��


   figure(1);
   
   imshow(uint8(original_image),'Border','tight');

   figure(2);
   
   imshow(uint8(reconstructed_image),'Border','tight');
%    
%    figure(3);
%    
%    imshow(uint8( X_detail),'Border','tight');








