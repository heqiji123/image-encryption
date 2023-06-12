%  基于加权策略的图像加密压缩算法

%  该算法的基本思想是通过加权策略，来对图像进行加密，同时提高压缩性能

%  编程人：张  波

%  时间：  2018年9月22日

%  单位：  重庆通信学院信息工程系


clc;

clear;

%  将需要的数据加入路径

addpath('C:\Users\zb\Desktop\code\images');    


%%  图像文件选择
 
filename = 'lenna';                                       %   图像文件名（这样做的好处在于实验时更改方便）

% filename = 'peppers';                               %   图像文件名（这样做的好处在于实验时更改方便）

% filename = 'barbara';                               %   图像文件名（这样做的好处在于实验时更改方便）

% filename = 'goldhill';                                %   图像文件名（这样做的好处在于实验时更改方便）

% filename = 'mandrill';                               %   图像文件名（这样做的好处在于实验时更改方便）

%%  参数设置区域

size_images =512;                       %   矩阵的行列数(本文只考虑方阵)

block_size = 32;                           %   子块的大小

subrate = 0.3;                              %   欠采样率

num_trial = 1;                              %   重复试验次数

num_levels = 3;                           %   小波分解水平数（Daubechies 9/7 wavelet transform）

allow = 30;                                   %   加权策略的阈值（大于该阈值进行加权）

length_of_blocks = block_size^2;                           %    块的长度（元素个数）

num_of_blocks = (size_images/block_size)^2;      %    总分块个数


%%  读取图像，并将图像调整至指定大小

original_filename = [ filename '.pgm'];                            %    存储图像变量

original_image = double(imread(original_filename));    %    读图像

original_image = double(imresize(original_image, [size_images,size_images]));     %    调整图像大小调整后的图像X大小为： size_images * size_images

% original_mean = mean(original_image(:));                    %    求图像的均值

% X = original_image - original_mean;                             %    原始图像减去均值

%%  小波变换（可供选择的两种小波变换：Daubechies 9/7 小波和正交小波变换）
%    小波变换对应论文中的步骤1

%  正交小波变换

% ww = DWT(size_images);                                                              %    构建正交小波变换矩阵(N*N)
% 
% wavelet_coeff = ww * original_image * ww';                               %    小波系数矩阵（若为X，表明是对均值处理后的图像变换）

% Daubechies 9/7 小波变换

wavelet_coeff = waveletcdf97(original_image, num_levels);          %    小波分解
 

%%  计算置换前的最大块稀疏度和平均块稀疏度

block_sparsity_wavelet_coeff = blocksparsity (wavelet_coeff, block_size, num_of_blocks, length_of_blocks, allow);     %  记录各块稀疏度

ave_block_sparsity = sum(block_sparsity_wavelet_coeff)/num_of_blocks;            %    计算平均块稀疏度

max_block_sparsity = max( block_sparsity_wavelet_coeff);                                    %     最大块稀疏度


%%  计算置零操作后可达的最大PSNR值

reconstructed_image_transform = waveletcdf97(wavelet_coeff, -num_levels);              %    小波逆变换

PSNRMAX = psnr(uint8(reconstructed_image_transform),uint8(original_image));        %    能够达到的最大PSNR值

 %%  zigzag置换预处理
 
%  wavelet_coeff = zigzagpermutation(wavelet_coeff, size_images);   


%%  加权采样（加权策略作为一种加密策略）

%  策略1：传统加权策略

  wavelet_coeff = reweight(wavelet_coeff);                            %     传统加权（仅仅提高采样效率）

%  策略2：带加密特征的加权
  
 reweight_matrix = reweight_matrix_fun(size_images, wavelet_coeff, allow);         %    加权矩阵（加密中采用的矩阵）
%  
 reweight_matrix1 = reweight_matrix_fun(size_images, wavelet_coeff, allow);      %    加权矩阵（解密中加权矩阵不匹配的解密效果）
%   
%   wavelet_coeff =  random_reweight( reweight_matrix, wavelet_coeff);                   %     随机加权


%%   矩阵置换操作

%  置换矩阵的生成（两种方式：交织置换和随机置换）

%%  策略1：随机置换策略
% 
% permutation1=randperm(size_images);                %   将数字随机排列
% 
% permutation2=randperm(size_images);                %   将数字随机排列
% 
% P_row = zeros(size_images, size_images);             %   初始化随机置换矩阵
% 
% P_column = zeros(size_images, size_images);       %   初始化随机置换矩阵
% 
% for i = 1:size_images
%     
%     P_row(i,:) = eyevec(size_images, permutation1(i));
%     
%     P_column(i,:) = eyevec(size_images, permutation2(i));
%     
% end    

%  策略2：交织置换

IP = ipm(size_images, block_size);                       %  交织矩阵的生成

P_row = IP;

P_column  = IP';

%% 交织置换操作

wavelet_coeff_ip =P_row * wavelet_coeff * P_column;                          %    矩阵置换

%%  压缩感知测量

col_wavelet_coeff_ip = im2col (wavelet_coeff_ip, [block_size block_size], 'distinct');                %   将子块排成列

for i = 1:num_trial
    
Phi = GenerateProjection(block_size, subrate, filename);           %   生成测量矩阵（orth处理后的矩阵）

y = Phi * col_wavelet_coeff_ip ;                                                     %    测量值计算

%% 迭代硬阈值算法重构（单独重构）

rec_col_wavelet_coeff_ip = zeros(length_of_blocks, num_of_blocks);              %  col_wavelet_coeff_mp的估计值

for j=1:num_of_blocks
    
    j
    
    if all( y ( :, j) )
        
    rec = separatedecoder( y ( : , j ) , Phi);                                                        %   临时变量（IHT算法，记录重构出来的列）
    
%     rec = nomp(y(:,j), Phi, length_of_blocks, ave_block_sparsity );             %   临时变量（OMP重构，记录重构出来的列）
    
    end 
    
         rec_col_wavelet_coeff_ip (: , j)=rec;
         
end

% 迭代硬阈值算法或OMP重构（单独重构）  END

%%   重新排列

% 列到图像

rec_wavelet_coeff_ip = col2im(rec_col_wavelet_coeff_ip, [block_size block_size], [size_images, size_images], 'distinct');          %  将列重排成矩阵（每列排成m*m的矩阵，排成后的矩阵维度为[m,N*N/m]）


%% 反变换恢复原图像，计算PSNR值

rec_wavelet_coeff = P_row' * rec_wavelet_coeff_ip * P_column';                     %    矩阵逆置换

rec_wavelet_coeff = unreweight(rec_wavelet_coeff);                                   %   普通加权的逆变换

% rec_wavelet_coeff  = random_unreweight(reweight_matrix, rec_wavelet_coeff);                   %   随机加权的逆变换
 
% reconstructed_image = ww' * rec_wavelet_coeff * ww;                           %   小波反变换(采用正交小波变换)

reconstructed_image = waveletcdf97(rec_wavelet_coeff, -num_levels);      %  小波反变换（Daubechies 9/7 小波）

% reconstructed_image = reconstructed_image + original_mean;             %  原始图像减去均值

reconstructed_image_wiener = wiener2(reconstructed_image, [3 3]);         %  维纳滤波后的图像矩阵（当压缩率比较小时，采用维纳滤波可改善性能）

PSNR (i)= psnr(uint8(reconstructed_image), uint8(original_image));                                   %   没有滤波的峰值信噪比

PSNR_wiener (i)= psnr(uint8(reconstructed_image_wiener), uint8(original_image));          %   滤波后的峰值信噪比
end
 
 PSNR_AV =  sum(PSNR)/num_trial
 
PSNR_wiener_AV =  sum(PSNR_wiener)/num_trial

%    
% % 截取图像局部
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


   
%  显示图像


   figure(1);
   
   imshow(uint8(original_image),'Border','tight');

   figure(2);
   
   imshow(uint8(reconstructed_image),'Border','tight');
%    
%    figure(3);
%    
%    imshow(uint8( X_detail),'Border','tight');








