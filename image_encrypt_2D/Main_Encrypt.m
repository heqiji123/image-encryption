%% 基于BCS，DNA的分布式灰度图像加密算法
clc;
clear;

% plain_image = double(imread('9.bmp'));
% plain_image = double(imread('2.tiff'));
% plain_image = double(imread('Lenna24_gray.bmp'));
% plain_image = double(imread('lenna.pgm','pgm'));
plain_image = double(imread('Lenna.tif','tif'));
% plain_image = double(imread('airplane.bmp','bmp'));
% plain_image = double(imread('couple.bmp','bmp'));
% plain_image = double(imread('fishingboat.bmp','bmp'));
% plain_image = double(imread('watch.bmp','bmp'));
figure;imshow(uint8(plain_image));title('原图片');
[~,N] = size(plain_image);

%-- 原始图像SHA_384值
meth = 'SHA-384';
K = hash(plain_image,meth);
K = reshape(K,48,8);

k = zeros(1,48);
for i=1:48
    k(i) = bin2dec(K(i,:));
end

%% 计算原始图像的三个通道的信息熵
%原始图片 
%R通道
T1=imhist(uint8(plain_image));                  %统计图像R通道灰度值从0~255的分布情况，存至T1
S1=sum(T1);                                     %计算整幅图像R通道的灰度值
entropy=0;                                      %原始图片R通道相关性

for i=1:length(T1)
    pp1=T1(i)/S1;                               %每个灰度值占比，即每个灰度值的概率
    if pp1~=0
        entropy=entropy-pp1*log2(pp1);          %求信息熵
    end
end


%% 计算四维混沌系统的初始值 x0 y0 z0 w0

%-- 密钥
t1 = 0.86;
t2 = 0.43;
t3 = 0.23;
t4 = 0.39;

sum_1 = k(1);
for i=2:6
    sum_1 = bitxor(sum_1,k(i));
end
x0 = (t1 * sum_1) /256;

sum_2 = k(7);
for i=8:12
    sum_2 = bitxor(sum_2,k(i));
end
y0 = (t2 * sum_2) / 256;

sum_3 = k(13);
for i=14:18
    sum_3 = bitxor(sum_3,k(i));
end
z0 = (t3 * sum_3) / 256;

sum_4 = k(19);
for i=20:24
    sum_4 = bitxor(sum_4,k(i));
end
w0 = (t4 * sum_4) / 256;


%% 计算LTS混沌系统的初始值 x01 r01 
% 初始参数  作为密钥
x01_ = 0.1;
r01_ = 2.5;
x02_ = 0.3;
r02_ = 3.3;


sum_1 = k(25)+k(26)+k(27)+k(28)+k(29)+k(30);
sum_2 = k(31)+k(32)+k(33)+k(34)+k(35)+k(36);
x01 = mod(((sum_1/256)+x01_),1);
r01 = mod(((sum_2/256)+r01_),4);


%% 计算LSS混沌系统的初始值 x02 r02
sum_3 = k(37);
for i=38:42
    sum_3 = bitxor(sum_3,k(i));
end
sum_4 = k(43);
for i=44:48
    sum_4 = bitxor(sum_4,k(i));
end

x02 = mod((sum_3/256)+x02_,1);
r02 = mod((sum_4/256)+r02_,4);



%% 生成测量矩阵

CR = 0.25;                                        %   压缩比率
image_size = N;                                   %   the number of rows (or columns) of the image ( square matrix only! )
block_size = 32;                                  %   block size
num_levels = 3;                                   %   Wavelet decomposition level    
max_iterations = 200;                             %   the maximum iteration number
length_of_blocks = block_size * block_size;       %   the length of block 
block_nums = (N/block_size) * (N/block_size);     %   块的数量

nums_col = length_of_blocks;                      %   分块后图像大小---row
nums_row = round(CR * nums_col);

Phi = create_MM(nums_row,nums_col,x01,r01);       %   LTS混沌序列生成测量矩阵

    

%% 压缩加密图像

image_blocked = im2col(plain_image, [block_size, block_size], 'distinct'); 
image_compress = Phi * image_blocked;



%% 量化压缩图像

image_compress_min = min(min(image_compress));
image_compress_max = max(max(image_compress));

image_compress_q = zeros(nums_row,block_nums);
for i=1:nums_row*block_nums
	image_compress_q(i) = round(255*((image_compress(i)-image_compress_min)/(image_compress_max-image_compress_min)));
end




%% 置乱——扩散(DNA)


% -- 生成用于DNA置乱扩散的四翼混沌序列X_、Y_、Z_、W_
L0 = 1500;
KK = -14;   %KK作为密钥之一
d = 10;
iterations = L0 + block_nums*nums_row*d;
iterate = Chaos_output(x0,y0,z0,w0,iterations,KK);
X = iterate(:,1);
X = X(L0+2:length(X));
Y = iterate(:,2);
Y = Y(L0+2:length(Y));
Z = iterate(:,3);
Z = Z(L0+2:length(Z));
W = iterate(:,4);
W = W(L0+2:length(W));

X_ = zeros(block_nums,nums_row);
Y_ = zeros(block_nums,nums_row);
Z_ = zeros(block_nums,nums_row);
W_ = zeros(block_nums,nums_row);
for i=1:block_nums*nums_row
    X_(i) = X(d*(i-1)+1);
    Y_(i) = Y(d*(i-1)+1);
    Z_(i) = Z(d*(i-1)+1);
    W_(i) = W(d*(i-1)+1);
end


% -- LSS生成用于DNA_encode DNA_decode的混沌序列
N0 = 1000;
Lss_DNA_encode = LSS(x02,r02,N0,block_nums);     %-- DNA_encode

x03 = mod((entropy)/256 + x02_,1);
r03 = mod((entropy)/256 + r02_,4);

Lss_DNA_decode = LSS(x03,r03,N0,block_nums);     %-- DNA_decode

for i=1:block_nums
    Lss_DNA_encode(i) = mod(round(Lss_DNA_encode(i)*1000),8) + 1;
    Lss_DNA_decode(i) = mod(round(Lss_DNA_decode(i)*1000),8) + 1;
end


image_encrypt = zeros(nums_row,block_nums);
for i=1:block_nums
    
    %-- 分块
    sub_block = reshape(image_compress_q(:,i),1,nums_row);
    sub_block_p = index_permutation(sub_block,X_(i,:));           %-- 基于混沌序列索引置乱


    %-- 基于混沌序列Lss_DNA_encode进行DNA编码,基于PZ进行DNA运算，基于W_进行DNA解码
    %-- DNA编码
    sub_block_DNA = DNA_encode_(sub_block_p,Lss_DNA_encode(i));

    
    %-- DNA运算
    PZ = [X_(i,:),Y_(i,:),Z_(i,:),W_(i,:)];
    for m=1:length(PZ)
            PZ(m) = mod(round(PZ(m)*1000),4)+1;
    end
    sub_block_DNA_C = DNA_compute(sub_block_DNA,PZ);

    
    %-- DNA解码
    sub_block = DNA_decode_(sub_block_DNA_C,Lss_DNA_decode(i));
    
    image_encrypt(:,i) = reshape(sub_block,nums_row,1);
end



%% 解密图像


image_decrypt = zeros(nums_row,block_nums);
for i=1:block_nums

    d_sub_block = reshape(image_encrypt(:,i),1,nums_row);
    
    %-- DNA编码
    d_sub_block_DNA = DNA_encode_(d_sub_block,Lss_DNA_decode(i));
    
    %-- DNA逆运算
    PZ_ic = [X_(i,:),Y_(i,:),Z_(i,:),W_(i,:)];
    for m=1:length(PZ_ic)
    	PZ_ic(m) = mod(round(PZ_ic(m)*1000),4)+1;
    end
    d_sub_block_DNA_ic = DNA_compute_inverse(d_sub_block_DNA,PZ_ic);
    
    %-- DNA解码
    d_sub_block_p = DNA_decode_(d_sub_block_DNA_ic,Lss_DNA_encode(i));
    
    %-- 逆置乱
    d_sub_block = index_permutation_inverse(d_sub_block_p,X_(i,:));
    
  
    image_decrypt(:,i) = reshape(d_sub_block,nums_row,1);
    
end


%% 逆量化解密图像


I_compress_qiv = zeros(nums_row,block_nums);
for i=1:nums_row*block_nums
	I_compress_qiv(i) = image_decrypt(i)*(image_compress_max - image_compress_min)/255 + image_compress_min;
end



%%  重构图像


image_reconstructed = BCS_ED(I_compress_qiv, Phi, image_size, image_size, num_levels );  
    
figure;imshow(uint8(image_encrypt));title('加密后的图像');
figure;imhist(uint8(image_encrypt));title('加密后的图像的直方图');
figure;imshow(uint8(image_reconstructed));title('重构的图像');

MSE = 0;
for i=1:N
    for j=1:N
        MSE = MSE+(double(plain_image(i,j))-double(image_reconstructed(i,j)))^2;
    end
end
MSE = MSE/(N*N);
PSNR = 10*log10((255*255)/MSE);  

% %% 嵌入过程
% D1_R = zeros(M,N);
% D2_R = zeros(M,N);
% D1_G = zeros(M,N);
% D2_G = zeros(M,N);
% D1_B = zeros(M,N);
% D2_B = zeros(M,N);
% for i=1:M
%     for j=1:N
%         D1_R(i,j) = mod(R4(i,j),10);
%         D2_R(i,j) = floor(R4(i,j)/10);
%         D1_G(i,j) = mod(G4(i,j),10);
%         D2_G(i,j) = floor(G4(i,j)/10);
%         D1_B(i,j) = mod(B4(i,j),10);
%         D2_B(i,j) = floor(B4(i,j)/10);
%     end
% end
% 
% D1_R = reshape(D1_R,N/2,N/2);
% D2_R = reshape(D2_R,N/2,N/2);
% D1_G = reshape(D1_G,N/2,N/2);
% D2_G = reshape(D2_G,N/2,N/2);
% D1_B = reshape(D1_B,N/2,N/2);
% D2_B = reshape(D2_B,N/2,N/2);
% 
% % 载体图片
% carrier_image = double(imread('House_carrier.png'));
% CA_R = carrier_image(:,:,1);
% CA_G = carrier_image(:,:,2);
% CA_B = carrier_image(:,:,3);
% 
% %嵌入D1_R D2_R D1_G D2_G D1_B D2_B
% [CA_R_A,CA_R_H,CA_R_V,CA_R_D] = IWT(CA_R);
% CA_R_A_mean = mean(CA_R_A,'all');
% Tag_R = zeros(N/2,N/2);
% [CA_G_A,CA_G_H,CA_G_V,CA_G_D] = IWT(CA_G);
% CA_G_A_mean = mean(CA_G_A,'all');
% Tag_G = zeros(N/2,N/2);
% [CA_B_A,CA_B_H,CA_B_V,CA_B_D] = IWT(CA_B);
% CA_B_A_mean = mean(CA_B_A,'all');
% Tag_B = zeros(N/2,N/2);
% for i=1:N/2
%     for j=1:N/2
%         if CA_R_A(i,j) >= CA_R_A_mean
%             Tag_R(i,j) = 1;
%             CA_R_H(i,j) = D1_R(i,j);
%             CA_R_V(i,j) = D2_R(i,j);
%         else
%             Tag_R(i,j) = 0;
%             CA_R_H(i,j) = D2_R(i,j);
%             CA_R_V(i,j) = D1_R(i,j);
%         end
%         if CA_G_A(i,j) >= CA_G_A_mean
%             Tag_G(i,j) = 1;
%             CA_G_H(i,j) = D1_G(i,j);
%             CA_G_V(i,j) = D2_G(i,j);
%         else
%             Tag_G(i,j) = 0;
%             CA_G_H(i,j) = D2_G(i,j);
%             CA_G_V(i,j) = D1_G(i,j);
%         end
%         if CA_B_A(i,j) >= CA_B_A_mean
%             Tag_B(i,j) = 1;
%             CA_B_H(i,j) = D1_B(i,j);
%             CA_B_V(i,j) = D2_B(i,j);
%         else
%             Tag_B(i,j) = 0;
%             CA_B_H(i,j) = D2_B(i,j);
%             CA_B_V(i,j) = D1_B(i,j);
%         end
%     end
% end
% 
% CA_R_carrier = IIWT(CA_R_A,CA_R_H,CA_R_V,CA_R_D);
% CA_G_carrier = IIWT(CA_G_A,CA_G_H,CA_G_V,CA_G_D);
% CA_B_carrier = IIWT(CA_B_A,CA_B_H,CA_B_V,CA_B_D);
% VM_image = zeros(N,N,3);
% VM_image(:,:,1) = CA_R_carrier;
% VM_image(:,:,2) = CA_G_carrier;
% VM_image(:,:,3) = CA_B_carrier;

% %-- 保存图像
% VM_image = uint8(VM_image);
% figure;imshow(VM_image);title('ciper image');
% imwrite(VM_image,'ciper_image.png','png');
% %-- 保存参数
% save('secret_key','k','t1','t2','t3','t4','t5','TS','d','CR','R_max','R_min','G_max','G_min','B_max','B_min','x01_','r01_','x02_','r02_','x03_','r03_');
% 
% image_encrypt = zeros(M,N,3);
% image_encrypt(:,:,1) = R8;
% image_encrypt(:,:,2) = G8;
% image_encrypt(:,:,3) = B8;
% image_encrypt = uint8(image_encrypt);
% figure;imhist(image_encrypt);



