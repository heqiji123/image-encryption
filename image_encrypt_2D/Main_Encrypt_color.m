%% 基于2DCS，DNA的分布式彩色图像加密算法
clc;
clear;

addpath('C:\Users\Thermos_121380\matlab\paper_code\2D compressive sensing with a multi-embedding strategy\BCS_image_encrypt\WaveletSoftware'); %  WaveletSoftware
tic;
plain_image = double(imread('peppers.bmp','bmp'));
imwrite(uint8(plain_image),'plain_image.bmp');
figure;imshow(uint8(plain_image));title('原图片');
[~,N,~] = size(plain_image);


%-- 原始图像SHA_384值
meth = 'SHA-384';
K = hash(plain_image,meth);
K = reshape(K,48,8);


k = zeros(1,48);
for i=1:48
    k(i) = bin2dec(K(i,:));
end

%--- 得到哈希值的十六进制
tmp = dec2hex(k);
[rows,~] = size(tmp);
res = char(1,2*rows);
for i = 1:rows
    res(2*i - 1) = tmp(i,1);
    res(2*i) = tmp(i,2);
end
res = reshape(res,1,length(res));

%% 计算原始图像的三个通道的信息熵
%原始图片 
%R通道
T1=imhist(uint8(plain_image(:,:,1)));           %统计图像R通道灰度值从0~255的分布情况，存至T1
S1=sum(T1);                                     %计算整幅图像R通道的灰度值
entropy_R=0;                                    %原始图片R通道相关性

for i=1:length(T1)
    pp1=T1(i)/S1;                               %每个灰度值占比，即每个灰度值的概率
    if pp1~=0
        entropy_R=entropy_R-pp1*log2(pp1);      %求信息熵
    end
end

T2=imhist(uint8(plain_image(:,:,2)));            %统计图像R通道灰度值从0~255的分布情况，存至T1
S2=sum(T2);                                      %计算整幅图像R通道的灰度值
entropy_G=0;                                     %原始图片R通道相关性

for i=1:length(T2)
    pp1=T2(i)/S2;                                %每个灰度值占比，即每个灰度值的概率
    if pp1~=0
        entropy_G=entropy_G-pp1*log2(pp1);       %求信息熵
    end
end


T3=imhist(uint8(plain_image(:,:,3)));            %统计图像R通道灰度值从0~255的分布情况，存至T1
S3=sum(T3);                                      %计算整幅图像R通道的灰度值
entropy_B=0;                                     %原始图片R通道相关性

for i=1:length(T3)
    pp1=T3(i)/S3;                                %每个灰度值占比，即每个灰度值的概率
    if pp1~=0
        entropy_B=entropy_B-pp1*log2(pp1);       %求信息熵
    end
end


%% 计算LTS混沌系统的初始值 x01 r01 
% 初始参数  作为密钥
t1 = 0.25;
t2 = 2.78;
t3 = 0.43;
t4 = 3.72;
t5 = 0.39;
t6 = 1.90;

sum_1 = k(1)+k(2)+k(3)+k(4)+k(5)+k(6);
sum_2 = k(7)+k(8)+k(9)+k(10)+k(11)+k(12);
x01 = mod(((sum_1/256)+t1),1);
r01 = mod(((sum_2/256)+t2),4);                  % -- 用于生成Phi1

x02 = mod((x01+(entropy_R+entropy_G)/2),1);
r02 = mod((r01+(entropy_R+entropy_B)/2),4);     % -- 用于生成Phi2


%% 计算LSS混沌系统的初始值 x03 r03
sum_3 = k(13);
for i=14:18
    sum_3 = bitxor(sum_3,k(i));
end
sum_4 = k(19);
for i=20:24
    sum_4 = bitxor(sum_4,k(i));
end

x03 = mod((sum_3/256)+t3,1);
r03 = mod((sum_4/256)+t4,4);

x04 = mod((x03+entropy_R+entropy_G+entropy_B),1);
r04 = mod((r03+entropy_R+entropy_G+entropy_B),4);



%% 计算TSS混沌系统的初始参数 x05 r05

sum_5 = bitxor(k(25)+k(26)+k(27),k(28)+k(29)+k(30));
x05 = mod((sum_5/256)+t5,1);

sum_6 = bitxor(k(31)+k(32)+k(33),k(34)+k(35)+k(36));
r05 = mod((sum_6/256)+t6,1);



%% Josephus 的参数  ??
x06 = (k(37)+k(38)+k(39)+k(40)+k(41)+k(42)) / 6;
r06 = bitxor(k(43)+k(44)+k(45),k(46)+k(47)+k(48));



%% 生成测量矩阵

CR = 0.25;                                             %   压缩比率
image_size = N;                                        %   the number of rows (or columns) of the image ( square matrix only! )
num_levels = 3;                                        %   Wavelet decomposition level   
M = round(sqrt(CR * image_size * image_size));   
block_size = 8;                                        %   块的大小
block_fill = block_size - mod(M,block_size);
if block_fill ~= block_size
    M = M + block_fill;
end
block_nums = (M/block_size) * (M/block_size);          %   块的数量

Phi1 = create_MM(M,N,x01,r01);                         %   LTS混沌序列生成测量矩阵
Phi2 = create_MM(M,N,x02,r02);

    


%% 分解图像
plain_image_R = plain_image(:,:,1);
plain_image_G = plain_image(:,:,2);
plain_image_B = plain_image(:,:,3);



%% 压缩图像

image_compress_R = Phi1 * plain_image_R * Phi2';
image_compress_G = Phi1 * plain_image_G * Phi2';
image_compress_B = Phi1 * plain_image_B * Phi2';



%% 量化压缩图像

image_compress_R_min = min(min(image_compress_R ));
image_compress_R_max = max(max(image_compress_R ));
image_compress_R_q = zeros(M,M);
image_compress_G_min = min(min(image_compress_G ));
image_compress_G_max = max(max(image_compress_G ));
image_compress_G_q = zeros(M,M);
image_compress_B_min = min(min(image_compress_B ));
image_compress_B_max = max(max(image_compress_B ));
image_compress_B_q = zeros(M,M);

for i=1:M
    for j=1:M
        image_compress_R_q(i,j) = round(255*((image_compress_R(i,j)-image_compress_R_min)/(image_compress_R_max-image_compress_R_min)));
        image_compress_G_q(i,j) = round(255*((image_compress_G(i,j)-image_compress_G_min)/(image_compress_G_max-image_compress_G_min)));
        image_compress_B_q(i,j) = round(255*((image_compress_B(i,j)-image_compress_B_min)/(image_compress_B_max-image_compress_B_min)));
    end
end


%% 置乱——扩散(DNA)

% -- LSS生成用于DNA_encode DNA_decode的混沌序列
N0 = 1000;
d = 10;
Lss_DNA_encode = LSS(x03,r03,N0,block_nums*d);     %-- DNA_encode
Lss_DNA_decode = LSS(x04,r04,N0,block_nums*d);     %-- DNA_decode

for i=1:block_nums
        Lss_DNA_encode(i) = mod(round(Lss_DNA_encode((i-1)*d+1)*1000),8) + 1;
        Lss_DNA_decode(i) = mod(round(Lss_DNA_decode((i-1)*d+1)*1000),8) + 1;  
end
Lss_DNA_encode = Lss_DNA_encode(1:block_nums);
Lss_DNA_decode = Lss_DNA_decode(1:block_nums);

% -- 生成用于DNA置乱扩散的混沌序列T
T = TSS(x05,r05,N0,block_nums*d);
for i=1:block_nums
    T(i) = mod(round(T((i-1)*d+1)*1000),4)+1;
end
T = T(1:block_nums);

% image_encrypt_R = zeros(M,M);
% image_encrypt_G = zeros(M,M);
% image_encrypt_B = zeros(M,M);

S = [image_compress_R_q, image_compress_G_q, image_compress_B_q];
start = mod(round(x06 * 1000),3*M);
step = mod(round(r06*1000),19)+1;
S_ = Josephus_permutation(S,start,step);                %-- 约瑟夫问题置乱  使得R、G、B三个分量相互扩散
image_compress_R_q_s = S_(:,1:M);
image_compress_G_q_s = S_(:,M+1:2*M);
image_compress_B_q_s = S_(:,2*M+1:3*M);

% -- DNA编码
% R_DNA_1 = DNA_encode_(image_compress_R_q_s, Lss_DNA_encode);
% G_DNA_1 = DNA_encode_(image_compress_G_q_s, Lss_DNA_encode);
% B_DNA_1 = DNA_encode_(image_compress_B_q_s, Lss_DNA_encode);

% -- DNA运算
% R_DNA_C1 = DNA_compute(R_DNA_1, T);
% G_DNA_C1 = DNA_compute(G_DNA_1, T);
% B_DNA_C1 = DNA_compute(B_DNA_1, T);

% -- R、G、B分量间进行加，减，异或
% R_DNA_C2 = DNA_operation(B_DNA_C1,R_DNA_C1,1);    %加法
% G_DNA_C2 = DNA_operation(R_DNA_C2,G_DNA_C1,2);    %减法
% B_DNA_C2 = DNA_operation(G_DNA_C2,B_DNA_C1,3);    %异或

% -- DNA解码
% image_encrypt_R = DNA_decode_(R_DNA_C1,Lss_DNA_decode);
% image_encrypt_G = DNA_decode_(G_DNA_C1,Lss_DNA_decode);
% image_encrypt_B = DNA_decode_(B_DNA_C1,Lss_DNA_decode);

for i=1:block_nums
  
	%-- 上一块
    if i==1
        sub_block_R_last = divide_blocks(block_size,image_compress_R_q_s,block_nums);
        sub_block_G_last = divide_blocks(block_size,image_compress_G_q_s,block_nums);
        sub_block_B_last = divide_blocks(block_size,image_compress_B_q_s,block_nums);

        sub_block_R_DNA_last = DNA_encode_(sub_block_R_last,Lss_DNA_encode(block_nums));
        sub_block_G_DNA_last = DNA_encode_(sub_block_G_last,Lss_DNA_encode(block_nums));
        sub_block_B_DNA_last = DNA_encode_(sub_block_B_last,Lss_DNA_encode(block_nums));
    else
        sub_block_R_last = divide_blocks(block_size,image_encrypt_R,i-1);
        sub_block_G_last = divide_blocks(block_size,image_encrypt_G,i-1);
        sub_block_B_last = divide_blocks(block_size,image_encrypt_B,i-1);

        sub_block_R_DNA_last = DNA_encode_(sub_block_R_last,Lss_DNA_decode(i-1));
        sub_block_G_DNA_last = DNA_encode_(sub_block_G_last,Lss_DNA_decode(i-1));
        sub_block_B_DNA_last = DNA_encode_(sub_block_B_last,Lss_DNA_decode(i-1));
    end
    
    %-- 分块
    sub_block_R = divide_blocks(block_size,image_compress_R_q_s,i);
    sub_block_G = divide_blocks(block_size,image_compress_G_q_s,i);
    sub_block_B = divide_blocks(block_size,image_compress_B_q_s,i);

    %-- 基于混沌序列Lss_DNA_encode进行DNA编码,基于PZ进行DNA运算，基于Lss_DNA_decode进行DNA解码
    %-- DNA编码
    sub_block_R_DNA = DNA_encode_(sub_block_R,Lss_DNA_encode(i));
    sub_block_G_DNA = DNA_encode_(sub_block_G,Lss_DNA_encode(i));
    sub_block_B_DNA = DNA_encode_(sub_block_B,Lss_DNA_encode(i));

    %-- 与上一块进行相加
    sub_block_R_DNA_C1 = DNA_operation(sub_block_R_DNA,sub_block_R_DNA_last,1);
    sub_block_G_DNA_C1 = DNA_operation(sub_block_G_DNA,sub_block_G_DNA_last,1);
    sub_block_B_DNA_C1 = DNA_operation(sub_block_B_DNA,sub_block_B_DNA_last,1);
    
    %-- DNA运算
    sub_block_R_DNA_C2 = DNA_compute(sub_block_R_DNA_C1,T(i));
    sub_block_G_DNA_C2 = DNA_compute(sub_block_G_DNA_C1,T(i));
    sub_block_B_DNA_C2 = DNA_compute(sub_block_B_DNA_C1,T(i));

    
    %-- DNA解码
    sub_block_R = DNA_decode_(sub_block_R_DNA_C2,Lss_DNA_decode(i));
    sub_block_G = DNA_decode_(sub_block_G_DNA_C2,Lss_DNA_decode(i));
    sub_block_B = DNA_decode_(sub_block_B_DNA_C2,Lss_DNA_decode(i));
    
    block_q = M / block_size;
    xx=floor(i/block_q)+1;      %第几大行
    yy=mod(i,block_q);          %第几大列
    if yy==0
        xx=xx-1;
        yy=block_q;
    end
    image_encrypt_R(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = sub_block_R;
    image_encrypt_G(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = sub_block_G;
    image_encrypt_B(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = sub_block_B;
end

%%%%%%%%%%%%%%%%%%%%%%%

% for i=1:block_nums
%   
%     %-- 分块
%     sub_block_R = divide_blocks(block_size,image_compress_R_q_s,i);
%     sub_block_G = divide_blocks(block_size,image_compress_G_q_s,i);
%     sub_block_B = divide_blocks(block_size,image_compress_B_q_s,i);
% 
%     %-- DNA编码
%     sub_block_R_DNA = DNA_encode_(sub_block_R,Lss_DNA_encode(i));
%     sub_block_G_DNA = DNA_encode_(sub_block_G,Lss_DNA_encode(i));
%     sub_block_B_DNA = DNA_encode_(sub_block_B,Lss_DNA_encode(i));
% 
%     %-- 三个分量之间进行异或
%     sub_block_R_DNA_C1 = DNA_operation(sub_block_B_DNA,sub_block_R_DNA,3);
%     sub_block_G_DNA_C1 = DNA_operation(sub_block_R_DNA_C1,sub_block_G_DNA,3);
%     sub_block_B_DNA_C1 = DNA_operation(sub_block_G_DNA_C1,sub_block_B_DNA,3);
%     
%     %-- DNA运算
%     sub_block_R_DNA_C2 = DNA_compute(sub_block_R_DNA_C1,T(i));
%     sub_block_G_DNA_C2 = DNA_compute(sub_block_G_DNA_C1,T(i));
%     sub_block_B_DNA_C2 = DNA_compute(sub_block_B_DNA_C1,T(i));
% 
%     
%     %-- DNA解码
%     sub_block_R = DNA_decode_(sub_block_R_DNA_C2,Lss_DNA_decode(i));
%     sub_block_G = DNA_decode_(sub_block_G_DNA_C2,Lss_DNA_decode(i));
%     sub_block_B = DNA_decode_(sub_block_B_DNA_C2,Lss_DNA_decode(i));
%     
%     block_q = M / block_size;
%     xx=floor(i/block_q)+1;      %第几大行
%     yy=mod(i,block_q);          %第几大列
%     if yy==0
%         xx=xx-1;
%         yy=block_q;
%     end
%     image_encrypt_R(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = sub_block_R;
%     image_encrypt_G(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = sub_block_G;
%     image_encrypt_B(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = sub_block_B;
% end

%%%%%%%%%%%%%%%%%%%%%%%


image_encrypt = zeros(M,M,3);
image_encrypt(:,:,1) = image_encrypt_R;
image_encrypt(:,:,2) = image_encrypt_G;
image_encrypt(:,:,3) = image_encrypt_B;


%% 嵌入图像
% D1_R = zeros(M,M);
% D2_R = zeros(M,M);
% D1_G = zeros(M,M);
% D2_G = zeros(M,M);
% D1_B = zeros(M,M);
% D2_B = zeros(M,M);
% for i=1:M
%     for j=1:M
%         D1_R(i,j) = mod(image_encrypt_R(i,j),10);
%         D2_R(i,j) = floor(image_encrypt_R(i,j)/10);
%         D1_G(i,j) = mod(image_encrypt_G(i,j),10);
%         D2_G(i,j) = floor(image_encrypt_G(i,j)/10);
%         D1_B(i,j) = mod(image_encrypt_B(i,j),10);
%         D2_B(i,j) = floor(image_encrypt_B(i,j)/10);
%     end
% end
% 
% 
% % 载体图片
% carrier_image = double(imread('Tiffany24.bmp'));
% imwrite(uint8(carrier_image),'carrier_image.bmp');
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
% 
% %-- 保存图像
% figure;imshow(uint8(VM_image));title('ciper image');
% imwrite(uint8(VM_image),'ciper_image.bmp','bmp');



%%  噪声攻击分析  阻塞攻击分析


%-- 椒盐噪声、乘性噪声、高斯噪声
% image_encrypt_R = double(imnoise(uint8(image_encrypt_R),'salt & pepper',0.0001));
% image_encrypt_G = double(imnoise(uint8(image_encrypt_G),'salt & pepper',0.0001));
% image_encrypt_B = double(imnoise(uint8(image_encrypt_B),'salt & pepper',0.0001));

% image_encrypt_R = double(imnoise(uint8(image_encrypt_R),'speckle',0.000001));
% image_encrypt_G = double(imnoise(uint8(image_encrypt_G),'speckle',0.000001));
% image_encrypt_B = double(imnoise(uint8(image_encrypt_B),'speckle',0.000001));

% image_encrypt_R = double(imnoise(uint8(image_encrypt_R),'gaussian',0,0.000001));
% image_encrypt_G = double(imnoise(uint8(image_encrypt_G),'gaussian',0,0.000001));
% image_encrypt_B = double(imnoise(uint8(image_encrypt_B),'gaussian',0,0.000001));

% -- 阻塞攻击             PSNR = 24.7662
% -- 16 * 16 -- 中心
% image_encrypt_R(121:136,121:136) = 0;
% image_encrypt_G(121:136,121:136) = 0;
% image_encrypt_B(121:136,121:136) = 0;

% -- 16 * 16 -- 四角
% image_encrypt_R(1:8,1:8) = 0;            %-- 左上
% image_encrypt_R(1:8,249:256) = 0;        %-- 右上
% image_encrypt_R(249:256,1:8) = 0;        %-- 左下
% image_encrypt_R(249:256,249:256) = 0;    %-- 右下
% 
% image_encrypt_G(1:8,1:8) = 0;            %-- 左上
% image_encrypt_G(1:8,249:256) = 0;        %-- 右上
% image_encrypt_G(249:256,1:8) = 0;        %-- 左下
% image_encrypt_G(249:256,249:256) = 0;    %-- 右下
% 
% image_encrypt_B(1:8,1:8) = 0;            %-- 左上
% image_encrypt_B(1:8,249:256) = 0;        %-- 右上
% image_encrypt_B(249:256,1:8) = 0;        %-- 左下
% image_encrypt_B(249:256,249:256) = 0;    %-- 右下

%-- 32 * 32  中心     PSNR = 19.56
% image_encrypt_R(113:144,113:144) = 0;
% image_encrypt_G(113:144,113:144) = 0;
% image_encrypt_B(113:144,113:144) = 0;

% -- 32 * 32 -- 四角
% image_encrypt_R(1:16,1:16) = 0;           %-- 左上
% image_encrypt_R(1:16,241:256) = 0;        %-- 右上
% image_encrypt_R(241:256,1:16) = 0;        %-- 左下
% image_encrypt_R(241:256,241:256) = 0;     %-- 右下
% 
% image_encrypt_G(1:16,1:16) = 0;           %-- 左上
% image_encrypt_G(1:16,241:256) = 0;        %-- 右上
% image_encrypt_G(241:256,1:16) = 0;        %-- 左下
% image_encrypt_G(241:256,241:256) = 0;     %-- 右下

% image_encrypt_B(1:16,1:16) = 0;           %-- 左上
% image_encrypt_B(1:16,241:256) = 0;        %-- 右上
% image_encrypt_B(241:256,1:16) = 0;        %-- 左下
% image_encrypt_B(241:256,241:256) = 0;     %-- 右下

%-- 48 * 48  中心
% image_encrypt_R(105:152,105:152) = 0;
% image_encrypt_G(105:152,105:152) = 0;
% image_encrypt_B(105:152,105:152) = 0;


%% 提取图像

% VM_image_R = VM_image(:,:,1);
% VM_image_G = VM_image(:,:,2);
% VM_image_B = VM_image(:,:,3);
% 
% [VM_R_A,VM_R_H,VM_R_V,VM_R_D] = IWT(VM_image_R);
% if exist('block_pollute_','var')
%     VM_R_A = block_pollute(VM_R_A,0);
%     VM_R_H = block_pollute(VM_R_H,0);
%     VM_R_V = block_pollute(VM_R_V,0);
%     VM_R_D = block_pollute(VM_R_D,0);
% end
%     
% E1_R = zeros(N/2,N/2);
% E2_R = zeros(N/2,N/2);
% 
% [VM_G_A,VM_G_H,VM_G_V,VM_G_D] = IWT(VM_image_G);
% if exist('block_pollute_','var')
%     VM_G_A = block_pollute(VM_G_A,0);
%     VM_G_H = block_pollute(VM_G_H,0);
%     VM_G_V = block_pollute(VM_G_V,0);
%     VM_G_D = block_pollute(VM_G_D,0);
% end
% 
% E1_G = zeros(N/2,N/2);
% E2_G = zeros(N/2,N/2);
% 
% [VM_B_A,VM_B_H,VM_B_V,VM_B_D] = IWT(VM_image_B);
% if exist('block_pollute_','var')
%     VM_B_A = block_pollute(VM_B_A,0);
%     VM_B_H = block_pollute(VM_B_H,0);
%     VM_B_V = block_pollute(VM_B_V,0);
%     VM_B_D = block_pollute(VM_B_D,0);
% end
% 
% E1_B = zeros(N/2,N/2);
% E2_B = zeros(N/2,N/2);
% 
% for i=1:N/2
%     for j=1:N/2
%         if Tag_R(i,j) == 1
%             E1_R(i,j) = VM_R_H(i,j);
%             E2_R(i,j) = VM_R_V(i,j);
%         else
%             E1_R(i,j) = VM_R_V(i,j);
%             E2_R(i,j) = VM_R_H(i,j);
%         end
%         if Tag_G(i,j) == 1
%             E1_G(i,j) = VM_G_H(i,j);
%             E2_G(i,j) = VM_G_V(i,j);
%         else
%             E1_G(i,j) = VM_G_V(i,j);
%             E2_G(i,j) = VM_G_H(i,j);
%         end
%         if Tag_B(i,j) == 1
%             E1_B(i,j) = VM_B_H(i,j);
%             E2_B(i,j) = VM_B_V(i,j);
%         else
%             E1_B(i,j) = VM_B_V(i,j);
%             E2_B(i,j) = VM_B_H(i,j);
%         end
%     end
% end
% 
% image_extract_R = zeros(M,M);
% image_extract_G = zeros(M,M);
% image_extract_B = zeros(M,M);
% for i=1:M
%     for j=1:M
%         image_extract_R(i,j) = E2_R(i,j) * 10 + E1_R(i,j);
%         image_extract_G(i,j) = E2_G(i,j) * 10 + E1_G(i,j);
%         image_extract_B(i,j) = E2_B(i,j) * 10 + E1_B(i,j);
%     end
% end

% image_extract_R(image_extract_R < 0) = 0;
% image_extract_G(image_extract_G < 0) = 0;
% image_extract_B(image_extract_B < 0) = 0; 


figure;imshow(uint8(image_encrypt));title('加密后的图像');
figure;imhist(uint8(image_encrypt));title('加密后的图像的直方图');
imwrite(uint8(image_encrypt),'image_encrypt.bmp');


% %-- 保存参数
save('secret_key','K','x01','r01','x02','r02','x03','r03','x04','r04','x05','r05','x06','r06','N0','d','CR','image_compress_R_min','image_compress_R_max','image_compress_G_min','image_compress_G_max','image_compress_B_min','image_compress_B_max','t1','t2','t3','t4','t5','t6','entropy_R','entropy_G','entropy_B','block_fill');




