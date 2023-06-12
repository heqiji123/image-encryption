clc;
clear;

load('secret_key','K','x01','r01','x02','r02','x03','r03','x04','r04','x05','r05','x06','r06','N0','d','CR','image_compress_R_min','image_compress_R_max','image_compress_G_min','image_compress_G_max','image_compress_B_min','image_compress_B_max','t1','t2','t3','t4','t5','t6','entropy_R','entropy_G','entropy_B','block_fill');
plain_image = imread('plain_image.bmp','bmp');
[~,N,~] = size(plain_image);
image_encrypt = double(imread('image_encrypt.bmp','bmp'));
image_encrypt_R = image_encrypt(:,:,1);
image_encrypt_G = image_encrypt(:,:,2);
image_encrypt_B = image_encrypt(:,:,3);


%-- 椒盐噪声、乘性噪声、高斯噪声
% image_encrypt_R = double(imnoise(uint8(image_encrypt_R),'salt & pepper',0.001));
% image_encrypt_G = double(imnoise(uint8(image_encrypt_G),'salt & pepper',0.001));
% image_encrypt_B = double(imnoise(uint8(image_encrypt_B),'salt & pepper',0.001));

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
% 
% image_encrypt_B(1:16,1:16) = 0;           %-- 左上
% image_encrypt_B(1:16,241:256) = 0;        %-- 右上
% image_encrypt_B(241:256,1:16) = 0;        %-- 左下
% image_encrypt_B(241:256,241:256) = 0;     %-- 右下

%-- 48 * 48  中心
% image_encrypt_R(105:152,105:152) = 0;
% image_encrypt_G(105:152,105:152) = 0;
% image_encrypt_B(105:152,105:152) = 0;

% -- 48 * 48 -- 四角
image_encrypt_R(1:24,1:24) = 0;           %-- 左上
image_encrypt_R(1:24,233:256) = 0;        %-- 右上
image_encrypt_R(233:256,1:24) = 0;        %-- 左下
image_encrypt_R(233:256,233:256) = 0;     %-- 右下

image_encrypt_G(1:24,1:24) = 0;           %-- 左上
image_encrypt_G(1:24,233:256) = 0;        %-- 右上
image_encrypt_G(233:256,1:24) = 0;        %-- 左下
image_encrypt_G(233:256,233:256) = 0;     %-- 右下

image_encrypt_B(1:24,1:24) = 0;           %-- 左上
image_encrypt_B(1:24,233:256) = 0;        %-- 右上
image_encrypt_B(233:256,1:24) = 0;        %-- 左下
image_encrypt_B(233:256,233:256) = 0;     %-- 右下

% K(48,8) = '1';
k = zeros(1,48);
for i=1:48
    k(i) = bin2dec(K(i,:));
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



%% Josephus 的参数
x06 = (k(37)+k(38)+k(39)+k(40)+k(41)+k(42)) / 6;
r06 = bitxor(k(43)+k(44)+k(45),k(46)+k(47)+k(48));



%% 生成测量矩阵

CR = 0.25;                                                %   压缩比率
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


% -- LSS生成用于DNA_encode DNA_decode的混沌序列
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



%% 解密图像
% image_decrypt_R_is = zeros(M,M);
% image_decrypt_G_is = zeros(M,M);
% image_decrypt_B_is = zeros(M,M);

% -- DNA编码
% R_DNA_1 = DNA_encode_(image_encrypt_R, Lss_DNA_decode);
% G_DNA_1 = DNA_encode_(image_encrypt_G, Lss_DNA_decode);
% B_DNA_1 = DNA_encode_(image_encrypt_B, Lss_DNA_decode);

% -- R、G、B、分量间的加，减，异或
% B_DNA_C1 = DNA_operation(G_DNA_1,B_DNA_1,3);   %异或
% G_DNA_C1 = DNA_operation(R_DNA_1,G_DNA_1,2);   %减法
% R_DNA_C1 = DNA_operation(B_DNA_C1,R_DNA_1,1);  %加法a

% -- DNA逆运算
% R_DNA_C2 = DNA_compute_inverse(R_DNA_1, T);
% G_DNA_C2 = DNA_compute_inverse(G_DNA_1, T);
% B_DNA_C2 = DNA_compute_inverse(B_DNA_1, T);

% -- DNA解码
% image_decrypt_R_is = DNA_decode_(R_DNA_C2,Lss_DNA_encode);
% image_decrypt_G_is = DNA_decode_(G_DNA_C2,Lss_DNA_encode);
% image_decrypt_B_is = DNA_decode_(B_DNA_C2,Lss_DNA_encode);

for i=block_nums:-1:1
    
    %-- 上一块
    if i==1
        d_sub_block_R_last = divide_blocks(block_size,image_decrypt_R_is,block_nums);
        d_sub_block_G_last = divide_blocks(block_size,image_decrypt_G_is,block_nums);
        d_sub_block_B_last = divide_blocks(block_size,image_decrypt_B_is,block_nums);

        d_sub_block_R_DNA_last = DNA_encode_(d_sub_block_R_last,Lss_DNA_encode(block_nums));
        d_sub_block_G_DNA_last = DNA_encode_(d_sub_block_G_last,Lss_DNA_encode(block_nums));
        d_sub_block_B_DNA_last = DNA_encode_(d_sub_block_B_last,Lss_DNA_encode(block_nums));
    else
        d_sub_block_R_last = divide_blocks(block_size,image_encrypt_R,i-1);
        d_sub_block_G_last = divide_blocks(block_size,image_encrypt_G,i-1);
        d_sub_block_B_last = divide_blocks(block_size,image_encrypt_B,i-1);

        d_sub_block_R_DNA_last = DNA_encode_(d_sub_block_R_last,Lss_DNA_decode(i-1));
        d_sub_block_G_DNA_last = DNA_encode_(d_sub_block_G_last,Lss_DNA_decode(i-1));
        d_sub_block_B_DNA_last = DNA_encode_(d_sub_block_B_last,Lss_DNA_decode(i-1));
    end

    d_sub_block_R = divide_blocks(block_size,image_encrypt_R,i);
    d_sub_block_G = divide_blocks(block_size,image_encrypt_G,i);
    d_sub_block_B = divide_blocks(block_size,image_encrypt_B,i);
    
    %-- DNA编码
    d_sub_block_R_DNA = DNA_encode_(d_sub_block_R,Lss_DNA_decode(i));
    d_sub_block_G_DNA = DNA_encode_(d_sub_block_G,Lss_DNA_decode(i));
    d_sub_block_B_DNA = DNA_encode_(d_sub_block_B,Lss_DNA_decode(i));
    
    %-- DNA逆运算
    d_sub_block_R_DNA_ic1 = DNA_compute_inverse(d_sub_block_R_DNA,T(i));
    d_sub_block_G_DNA_ic1 = DNA_compute_inverse(d_sub_block_G_DNA,T(i));
    d_sub_block_B_DNA_ic1 = DNA_compute_inverse(d_sub_block_B_DNA,T(i));
    
    %-- 与上一块DNA相加
    d_sub_block_R_DNA_ic2 = DNA_operation(d_sub_block_R_DNA_ic1,d_sub_block_R_DNA_last,1);
    d_sub_block_G_DNA_ic2 = DNA_operation(d_sub_block_G_DNA_ic1,d_sub_block_G_DNA_last,1);
    d_sub_block_B_DNA_ic2 = DNA_operation(d_sub_block_B_DNA_ic1,d_sub_block_B_DNA_last,1);
    
    %-- DNA解码
    d_sub_block_R = DNA_decode_(d_sub_block_R_DNA_ic2,Lss_DNA_encode(i));
    d_sub_block_G = DNA_decode_(d_sub_block_G_DNA_ic2,Lss_DNA_encode(i));
    d_sub_block_B = DNA_decode_(d_sub_block_B_DNA_ic2,Lss_DNA_encode(i));
    
    block_q = M / block_size;
    xx=floor(i/block_q)+1;      %第几大行
    yy=mod(i,block_q);          %第几大列
    if yy==0
        xx=xx-1;
        yy=block_q;
    end
    image_decrypt_R_is(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = d_sub_block_R;
    image_decrypt_G_is(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = d_sub_block_G;
    image_decrypt_B_is(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = d_sub_block_B;
    
end

%%%%%%%%%%%%%%%%%%%%%%%

% for i=block_nums:-1:1
%     
%     %-- 分块
%     d_sub_block_R = divide_blocks(block_size,image_encrypt_R,i);
%     d_sub_block_G = divide_blocks(block_size,image_encrypt_G,i);
%     d_sub_block_B = divide_blocks(block_size,image_encrypt_B,i);
%     
%     %-- DNA编码
%     d_sub_block_R_DNA = DNA_encode_(d_sub_block_R,Lss_DNA_decode(i));
%     d_sub_block_G_DNA = DNA_encode_(d_sub_block_G,Lss_DNA_decode(i));
%     d_sub_block_B_DNA = DNA_encode_(d_sub_block_B,Lss_DNA_decode(i));
%     
%     %-- DNA逆运算
%     d_sub_block_R_DNA_ic1 = DNA_compute_inverse(d_sub_block_R_DNA,T(i));
%     d_sub_block_G_DNA_ic1 = DNA_compute_inverse(d_sub_block_G_DNA,T(i));
%     d_sub_block_B_DNA_ic1 = DNA_compute_inverse(d_sub_block_B_DNA,T(i));
%     
%     %-- 与上一块DNA相加
%     d_sub_block_B_DNA_ic2 = DNA_operation(d_sub_block_B_DNA_ic1,d_sub_block_G_DNA_ic1,3);
%     d_sub_block_G_DNA_ic2 = DNA_operation(d_sub_block_G_DNA_ic1,d_sub_block_R_DNA_ic1,3);
%     d_sub_block_R_DNA_ic2 = DNA_operation(d_sub_block_R_DNA_ic1,d_sub_block_B_DNA_ic2,3);
%     
%     %-- DNA解码
%     d_sub_block_R = DNA_decode_(d_sub_block_R_DNA_ic2,Lss_DNA_encode(i));
%     d_sub_block_G = DNA_decode_(d_sub_block_G_DNA_ic2,Lss_DNA_encode(i));
%     d_sub_block_B = DNA_decode_(d_sub_block_B_DNA_ic2,Lss_DNA_encode(i));
%     
%     block_q = M / block_size;
%     xx=floor(i/block_q)+1;      %第几大行
%     yy=mod(i,block_q);          %第几大列
%     if yy==0
%         xx=xx-1;
%         yy=block_q;
%     end
%     image_decrypt_R_is(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = d_sub_block_R;
%     image_decrypt_G_is(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = d_sub_block_G;
%     image_decrypt_B_is(block_size*(xx-1)+1:block_size*xx,block_size*(yy-1)+1:block_size*yy) = d_sub_block_B;
%     
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%

%-- 逆置乱
SS = [image_decrypt_R_is,image_decrypt_G_is,image_decrypt_B_is];
start = mod(round(x06 * 1000),3*M);
step = mod(round(r06*1000),19) + 1;
SS_ = Josephus_permutation_inverse(SS,start,step);
image_decrypt_R = SS_(:,1:M);
image_decrypt_G = SS_(:,M+1:2*M);
image_decrypt_B = SS_(:,2*M+1:3*M);


%% 逆量化解密图像


I_compress_R_qiv = zeros(M,M);
I_compress_G_qiv = zeros(M,M);
I_compress_B_qiv = zeros(M,M);
for i=1:M
    for j=1:M
        I_compress_R_qiv(i,j) = (image_decrypt_R(i,j)*(image_compress_R_max-image_compress_R_min)/255) + image_compress_R_min;
        I_compress_G_qiv(i,j) = (image_decrypt_G(i,j)*(image_compress_G_max-image_compress_G_min)/255) + image_compress_G_min;
        I_compress_B_qiv(i,j) = (image_decrypt_B(i,j)*(image_compress_B_max-image_compress_B_min)/255) + image_compress_B_min;
    end
end


%%  重构图像

image_reconstructed_R = TDPG (I_compress_R_qiv, Phi1, Phi2, image_size, image_size, num_levels);  
image_reconstructed_G = TDPG (I_compress_G_qiv, Phi1, Phi2, image_size, image_size, num_levels);  
image_reconstructed_B = TDPG (I_compress_B_qiv, Phi1, Phi2, image_size, image_size, num_levels);  

image_reconstructed = zeros(image_size,image_size,3);
image_reconstructed(:,:,1) = image_reconstructed_R;
image_reconstructed(:,:,2) = image_reconstructed_G;
image_reconstructed(:,:,3) = image_reconstructed_B;


figure;imshow(uint8(image_reconstructed));title('重构的图像');
imwrite(uint8(image_reconstructed),'image_reconstructed.bmp');

MSE = 0;
for i=1:N
    for j=1:N
        MSE = MSE+(double(plain_image(i,j))-double(image_reconstructed(i,j)))^2;
    end
end
MSE = MSE/(N*N);
PSNR = 10*log10((255*255)/MSE);  

% test_image_NPCR = double(imread('image_encrypt_lena_with_t1.bmp','bmp'));
% test_image_1 = double(imread('image_encrypt_lena.bmp','bmp'));
% D = 0;
% for i = 1:N
%     for j = 1:N
%         if test_image_1(i,j) ~= test_image_NPCR(i,j)
%             D = D + 1;
%         end
%     end
% end
% NPCR = D / (N * N);