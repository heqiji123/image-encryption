%% 图像质量分析

clear;
clc;

% -- 1、密钥空间分析
% -- 2、密钥敏感性分析
% -- 3、像素相关性
% -- 4、信息熵
% -- 5、噪声攻击分析
% -- 6、阻塞攻击分析
% -- 7、已知明文和选择明文攻击 (做文字分析)
% -- 8、时间分析

plain_image = double(imread('plain_image.bmp'));      % -- 明文图像
image_encrypt = double(imread('image_encrypt.bmp'));  % -- 加密图像
carrier_image = double(imread('carrier_image.bmp'));  % -- 载体图像
ciper_image = double(imread('ciper_image.bmp'));      % -- VM图像


%% 密钥空间大小

%-- 1、 k  2^48  sha-384值  用于生成四翼混沌系统，LTS，LSS, TSS的初始值 x01 r01 x02 r02 x03 r03
%-- 2、 entropy_R entropy_G entropy_B  3*10^14  明文信息熵用来更新LSS混沌系统初始值x02 r02得到 { x03 r03 }  生成用与DNA解密的混沌序列
%-- 3、 t1 t2 t3 t4  4 * 10^14  用于生成四翼混沌序列的初始值 { x0 y0 z0 w0 }
%-- 4、 x01_ r01_ x02_ r02_ x03_ r03_   6 * 10^14    用于生成LTS LSS TSS的初始值 { x01 r01 x02 r02 x03 r03 }
%-- 5、 量化参数 6 * 10^14
%--     {image_compress_R_min、image_compress_R_max、image_compress_G_min、image_compress_G_max、image_compress_B_min、image_compress_B_max}

%-- 密钥空间大小  (10^14)^9 >> 2^100 


%% 图像的直方图
figure;imhist(uint8(plain_image(:,:,1)));title('明文图像的R通道的直方图');
figure;imhist(uint8(plain_image(:,:,2)));title('明文图像的G通道的直方图');
figure;imhist(uint8(plain_image(:,:,3)));title('明文图像的B通道的直方图');

figure;imhist(uint8(image_encrypt(:,:,1)));title('加密图像的R通道的直方图');
figure;imhist(uint8(image_encrypt(:,:,2)));title('加密图像的G通道的直方图');
figure;imhist(uint8(image_encrypt(:,:,3)));title('加密图像的B通道的直方图');



%% 信息熵分析（plain_image image_encrypt）
plain_image_entropy_R = entropy(plain_image(:,:,1));
plain_image_entropy_G = entropy(plain_image(:,:,2));
plain_image_entropy_B = entropy(plain_image(:,:,3));

image_encrypt_entropy_R = entropy(image_encrypt(:,:,1));
image_encrypt_entropy_G = entropy(image_encrypt(:,:,2));
image_encrypt_entropy_B = entropy(image_encrypt(:,:,3));



%% 原始图片图片像素相关性分析

[M,N,k] = size(plain_image);
plain_image_R = plain_image(:,:,1);
plain_image_G = plain_image(:,:,2);
plain_image_B = plain_image(:,:,3);

%{
先随机在0~M-1行和0~N-1列选中5000个像素点，
计算水平相关性时，选择每个点的相邻的右边的点；
计算垂直相关性时，选择每个点的相邻的下方的点；
计算对角线相关性时，选择每个点的相邻的右下方的点。
%}
Num_rand=5000;    %随机取5000对像素点
x1=ceil(rand(1,Num_rand)*(M-1));      %生成5000个1~M-1的随机整数作为行  ceil()为向上取整函数 rand生成1内的小数
y1=ceil(rand(1,Num_rand)*(N-1));      %生成5000个1~N-1的随机整数作为列
k1=ceil(rand(1,Num_rand)*(k));      %生成5000个1~k-1的随机整数作为分量
%预分配内存
XX_SP1=zeros(1,Num_rand);YY_SP1=zeros(1,Num_rand);            %水平 R G B三个分量上取1个像素点和这个像素点水平相邻的像素点
XX_CZ1=zeros(1,Num_rand);YY_CZ1=zeros(1,Num_rand);        %垂直 R G B三个分量上取1个像素点和这个像素点垂直相邻的像素点
XX_DJX1=zeros(1,Num_rand);YY_DJX1=zeros(1,Num_rand);      %对角线 R G B三个分量上取1个像素点和这个像素点对角线相邻的像素点

for i=1:Num_rand
    %水平
    XX_SP1(i)=plain_image(x1(i),y1(i),k1(i));
    YY_SP1(i)=plain_image(x1(i),y1(i)+1,k1(i)); %水平方向上相邻元素的为某个元素的右边
    
    %垂直
    XX_CZ1(i)=plain_image(x1(i),y1(i),k1(i));
    YY_CZ1(i)=plain_image(x1(i)+1,y1(i),k1(i)); %垂直方向上相邻元素的为某个元素的下边

    %对角线
    XX_DJX1(i)=plain_image(x1(i),y1(i),k1(i));
    YY_DJX1(i)=plain_image(x1(i)+1,y1(i)+1,k1(i)); %对角线方向上相邻元素的为某个元素的右下边

end
%水平   scatter用来绘制散点图
figure;scatter(XX_SP1,YY_SP1,18,'filled');xlabel(' ');ylabel(' ');title('原始图像水平相邻元素相关性点图');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);

%垂直
figure;scatter(XX_CZ1,YY_CZ1,18,'filled');xlabel(' ');ylabel(' ');title('原始图像垂直相邻元素相关性点图');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);

%对角线
figure;scatter(XX_DJX1,YY_DJX1,18,'filled');xlabel(' ');ylabel(' ');title('原始图像对角线相邻元素相关性点图');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);

%R通道
EX1_R = 0;EY1_SP_R = 0; DX1_R=0; DY1_SP_R = 0;COVXY1_SP_R = 0;    %计算水平相关性时需要的变量
EY1_CZ_R = 0; DY1_CZ_R = 0;COVXY1_CZ_R = 0;                %垂直
EY1_DJX_R = 0;DY1_DJX_R =0;COVXY1_DJX_R=0;             %对角线
%G通道
EX1_G = 0;EY1_SP_G = 0; DX1_G=0; DY1_SP_G = 0;COVXY1_SP_G = 0;    %计算水平相关性时需要的变量
EY1_CZ_G = 0; DY1_CZ_G = 0;COVXY1_CZ_G = 0;                %垂直
EY1_DJX_G = 0;DY1_DJX_G =0;COVXY1_DJX_G=0;             %对角线
%R通道
EX1_B = 0;EY1_SP_B = 0; DX1_B=0; DY1_SP_B = 0;COVXY1_SP_B = 0;    %计算水平相关性时需要的变量
EY1_CZ_B = 0; DY1_CZ_B = 0;COVXY1_CZ_B = 0;                %垂直
EY1_DJX_B = 0;DY1_DJX_B =0;COVXY1_DJX_B=0;             %对角线
I1_c = double(plain_image_R); %将像素用浮点数来表示
I2_c = double(plain_image_G); 
I3_c = double(plain_image_B); 

for i=1:Num_rand
    %统计某个随机像素点的和，用EX1来表示
    EX1_R = EX1_R + I1_c(x1(i),y1(i)); 
    EX1_G = EX1_G + I2_c(x1(i),y1(i)); 
    EX1_B = EX1_B + I3_c(x1(i),y1(i)); 
    
    %统计相邻像素的和，包括水平、垂直、对角线上的相邻像素、分别对应EY1_SP、EY1_CZ、EY1_DJX
    %R通道
    EY1_SP_R=EY1_SP_R+I1_c(x1(i)+1,y1(i));
    EY1_CZ_R=EY1_CZ_R+I1_c(x1(i),y1(i)+1);
    EY1_DJX_R=EY1_DJX_R+I1_c(x1(i)+1,y1(i)+1);
    %G通道
    EY1_SP_G=EY1_SP_G+I2_c(x1(i)+1,y1(i));
    EY1_CZ_G=EY1_CZ_G+I2_c(x1(i),y1(i)+1);
    EY1_DJX_G=EY1_DJX_G+I2_c(x1(i)+1,y1(i)+1);
    %B通道
    EY1_SP_B=EY1_SP_B+I3_c(x1(i)+1,y1(i));
    EY1_CZ_B=EY1_CZ_B+I3_c(x1(i),y1(i)+1);
    EY1_DJX_B=EY1_DJX_B+I3_c(x1(i)+1,y1(i)+1);
end
%求均值
%统一在循环外除以像素点对数5000，可减少运算次数
% R通道
EX1_R=EX1_R/Num_rand;
EY1_SP_R=EY1_SP_R/Num_rand;
EY1_CZ_R=EY1_CZ_R/Num_rand;
EY1_DJX_R=EY1_DJX_R/Num_rand;
% G通道
EX1_G=EX1_G/Num_rand;
EY1_SP_G=EY1_SP_G/Num_rand;
EY1_CZ_G=EY1_CZ_G/Num_rand;
EY1_DJX_G=EY1_DJX_G/Num_rand;
% B通道
EX1_B=EX1_B/Num_rand;
EY1_SP_B=EY1_SP_B/Num_rand;
EY1_CZ_B=EY1_CZ_B/Num_rand;
EY1_DJX_B=EY1_DJX_B/Num_rand;
for i=1:Num_rand  
    %求方差
    %第一个像素点的D，水平、垂直、对角线时计算得出第一个像素点的D相同，统一用DX表示
    DX1_R=DX1_R+(I1_c(x1(i),y1(i))-EX1_R)^2;
    DX1_G=DX1_G+(I2_c(x1(i),y1(i))-EX1_G)^2;
    DX1_B=DX1_B+(I3_c(x1(i),y1(i))-EX1_B)^2;
    %第二个像素点的D，水平、垂直、对角线的D分别对应DY1_SP、DY1_CZ、DY1_DJX
    %R通道
    DY1_SP_R=DY1_SP_R+(I1_c(x1(i)+1,y1(i))-EY1_SP_R)^2;
    DY1_CZ_R=DY1_CZ_R+(I1_c(x1(i),y1(i)+1)-EY1_CZ_R)^2;
    DY1_DJX_R=DY1_DJX_R+(I1_c(x1(i)+1,y1(i)+1)-EY1_DJX_R)^2;
    %G通道
    DY1_SP_G=DY1_SP_G+(I2_c(x1(i)+1,y1(i))-EY1_SP_G)^2;
    DY1_CZ_G=DY1_CZ_G+(I2_c(x1(i),y1(i)+1)-EY1_CZ_G)^2;
    DY1_DJX_G=DY1_DJX_G+(I2_c(x1(i)+1,y1(i)+1)-EY1_DJX_G)^2;
    %B通道
    DY1_SP_B=DY1_SP_B+(I3_c(x1(i)+1,y1(i))-EY1_SP_B)^2;
    DY1_CZ_B=DY1_CZ_B+(I3_c(x1(i),y1(i)+1)-EY1_CZ_B)^2;
    DY1_DJX_B=DY1_DJX_B+(I3_c(x1(i)+1,y1(i)+1)-EY1_DJX_B)^2;
    %两个相邻像素点相关函数的计算，水平、垂直、对角线
    
    %求协方差
    %R通道
    COVXY1_SP_R=COVXY1_SP_R+(I1_c(x1(i),y1(i))-EX1_R)*(I1_c(x1(i)+1,y1(i))-EY1_SP_R);
    COVXY1_CZ_R=COVXY1_CZ_R+(I1_c(x1(i),y1(i))-EX1_R)*(I1_c(x1(i),y1(i)+1)-EY1_CZ_R);
    COVXY1_DJX_R=COVXY1_DJX_R+(I1_c(x1(i),y1(i))-EX1_R)*(I1_c(x1(i)+1,y1(i)+1)-EY1_DJX_R);
    %G通道
    COVXY1_SP_G=COVXY1_SP_G+(I2_c(x1(i),y1(i))-EX1_G)*(I2_c(x1(i)+1,y1(i))-EY1_SP_G);
    COVXY1_CZ_G=COVXY1_CZ_G+(I2_c(x1(i),y1(i))-EX1_G)*(I2_c(x1(i),y1(i)+1)-EY1_CZ_G);
    COVXY1_DJX_G=COVXY1_DJX_G+(I2_c(x1(i),y1(i))-EX1_G)*(I2_c(x1(i)+1,y1(i)+1)-EY1_DJX_G);
    %B通道
    COVXY1_SP_B=COVXY1_SP_B+(I3_c(x1(i),y1(i))-EX1_B)*(I3_c(x1(i)+1,y1(i))-EY1_SP_B);
    COVXY1_CZ_B=COVXY1_CZ_B+(I3_c(x1(i),y1(i))-EX1_B)*(I3_c(x1(i),y1(i)+1)-EY1_CZ_B);
    COVXY1_DJX_B=COVXY1_DJX_B+(I3_c(x1(i),y1(i))-EX1_B)*(I3_c(x1(i)+1,y1(i)+1)-EY1_DJX_B);
end
%统一在循环外除以像素点对数5000，可减少运算次数
%求方差和协方差
%R通道
DX1_R=DX1_R/Num_rand;
DY1_SP_R=DY1_SP_R/Num_rand;
DY1_CZ_R=DY1_CZ_R/Num_rand;
DY1_DJX_R=DY1_DJX_R/Num_rand;
COVXY1_SP_R=COVXY1_SP_R/Num_rand;
COVXY1_CZ_R=COVXY1_CZ_R/Num_rand;
COVXY1_DJX_R=COVXY1_DJX_R/Num_rand;
%G通道
DX1_G=DX1_G/Num_rand;
DY1_SP_G=DY1_SP_G/Num_rand;
DY1_CZ_G=DY1_CZ_G/Num_rand;
DY1_DJX_G=DY1_DJX_G/Num_rand;
COVXY1_SP_G=COVXY1_SP_G/Num_rand;
COVXY1_CZ_G=COVXY1_CZ_G/Num_rand;
COVXY1_DJX_G=COVXY1_DJX_G/Num_rand;
%B通道
DX1_B=DX1_B/Num_rand;
DY1_SP_B=DY1_SP_B/Num_rand;
DY1_CZ_B=DY1_CZ_B/Num_rand;
DY1_DJX_B=DY1_DJX_B/Num_rand;
COVXY1_SP_B=COVXY1_SP_B/Num_rand;
COVXY1_CZ_B=COVXY1_CZ_B/Num_rand;
COVXY1_DJX_B=COVXY1_DJX_B/Num_rand;

%求相关性系数
%水平、垂直、对角线的相关性
%R通道
RXY1_SP_R=COVXY1_SP_R/sqrt(DX1_R*DY1_SP_R);
RXY1_CZ_R=COVXY1_CZ_R/sqrt(DX1_R*DY1_CZ_R);
RXY1_DJX_R=COVXY1_DJX_R/sqrt(DX1_R*DY1_DJX_R);
%G通道
RXY1_SP_G=COVXY1_SP_G/sqrt(DX1_G*DY1_SP_G);
RXY1_CZ_G=COVXY1_CZ_G/sqrt(DX1_G*DY1_CZ_G);
RXY1_DJX_G=COVXY1_DJX_G/sqrt(DX1_G*DY1_DJX_G);
%B通道
RXY1_SP_B=COVXY1_SP_B/sqrt(DX1_B*DY1_SP_B);
RXY1_CZ_B=COVXY1_CZ_B/sqrt(DX1_B*DY1_CZ_B);
RXY1_DJX_B=COVXY1_DJX_B/sqrt(DX1_B*DY1_DJX_B);


disp('加密前R、G、B通道相关性：');
disp(['加密前图片R通道相关性：','  水平相关性=',num2str(RXY1_SP_R),'  垂直相关性=',num2str(RXY1_CZ_R),'  对角线相关性=',num2str(RXY1_DJX_R)]);
disp(['加密前图片G通道相关性：','  水平相关性=',num2str(RXY1_SP_G),'  垂直相关性=',num2str(RXY1_CZ_G),'  对角线相关性=',num2str(RXY1_DJX_G)]);
disp(['加密前图片B通道相关性：','  水平相关性=',num2str(RXY1_SP_B),'  垂直相关性=',num2str(RXY1_CZ_B),'  对角线相关性=',num2str(RXY1_DJX_B)]);
disp('');


%% 加密图片的相关像素相关性
[M,N,k] = size(image_encrypt);
image_encrypt_R = image_encrypt(:,:,1);
image_encrypt_G = image_encrypt(:,:,2);
image_encrypt_B = image_encrypt(:,:,3);

%{
先随机在0~M-1行和0~N-1列选中5000个像素点，
计算水平相关性时，选择每个点的相邻的右边的点；
计算垂直相关性时，选择每个点的相邻的下方的点；
计算对角线相关性时，选择每个点的相邻的右下方的点。
%}
Num_rand=5000;    %随机取5000对像素点
x2=ceil(rand(1,Num_rand)*(M-1));      %生成5000个1~M-1的随机整数作为行  ceil()为向上取整函数 rand生成1内的小数
y2=ceil(rand(1,Num_rand)*(N-1));      %生成5000个1~N-1的随机整数作为列
k2=ceil(rand(1,Num_rand)*k); 
%预分配内存
XX_SP2=zeros(1,Num_rand);YY_SP2=zeros(1,Num_rand);        %水平 R G B三个分量上取1个像素点和这个像素点水平相邻的像素点

XX_CZ2=zeros(1,Num_rand);YY_CZ2=zeros(1,Num_rand);        %垂直 R G B三个分量上取1个像素点和这个像素点垂直相邻的像素点

XX_DJX2=zeros(1,Num_rand);YY_DJX2=zeros(1,Num_rand);      %对角线 R G B三个分量上取1个像素点和这个像素点对角线相邻的像素点

for i=1:Num_rand
    %水平
    XX_SP2(i)=image_encrypt(x2(i),y2(i),k2(i));
    YY_SP2(i)=image_encrypt(x2(i),y2(i)+1,k2(i));          %水平方向上相邻元素的为某个元素的右边
    %垂直
    XX_CZ2(i)=image_encrypt(x2(i),y2(i),k2(i));
    YY_CZ2(i)=image_encrypt(x2(i)+1,y2(i),k2(i));          %垂直方向上相邻元素的为某个元素的下边
    %对角线
    XX_DJX2(i)=image_encrypt(x2(i),y2(i),k2(i));
    YY_DJX2(i)=image_encrypt(x2(i)+1,y2(i)+1,k2(i));       %对角线方向上相邻元素的为某个元素的右下边
end
%水平   scatter用来绘制散点图
figure;scatter(XX_SP2,YY_SP2,18,'filled');xlabel(' ');ylabel(' ');title('');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);

%垂直
figure;scatter(XX_CZ2,YY_CZ2,18,'filled');xlabel(' ');ylabel(' ');title('');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);

%对角线
figure;scatter(XX_DJX2,YY_DJX2,18,'filled');xlabel(' ');ylabel(' ');title('');axis([0 255,0 255]);set(gca,'XTick',0:15:255);set(gca,'YTick',0:15:255);

%R通道
EX2_R=0;
EY2_SP_R=0;DX2_R=0;DY2_SP_R=0;COVXY2_SP_R=0;        %计算水平相关性时需要的变量
EY2_CZ_R=0;DY2_CZ_R=0;COVXY2_CZ_R=0;                %垂直
EY2_DJX_R=0;DY2_DJX_R=0;COVXY2_DJX_R=0;             %对角线
%G通道
EX2_G=0;
EY2_SP_G=0;DX2_G=0;DY2_SP_G=0;COVXY2_SP_G=0;
EY2_CZ_G=0;DY2_CZ_G=0;COVXY2_CZ_G=0;
EY2_DJX_G=0;DY2_DJX_G=0;COVXY2_DJX_G=0;
%B通道
EX2_B=0;
EY2_SP_B=0;DX2_B=0;DY2_SP_B=0;COVXY2_SP_B=0;
EY2_CZ_B=0;DY2_CZ_B=0;COVXY2_CZ_B=0;
EY2_DJX_B=0;DY2_DJX_B=0;COVXY2_DJX_B=0;

R_c=double(image_encrypt_R); %将像素用浮点数来表示
G_c=double(image_encrypt_G);
B_c=double(image_encrypt_B);

for i=1:Num_rand
    %统计某个随机像素点在加密通道R、G、B上的和，用EX2_R、EX2_G、EX2_B来表示
    EX2_R=EX2_R+R_c(x2(i),y2(i)); 
    EX2_G=EX2_G+G_c(x2(i),y2(i)); 
    EX2_B=EX2_B+B_c(x2(i),y2(i)); 
    %统计相邻像素的和，包括水平、垂直、对角线上的相邻像素、分别对应EY2_SP、EY2_CZ、EY2_DJX
    %R通道
    EY2_SP_R=EY2_SP_R+R_c(x2(i)+1,y2(i));
    EY2_CZ_R=EY2_CZ_R+R_c(x2(i),y2(i)+1);
    EY2_DJX_R=EY2_DJX_R+R_c(x2(i)+1,y2(i)+1);
    %G通道
    EY2_SP_G=EY2_SP_G+G_c(x2(i)+1,y2(i));
    EY2_CZ_G=EY2_CZ_G+G_c(x2(i),y2(i)+1);
    EY2_DJX_G=EY2_DJX_G+G_c(x2(i)+1,y2(i)+1);
    %B通道
    EY2_SP_B=EY2_SP_B+B_c(x2(i)+1,y2(i));
    EY2_CZ_B=EY2_CZ_B+B_c(x2(i),y2(i)+1);
    EY2_DJX_B=EY2_DJX_B+B_c(x2(i)+1,y2(i)+1);
end
%求均值
%统一在循环外除以像素点对数5000，可减少运算次数
% R通道
EX2_R=EX2_R/Num_rand;
EY2_SP_R=EY2_SP_R/Num_rand;
EY2_CZ_R=EY2_CZ_R/Num_rand;
EY2_DJX_R=EY2_DJX_R/Num_rand;
% G通道
EX2_G=EX2_G/Num_rand;
EY2_SP_G=EY2_SP_G/Num_rand;
EY2_CZ_G=EY2_CZ_G/Num_rand;
EY2_DJX_G=EY2_DJX_G/Num_rand;
% B通道
EX2_B=EX2_B/Num_rand;
EY2_SP_B=EY2_SP_B/Num_rand;
EY2_CZ_B=EY2_CZ_B/Num_rand;
EY2_DJX_B=EY2_DJX_B/Num_rand;
for i=1:Num_rand  
    %求方差
    %第一个像素点的D，水平、垂直、对角线时计算得出第一个像素点的D相同，统一用DX表示
    DX2_R=DX2_R+(R_c(x2(i),y2(i))-EX2_R)^2;
    DX2_G=DX2_G+(G_c(x2(i),y2(i))-EX2_G)^2;
    DX2_B=DX2_B+(B_c(x2(i),y2(i))-EX2_B)^2;
    %第二个像素点的D，水平、垂直、对角线的D分别对应DY1_SP、DY1_CZ、DY1_DJX
    %R通道
    DY2_SP_R=DY2_SP_R+(R_c(x2(i)+1,y2(i))-EY2_SP_R)^2;
    DY2_CZ_R=DY2_CZ_R+(R_c(x2(i),y2(i)+1)-EY2_CZ_R)^2;
    DY2_DJX_R=DY2_DJX_R+(R_c(x2(i)+1,y2(i)+1)-EY2_DJX_R)^2;
    %G通道
    DY2_SP_G=DY2_SP_G+(G_c(x2(i)+1,y2(i))-EY2_SP_G)^2;
    DY2_CZ_G=DY2_CZ_G+(G_c(x2(i),y2(i)+1)-EY2_CZ_G)^2;
    DY2_DJX_G=DY2_DJX_G+(G_c(x2(i)+1,y2(i)+1)-EY2_DJX_G)^2;
    %B通道
    DY2_SP_B=DY2_SP_B+(B_c(x2(i)+1,y2(i))-EY2_SP_B)^2;
    DY2_CZ_B=DY2_CZ_B+(B_c(x2(i),y2(i)+1)-EY2_CZ_B)^2;
    DY2_DJX_B=DY2_DJX_B+(B_c(x2(i)+1,y2(i)+1)-EY2_DJX_B)^2;
    %两个相邻像素点相关函数的计算，水平、垂直、对角线
    
    %求协方差
    %R通道
    COVXY2_SP_R=COVXY2_SP_R+(R_c(x2(i),y2(i))-EX2_R)*(R_c(x2(i)+1,y2(i))-EY2_SP_R);
    COVXY2_CZ_R=COVXY2_CZ_R+(R_c(x2(i),y2(i))-EX2_R)*(R_c(x2(i),y2(i)+1)-EY2_CZ_R);
    COVXY2_DJX_R=COVXY2_DJX_R+(R_c(x2(i),y2(i))-EX2_R)*(R_c(x2(i)+1,y2(i)+1)-EY2_DJX_R);
    %G通道
    COVXY2_SP_G=COVXY2_SP_G+(G_c(x2(i),y2(i))-EX2_G)*(G_c(x2(i)+1,y2(i))-EY2_SP_G);
    COVXY2_CZ_G=COVXY2_CZ_G+(G_c(x2(i),y2(i))-EX2_G)*(G_c(x2(i),y2(i)+1)-EY2_CZ_G);
    COVXY2_DJX_G=COVXY2_DJX_G+(G_c(x2(i),y2(i))-EX2_G)*(G_c(x2(i)+1,y2(i)+1)-EY2_DJX_G);
    %B通道
    COVXY2_SP_B=COVXY2_SP_B+(B_c(x2(i),y2(i))-EX2_B)*(B_c(x2(i)+1,y2(i))-EY2_SP_B);
    COVXY2_CZ_B=COVXY2_CZ_B+(B_c(x2(i),y2(i))-EX2_B)*(B_c(x2(i),y2(i)+1)-EY2_CZ_B);
    COVXY2_DJX_B=COVXY2_DJX_B+(B_c(x2(i),y2(i))-EX2_B)*(B_c(x2(i)+1,y2(i)+1)-EY2_DJX_B);
end
%统一在循环外除以像素点对数5000，可减少运算次数
%求方差和协方差
%R通道
DX2_R=DX2_R/Num_rand;
DY2_SP_R=DY2_SP_R/Num_rand;
DY2_CZ_R=DY2_CZ_R/Num_rand;
DY2_DJX_R=DY2_DJX_R/Num_rand;
COVXY2_SP_R=COVXY2_SP_R/Num_rand;
COVXY2_CZ_R=COVXY2_CZ_R/Num_rand;
COVXY2_DJX_R=COVXY2_DJX_R/Num_rand;
%G通道
DX2_G=DX2_G/Num_rand;
DY2_SP_G=DY2_SP_G/Num_rand;
DY2_CZ_G=DY2_CZ_G/Num_rand;
DY2_DJX_G=DY2_DJX_G/Num_rand;
COVXY2_SP_G=COVXY2_SP_G/Num_rand;
COVXY2_CZ_G=COVXY2_CZ_G/Num_rand;
COVXY2_DJX_G=COVXY2_DJX_G/Num_rand;
%B通道
DX2_B=DX2_B/Num_rand;
DY2_SP_B=DY2_SP_B/Num_rand;
DY2_CZ_B=DY2_CZ_B/Num_rand;
DY2_DJX_B=DY2_DJX_B/Num_rand;
COVXY2_SP_B=COVXY2_SP_B/Num_rand;
COVXY2_CZ_B=COVXY2_CZ_B/Num_rand;
COVXY2_DJX_B=COVXY2_DJX_B/Num_rand;

%求相关性系数
%水平、垂直、对角线的相关性
%R通道
RXY2_SP_R=COVXY2_SP_R/sqrt(DX2_R*DY2_SP_R);
RXY2_CZ_R=COVXY2_CZ_R/sqrt(DX2_R*DY2_CZ_R);
RXY2_DJX_R=COVXY2_DJX_R/sqrt(DX2_R*DY2_DJX_R);
%G通道
RXY2_SP_G=COVXY2_SP_G/sqrt(DX2_G*DY2_SP_G);
RXY2_CZ_G=COVXY2_CZ_G/sqrt(DX2_G*DY2_CZ_G);
RXY2_DJX_G=COVXY2_DJX_G/sqrt(DX2_G*DY2_DJX_G);
%B通道
RXY2_SP_B=COVXY2_SP_B/sqrt(DX2_B*DY2_SP_B);
RXY2_CZ_B=COVXY2_CZ_B/sqrt(DX2_B*DY2_CZ_B);
RXY2_DJX_B=COVXY2_DJX_B/sqrt(DX2_B*DY2_DJX_B);

disp('加密后R、G、B通道相关性：');
disp(['加密后图片R通道相关性：','  水平相关性=',num2str(RXY2_SP_R),'  垂直相关性=',num2str(RXY2_CZ_R),'  对角线相关性=',num2str(RXY2_DJX_R)]);
disp(['加密后图片G通道相关性：','  水平相关性=',num2str(RXY2_SP_G),'  垂直相关性=',num2str(RXY2_CZ_G),'  对角线相关性=',num2str(RXY2_DJX_G)]);
disp(['加密后图片B通道相关性：','  水平相关性=',num2str(RXY2_SP_B),'  垂直相关性=',num2str(RXY2_CZ_B),'  对角线相关性=',num2str(RXY2_DJX_B)]);
disp('');






