%% 计算密文图像相邻相关系数

image_encrypt = double(imread('image_encrypt.bmp'));  % -- 加密图像

[M,N,k] = size(image_encrypt);

Num_rand=5000;    %随机取5000对像素点
x2=ceil(rand(1,Num_rand)*(M-1));      %生成5000个1~M-1的随机整数作为行  ceil()为向上取整函数 rand生成1内的小数
y2=ceil(rand(1,Num_rand)*(N-1));      %生成5000个1~N-1的随机整数作为列
k2=ceil(rand(1,Num_rand)*(k-1)); 
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


EX2=0;
EY2_SP=0;DX2=0;DY2_SP=0;COVXY2_SP=0;        %计算水平相关性时需要的变量
EY2_CZ=0;DY2_CZ=0;COVXY2_CZ=0;                %垂直
EY2_DJX=0;DY2_DJX=0;COVXY2_DJX=0;             %对角线

for i=1:Num_rand
    %统计某个随机像素点在加密通道上的和，用EX2来表示
    EX2=EX2+image_encrypt(x2(i),y2(i),k2(i)); 
    %统计相邻像素的和，包括水平、垂直、对角线上的相邻像素、分别对应EY2_SP、EY2_CZ、EY2_DJX
    EY2_SP=EY2_SP+image_encrypt(x2(i)+1,y2(i),k2(i));
    EY2_CZ=EY2_CZ+image_encrypt(x2(i),y2(i)+1,k2(i));
    EY2_DJX=EY2_DJX+image_encrypt(x2(i)+1,y2(i)+1,k2(i));
end
%求均值
%统一在循环外除以像素点对数5000，可减少运算次数
EX2=EX2/Num_rand;
EY2_SP=EY2_SP/Num_rand;
EY2_CZ=EY2_CZ/Num_rand;
EY2_DJX=EY2_DJX/Num_rand;

for i=1:Num_rand  
    %求方差
    %第一个像素点的D，水平、垂直、对角线时计算得出第一个像素点的D相同，统一用DX表示
    DX2=DX2+(image_encrypt(x2(i),y2(i),k2(i))-EX2)^2;

    %第二个像素点的D，水平、垂直、对角线的D分别对应DY1_SP、DY1_CZ、DY1_DJX
    DY2_SP=DY2_SP+(image_encrypt(x2(i)+1,y2(i),k2(i))-EY2_SP)^2;
    DY2_CZ=DY2_CZ+(image_encrypt(x2(i),y2(i)+1,k2(i))-EY2_CZ)^2;
    DY2_DJX=DY2_DJX+(image_encrypt(x2(i)+1,y2(i)+1,k2(i))-EY2_DJX)^2;
   
    %两个相邻像素点相关函数的计算，水平、垂直、对角线
    
    %求协方差
    %R通道
    COVXY2_SP=COVXY2_SP+(image_encrypt(x2(i),y2(i),k2(i))-EX2)*(image_encrypt(x2(i)+1,y2(i),k2(i))-EY2_SP);
    COVXY2_CZ=COVXY2_CZ+(image_encrypt(x2(i),y2(i),k2(i))-EX2)*(image_encrypt(x2(i),y2(i)+1,k2(i))-EY2_CZ);
    COVXY2_DJX=COVXY2_DJX+(image_encrypt(x2(i),y2(i),k2(i))-EX2)*(image_encrypt(x2(i)+1,y2(i)+1,k2(i))-EY2_DJX);
   
end
%统一在循环外除以像素点对数5000，可减少运算次数
%求方差和协方差
%R通道
DX2=DX2/Num_rand;
DY2_SP=DY2_SP/Num_rand;
DY2_CZ=DY2_CZ/Num_rand;
DY2_DJX=DY2_DJX/Num_rand;
COVXY2_SP=COVXY2_SP/Num_rand;
COVXY2_CZ=COVXY2_CZ/Num_rand;
COVXY2_DJX=COVXY2_DJX/Num_rand;

%求相关性系数
%水平、垂直、对角线的相关性
RXY2_SP=COVXY2_SP/sqrt(DX2*DY2_SP);
RXY2_CZ=COVXY2_CZ/sqrt(DX2*DY2_CZ);
RXY2_DJX=COVXY2_DJX/sqrt(DX2*DY2_DJX);

disp('加密后R、G、B通道相关性：');
disp(['加密后图片相关性：','  水平相关性=',num2str(RXY2_SP),'  垂直相关性=',num2str(RXY2_CZ),'  对角线相关性=',num2str(RXY2_DJX)]);
disp('');