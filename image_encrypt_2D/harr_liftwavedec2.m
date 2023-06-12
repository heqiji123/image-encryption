%文件名：liftwavedec2.m
%函数功能：二代小波harr变换，整数小波变换
%输入格式举例：imgwave=harr_liftwavedec2(image,256,256,3)
%参数说明：
% image--输入的图像矩阵，要为方阵
% m n   --输入的图像矩阵大小
% count    --小波变换次数
%测试用例：
% img=imread('lena.jpg');
% [m,n]=size(img);
% imgwave=harr_liftwavedec2(img,m,n,3);
% imshow(imgwave);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imgwave=harr_liftwavedec2(image,m,n,count)
img=double(image);
M=m;
N=n;
for i=1:count
     imgwave1=lwavedec2(img,M,N);
     imgwave(1:M,1:N)=imgwave1;
     M=M/2;
	 N=N/2;
     img=imgwave1(1:M,1:N);
end
end
%
% 二代小波harr变换，整数小波变换
%
function f_row=lwavedec2(image,M,N)
f=image;
S=M/2;
T=N/2;               %  子图像维数


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   1.列变换

%  A.分裂（奇偶分开）

f1=f([1:2:M-1],:);  %  奇数
f2=f([2:2:M],:);    %  偶数

% f1(:,T+1)=f1(:,1);  %  补列
% f2(T+1,:)=f2(1,:);  %  补行

%  B.预测

for i_hc=1:S;
    high_frequency_column(i_hc,:)=f1(i_hc,:)-f2(i_hc,:);
end;

% high_frequency_column(T+1,:)=high_frequency_column(1,:);  %  补行

%  C.更新

for i_lc=1:S;
    low_frequency_column(i_lc,:)=f2(i_lc,:)+1/2*high_frequency_column(i_lc,:);
end;

%  D.合并
f_column([1:1:S],:)=low_frequency_column([1:S],:);
f_column([S+1:1:M],:)=high_frequency_column([1:S],:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   2.行变换

%  A.分裂（奇偶分开）

f1=f_column(:,[1:2:N-1]);  %  奇数
f2=f_column(:,[2:2:N]);    %  偶数


% f2(:,T+1)=f2(:,1);    %  补行

%  B.预测

for i_hr=1:T;
    high_frequency_row(:,i_hr)=f1(:,i_hr)-f2(:,i_hr);
end;

% high_frequency_row(:,T+1)=high_frequency_row(:,1);  %  补行

%  C.更新

for i_lr=1:T;
    low_frequency_row(:,i_lr)=f2(:,i_lr)+1/2*high_frequency_row(:,i_lr);
end;

%  D.合并
f_row(:,[1:1:T])=low_frequency_row(:,[1:T]);
f_row(:,[T+1:1:N])=high_frequency_row(:,[1:T]);
end