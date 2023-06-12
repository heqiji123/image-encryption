%�ļ�����liftwavedec2.m
%�������ܣ�����С��harr�任������С���任
%�����ʽ������imgwave=harr_liftwavedec2(image,256,256,3)
%����˵����
% image--�����ͼ�����ҪΪ����
% m n   --�����ͼ������С
% count    --С���任����
%����������
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
% ����С��harr�任������С���任
%
function f_row=lwavedec2(image,M,N)
f=image;
S=M/2;
T=N/2;               %  ��ͼ��ά��


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   1.�б任

%  A.���ѣ���ż�ֿ���

f1=f([1:2:M-1],:);  %  ����
f2=f([2:2:M],:);    %  ż��

% f1(:,T+1)=f1(:,1);  %  ����
% f2(T+1,:)=f2(1,:);  %  ����

%  B.Ԥ��

for i_hc=1:S;
    high_frequency_column(i_hc,:)=f1(i_hc,:)-f2(i_hc,:);
end;

% high_frequency_column(T+1,:)=high_frequency_column(1,:);  %  ����

%  C.����

for i_lc=1:S;
    low_frequency_column(i_lc,:)=f2(i_lc,:)+1/2*high_frequency_column(i_lc,:);
end;

%  D.�ϲ�
f_column([1:1:S],:)=low_frequency_column([1:S],:);
f_column([S+1:1:M],:)=high_frequency_column([1:S],:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   2.�б任

%  A.���ѣ���ż�ֿ���

f1=f_column(:,[1:2:N-1]);  %  ����
f2=f_column(:,[2:2:N]);    %  ż��


% f2(:,T+1)=f2(:,1);    %  ����

%  B.Ԥ��

for i_hr=1:T;
    high_frequency_row(:,i_hr)=f1(:,i_hr)-f2(:,i_hr);
end;

% high_frequency_row(:,T+1)=high_frequency_row(:,1);  %  ����

%  C.����

for i_lr=1:T;
    low_frequency_row(:,i_lr)=f2(:,i_lr)+1/2*high_frequency_row(:,i_lr);
end;

%  D.�ϲ�
f_row(:,[1:1:T])=low_frequency_row(:,[1:T]);
f_row(:,[T+1:1:N])=high_frequency_row(:,[1:T]);
end