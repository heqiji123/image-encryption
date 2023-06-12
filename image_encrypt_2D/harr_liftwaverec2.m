function imgwave=harr_liftwaverec2(image,m,n,count)
M=m/(2^count);
N=n/(2^count);
imgwave=image;
for i=count:-1:1
     M=M*2;
     N=N*2;
     img=imgwave(1:M,1:N);
     imgwave1=lwaverec2(img,M,N);
     imgwave(1:M,1:N)=imgwave1;
end
imgwave = uint8(imgwave);
end

function f_column=lwaverec2(f_row,M,N)
S=M/2;
T=N/2;
%figure(2);
%subplot(221),imshow(f_row),title('变换后图像');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%反变换%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   1.行变换

%   A.提取（低频高频分开）

f1=f_row(:,[T+1:1:N]);  %  奇数
f2=f_row(:,[1:1:T]);    %  偶数


% f2(:,T+1)=f2(:,1);    %  补行

%  B.更新

for i_lr=1:T;
    low_frequency_row(:,i_lr)=f2(:,i_lr)-1/2*f1(:,i_lr);
end;

%  C.预测

for i_hr=1:T;
    high_frequency_row(:,i_hr)=f1(:,i_hr)+low_frequency_row(:,i_hr);
end;

% high_frequency_row(:,T+1)=high_frequency_row(:,1);  %  补行


%  D.合并(奇偶分开合并)
f_row(:,[2:2:N])=low_frequency_row(:,[1:T]);
f_row(:,[1:2:N-1])=high_frequency_row(:,[1:T]);   
%subplot(223),imshow(f_row),title('行变换图像');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   2.列变换

%  A.提取（低频高频分开）

f1=f_row([S+1:1:M],:);  %  奇数
f2=f_row([1:1:S],:);    %  偶数

% f1(:,T+1)=f1(:,1);  %  补列
% f2(T+1,:)=f2(1,:);  %  补行

%  B.更新

for i_lc=1:S;
    low_frequency_column(i_lc,:)=f2(i_lc,:)-1/2*f1(i_lc,:);
end;

%  C.预测

for i_hc=1:S;
    high_frequency_column(i_hc,:)=f1(i_hc,:)+low_frequency_column(i_hc,:);
end;

% high_frequency_column(T+1,:)=high_frequency_column(1,:);  %  补行

%  D.合并(奇偶分开合并）
f_column([2:2:M],:)=low_frequency_column([1:S],:);
f_column([1:2:M-1],:)=high_frequency_column([1:S],:);
    
%subplot(224),imshow(f_column),title('行变换图像');
end
