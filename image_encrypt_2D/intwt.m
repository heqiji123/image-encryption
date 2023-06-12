% integar wavelet transform
%using 9-7 filter
%signal symmetry extention mathod:(2,2)  ...CBAABCDEFFED...
%ss  dd
%sd  ds 
%clear all;
%fid=fopen('f:\data\aviris\aviris69.raw','rb');
%status=fseek(fid,512*512*2*40,'bof');
row=256;
col=256;
%orin=double(fread(fid,[row col],'uint16'));
orin=ss;
%%%%%%%%%%%%%%%%%% COEFFECIENT DEFINITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=-1.586134342;
b=-0.05298011854;
r=0.8829110762;
d=0.4435068522;
%k=1.149604398;
%%%%%%%%%%%%%%%%%%%========================%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%init
d1_1=zeros(row,col/2);                    %high frequency components
d1=zeros(row,col/2);
s1_1=zeros(row,col/2);                    %low frequency components
s1=zeros(row,col/2);

sd_1=zeros(row/2,col/2);                    %high frequency components
sd=zeros(row/2,col/2);
dd_1=zeros(row/2,col/2);
dd=zeros(row/2,col/2);
ds_1=zeros(row/2,col/2);
ds=zeros(row/2,col/2);
ss_1=zeros(row/2,col/2);                    %low frequency components
ss=zeros(row/2,col/2);
%%%%%%%%%%%%%%%%%%==============================%%%%%%%%%%%%%%%%%%%%%%%%
%%%border extension

%LINE TRANSFORM 
d1_1(:,col/2)=orin(:,col-1)+round(a*(orin(:,col)+orin(:,col-1)));
for j=1:(col/2-1)
   d1_1(:,j)=orin(:,2*j-1)+round(a*(orin(:,2*j)+orin(:,2*j+2)));
end

s1_1(:,1)=orin(:,2)+round(b*(d1_1(:,1)+d1_1(:,1)));
for j=2:col/2
   s1_1(:,j)=orin(:,2*j)+round(b*(d1_1(:,j)+d1_1(:,j-1)));
end

d1(:,col/2)=d1_1(:,col/2)+round(r*(s1_1(:,col/2)+s1_1(:,col/2)));
for j=1:(col/2-1)
   d1(:,j)=d1_1(:,j)+round(r*(s1_1(:,j)+s1_1(:,j+1)));
end
s1(:,1)=s1_1(:,1)+round(d*d1(:,1)*2);
for j=2:col/2
   s1(:,j)=s1_1(:,j)+round(d*(d1(:,j)+d1(:,j-1)));
end
%s1=round(s1*k);
%d1=round(d1/k);
figure;imshow(uint8(s1'))
figure;imshow(uint8(d1'))

%col TRANSFORM
sd_1(row/2,:)=s1(row-1,:)+round(a*(s1(row,:)+s1(row-1,:)));
for i=1:(row/2-1)
    sd_1(i,:)=s1(2*i-1,:)+round(a*(s1(2*i,:)+s1(2*i+2,:)));
end
ss_1(1,:)=s1(2,:)+round(b*(sd_1(1,:)+sd_1(1,:)));
for i=2:row/2
   ss_1(i,:)=s1(2*i,:)+round(b*(sd_1(i,:)+sd_1(i-1,:)));
end
sd(row/2,:)=sd_1(row/2,:)+round(r*(ss_1(row/2,:)+ss_1(row/2,:)));
for i=1:(row/2-1)
   sd(i,:)=sd_1(i,:)+round(r*(ss_1(i,:)+ss_1(i+1,:)));
end
ss(1,:)=ss_1(1,:)+round(d*sd(1,:)*2);
for i=2:row/2
   ss(i,:)=ss_1(i,:)+round(d*(sd(i,:)+sd(i-1,:)));
end
%ss=round(ss*k);
%sd=round(sd/k);

dd_1(row/2,:)=d1(row-1,:)+round(a*(d1(row,:)+d1(row-1,:)));
for i=1:(row/2-1)
    dd_1(i,:)=d1(2*i-1,:)+round(a*(d1(2*i,:)+d1(2*i+2,:)));
end
ds_1(1,:)=d1(2,:)+round(b*(dd_1(1,:)+dd_1(1,:)));
for i=2:row/2
   ds_1(i,:)=d1(2*i,:)+round(b*(dd_1(i,:)+dd_1(i-1,:)));
end
dd(row/2,:)=dd_1(row/2,:)+round(r*(ds_1(row/2,:)+ds_1(row/2,:)));
for i=1:(row/2-1)
   dd(i,:)=dd_1(i,:)+round(r*(ds_1(i,:)+ds_1(i+1,:)));
end
ds(1,:)=ds_1(1,:)+round(d*dd(1,:)*2);
for i=2:row/2
   ds(i,:)=ds_1(i,:)+round(d*(dd(i,:)+dd(i-1,:)));
end
%dd=round(dd*k);
%ds=round(ds/k);
Max1=max(max(ss))
Min1=min(min(ss))
b=255*(ss-Min1)/(Max1-Min1);
figure;imshow(uint8(b'))
Max2=max(max(sd))
Min2=min(min(sd))
Max3=max(max(dd))
Min3=min(min(dd))
Max4=max(max(ds))
Min4=min(min(ds))
b=255*(sd-Min2)/(Max2-Min2);
figure;imshow(uint8(b'))
b=255*(sd-Min3)/(Max3-Min3);
figure;imshow(uint8(b'))
b=255*(dd-Min4)/(Max4-Min4);
figure;imshow(uint8(b'))

%figure;imshow(uint8(sd'),[Min2,Max2])
%figure;imshow(uint8(dd'),[Min3,Max3])
%figure;imshow(uint8(ds'),[Min4,Max4])
%fclose(fid);





   
