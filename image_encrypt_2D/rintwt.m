%%%%---------process of the reconstruction of the image-----------%%%%
%using 9-7 filter
%signal symmetry extention mathod:(2,2)  ...CBAABCDEFFED...
%ss  dd
%sd  ds 

%init
clear all;
row=512;
column=512;
a=-1.586134342;
b=-0.05298011854;
r=0.8829110762;
d=0.4435068522;

load intwcoef;
reconst=zeros(512,512);

d2=zeros(512,256);                 %high frequency components
d2_2=zeros(512,256);
s2=zeros(512,256);                 %low frequency components
s2_2=zeros(512,256);

sd_2=zeros(256,256);                    %high frequency components
dd_2=zeros(256,256);
ds_2=zeros(256,256);
ss_2=zeros(256,256);                    %low frequency components


for i=2:row/2
   ds_2(i,:)=ds(i,:)-round(d*(dd(i,:)+dd(i-1,:)));
end
ds_2(1,:)=ds(1,:)-round(d*dd(1,:)*2);

for i=1:(row/2-1)
   dd_2(i,:)=dd(i,:)-round(r*(ds_2(i,:)+ds_2(i+1,:)));
end
dd_2(row/2,:)=dd(row/2,:)-round(r*(ds_2(row/2,:)+ds_2(row/2,:)));

for i=2:row/2
   d2_2(2*i,:)=ds_2(i,:)-round((dd_2(i,:)+dd_2(i-1,:))*b);
end
d2_2(2,:)=ds_2(1,:)-round((dd_2(1,:)+dd_2(1,:))*b);

for i=1:(row/2-1)
    d2_2(2*i-1,:)=dd_2(i,:)-round((d2_2(2*i,:)+d2_2(2*i+2,:))*a);
 end
d2_2(row-1,:)=dd_2(row/2,:)-round((d2_2(row,:)+d2_2(row-1,:))*a);
 
 

ss_2(1,:)=ss(1,:)-round(d*sd(1,:)*2);
for i=2:row/2
   ss_2(i,:)=ss(i,:)-round(d*(sd(i,:)+sd(i-1,:)));
end
sd_2(row/2,:)=sd(row/2,:)-round(r*(ss_2(row/2,:)+ss_2(row/2,:)));
for i=1:(row/2-1)
   sd_2(i,:)=sd(i,:)-round(r*(ss_2(i,:)+ss_2(i+1,:)));
end
s2_2(2,:)=ss_2(1,:)-round((sd_2(1,:)+sd_2(1,:))*b);
for i=2:row/2
   s2_2(2*i,:)=ss_2(i,:)-round((sd_2(i,:)+sd_2(i-1,:))*b);
end
s2_2(row-1,:)=sd_2(row/2,:)-round((s2_2(row,:)+s2_2(row-1,:))*a);
for i=1:(row/2-1)
    s2_2(2*i-1,:)=sd_2(i,:)-round((s2_2(2*i,:)+s2_2(2*i+2,:))*a);
 end
 figure;imshow(uint8(s2_2))
 figure;imshow(uint8(d2_2))
 %%%%---------------------------------------------------%%%%% 
s2(:,1)=s2_2(:,1)-round(d*d2_2(:,1)*2);
for j=2:column/2
   s2(:,j)=s2_2(:,j)-round(d*(d2_2(:,j)+d2_2(:,j-1)));
end
d2(:,column/2)=d2_2(:,column/2)-round(r*(s2(:,column/2)+s2(:,column/2)));
for j=1:(column/2-1)
   d2(:,j)=d2_2(:,j)-round(r*(s2(:,j)+s2(:,j+1)));
end
reconst(:,2)=s2(:,1)-round(b*(d2(:,1)+d2(:,1)));
for j=2:column/2
   reconst(:,2*j)=s2(:,j)-round(b*(d2(:,j)+d2(:,j-1)));
end
reconst(:,column-1)=d2(:,column/2)+round(a*(reconst(:,column)+reconst(:,column-1)));
for j=1:(column/2-1)
   reconst(:,2*j-1)=d2(:,j)-round(a*(reconst(:,2*j)+reconst(:,2*j+2)));
end
Max1=max(max(reconst))
Min1=min(min(reconst))
b=255*(reconst-Min1)/(Max1-Min1);
figure;imshow(uint8(b'))

%figure;imshow(uint8(reconst'))
