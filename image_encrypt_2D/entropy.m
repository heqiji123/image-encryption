function [entropy] = entropy(image)

T1=imhist(uint8(image));                        %统计图像R通道灰度值从0~255的分布情况，存至T1
S1=sum(T1);                                     %计算整幅图像R通道的灰度值
entropy=0;                                      %原始图片R通道相关性

for i=1:length(T1)
    pp1=T1(i)/S1;                               %每个灰度值占比，即每个灰度值的概率
    if pp1~=0
        entropy=entropy-pp1*log2(pp1);          %求信息熵
    end
end

