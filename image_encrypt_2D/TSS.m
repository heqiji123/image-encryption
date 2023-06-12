function T = TSS(x0,r0,N0,N)
%--- TSS混沌系统
%--- x0 r0初始值
%--- N0需要丢弃的前N0个值
%--- N迭代的长度
T = zeros(1,N0+1+N);
T(1,1) = x0;

for i=2:(N0+1+N)
    if T(i-1)<0.5
        T(i) = mod(r0*T(i-1)/2+(4-r0)*sin(pi*T(i-1))/4,1);
    else
        T(i) = mod(r0*(1-T(i-1))/2+(4-r0)*sin(pi*T(i-1))/4,1);
    end
end
T = T(1,N0+2:length(T));
end

