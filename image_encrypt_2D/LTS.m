function T = LTS(x0,r0,N0,N)
% LTS混沌系统
% x0 r0初始值  
% N0 舍去的长度
% N 混沌系统的长度
T = zeros(1,N0+1+N);
T(1) = x0;
for i=2:(N0+N)+1
    if T(i-1)<0.5
        T(i) = mod(r0*T(i-1)*(1-T(i-1))+(4-r0)*T(i-1)/2,1);
    end
    if T(i-1)>=0.5
        T(i) = mod(r0*T(i-1)*(1-T(i-1))+(4-r0)*(1-T(i-1))/2,1);
    end
end
T = T(1,N0+2:length(T));
end

