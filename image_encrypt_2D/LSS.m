function P = LSS(x0,r0,N0,N)
P = zeros(1,N0+1+N);
P(1) = x0;
for i=2:(N0+N+1)
    P(i) = mod((r0*P(i-1)*(1-P(i-1))+(4-r0)*sin(pi*P(i-1))/4),1);
end
P = P(1,N0+2:length(P));
end
