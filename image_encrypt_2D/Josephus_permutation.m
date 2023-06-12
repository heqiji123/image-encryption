function S_ = Josephus_permutation(S,start,step)

[M,N] = size(S);
S_ = zeros(M,N);
T = Josephus(N,start,step);
for i = 1:N
    S_(:,i) = S(:,T(i));
end
end

