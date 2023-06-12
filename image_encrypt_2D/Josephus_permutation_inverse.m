function S_ = Josephus_permutation_inverse(S,start,step)

[M,N] = size(S);
S_ = zeros(M,N);
T = Josephus(N,start,step);
for i = 1:N
    S_(:,T(i)) = S(:,i);
end
end

