function PS = index_permutation_inverse(S,X)
%S -- block_size * block_size
%X -- 1*(block_size * block_size)
%PS -- block_size * block_size

    [M,N] = size(S);
    S = reshape(S,1,M*N);
    [~,T] = sort(X);
    PS = zeros(1,M*N);
    for i=1:M*N
        PS(T(i)) = S(i);
    end
    
    PS = reshape(PS,M,N);

end

