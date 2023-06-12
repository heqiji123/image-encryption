function S = index_permutation(PS,X)
%-- PS - block_size*block_size
%-- X - 1*(block_size*block_size) 
%-- S - block_size*block_size 

    [~,T] = sort(X);
    [M,N] = size(PS);
    PS = reshape(PS,1,M*N);
    S = zeros(1,M*N);
    for i=1:M*N
        S(i) = PS(T(i));
    end
    
    S = reshape(S,M,N);
end

