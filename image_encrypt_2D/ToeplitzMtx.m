function [ Phi ] = ToeplitzMtx( M,N )
%ToeplitzMtx Summary of this function goes here
%   Generate Toeplitz matrix 
%   M -- RowNumber
%   N -- ColumnNumber
%   Phi -- The Toeplitz matrix

%% Generate a random vector
%     %(1)Gauss
%     u = randn(1,2*N-1);
    %(2)Bernoulli
    u = randi([0,1],1,2*N-1);
    u(u==0) = -1;
%% Generate Toeplitz matrix   托普利兹矩阵
    Phi_t = toeplitz(u(N:end),fliplr(u(1:N)));
    Phi = Phi_t(1:M,:);
end