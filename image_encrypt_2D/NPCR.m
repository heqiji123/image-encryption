function npcr = NPCR(image1,image2)
[~,N,n] = size(image1);
D = 0;
for i = 1:N
    for j = 1:N
        for k = 1:n
            if image1(i,j,k) ~= image2(i,j,k)
                D = D + 1;
            end
        end
    end
end
npcr = D / (N * N * n);
end

