function image_pollute = block_pollute(image,type)

%-- 丢失左下角1/16
if type == 1
    
    [M,N,~] = size(image);
    image_R = image(:,:,1);
    image_G = image(:,:,2);
    image_B = image(:,:,3);
    m = M/4;
    n = N/4;
    for i=1:M
        for j=1:N
            if (i>3*m)&&(j<=n)
                image_R(i,j) = 0;
                image_G(i,j) = 0;
                image_B(i,j) = 0;
            end
        end
    end
    
    image_pollute = zeros(M,N,3);
    image_pollute(:,:,1) = image_R;
    image_pollute(:,:,2) = image_G;
    image_pollute(:,:,3) = image_B;
    
else
    [M,N] = size(image);
    m = M/4;
    n = N/4;
    for i=1:M
        for j=1:N
            if (i>3*m)&&(j<=n)
                image(i,j) = 0;
            end
        end
    end
    
    image_pollute = image;

end

