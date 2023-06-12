
function reconstructed_image = BCS_ED(y, Phi,num_rows, num_cols, num_levels)

lambda = 15;                        
max_iterations = 200;           
[~,N] = size(Phi);                  
block_size = sqrt(N);              
TOL = 0.00001;                        
D_prev = 0;    

x = Phi' * y;                          

for i = 1:max_iterations
           
  [x, D] = SPLIteration(y,x,Phi,block_size,num_rows,num_cols,lambda,num_levels);

if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))
          break;   
end  
  D_prev = D; 
end

[x,~] = SPLIteration(y,x,Phi,block_size,num_rows,num_cols,lambda,num_levels);

x = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distict');                                

reconstructed_image = x;                              

end



function [x, D] = SPLIteration(y,x,Phi,block_size,num_rows,num_cols,lambda,num_levels)

[af, sf] = farras;
L = 12;


x = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct');                                                  

x_hat = wiener2(x, [3, 3]);                                   %-- winner检测
x_hat = im2col(x_hat, [block_size block_size], 'distinct');
x_hat = x_hat + Phi' * (y - Phi * x_hat);                     %最小二乘法 

x1 = col2im(x_hat, [block_size block_size], ...
    [num_rows num_cols], 'distinct');                        

% x_check = waveletcdf97(x1, num_levels);                                             
% threshold = lambda  * sqrt(2 * log(num_rows * num_cols)) * (median(abs(x_check(:))) / 0.6745);
% x_check(abs(x_check) < threshold) = 0;                    
% x_bar = waveletcdf97(x_check , -num_levels);                  %-- 进行阈值处理
% x_bar = im2col(x_bar, [block_size block_size], 'distinct');       


x_check = dwt2D(symextend(x1, L * 2^(num_levels - 1)), ...
    num_levels, af);

if (nargin == 9)
  end_level = 1;
else
  end_level = num_levels - 1;
end
x_check = SPLBivariateShrinkage(x_check, end_level, lambda);

x_bar = idwt2D(x_check, num_levels, sf);
Irow = (L * 2^(num_levels - 1) + 1):(L * 2^(num_levels - 1) + num_rows);
Icol = (L * 2^(num_levels - 1) + 1):(L * 2^(num_levels - 1) + num_cols);
x_bar = x_bar(Irow, Icol);
x_bar = im2col(x_bar, [block_size block_size], 'distinct');


x = x_bar + Phi' * (y - Phi * x_bar);                         %最小二乘法

x2 = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct');        
  
D = RMS(x1, x2);                   

end
  
function r = RMS(x1, x2)

error = x1 - x2;

r = sqrt(mean(error(:).^2));
end


function x_check = SPLBivariateShrinkage(x_check, end_level, lambda)

windowsize  = 3;
windowfilt = ones(1, windowsize)/windowsize;

tmp = x_check{1}{3};
Nsig = median(abs(tmp(:)))/0.6745;

for scale = 1:end_level
  for dir = 1:3
    Y_coefficient = x_check{scale}{dir};
    
    Y_parent = x_check{scale+1}{dir};
    
    Y_parent = expand(Y_parent);
    
    Wsig = conv2(windowfilt, windowfilt, (Y_coefficient).^2, 'same');
    Ssig = sqrt(max(Wsig-Nsig.^2, eps));
    
    T = sqrt(3)*Nsig^2./Ssig;
    
    x_check{scale}{dir} = bishrink(Y_coefficient, ...
	Y_parent, T*lambda);
  end
end

end




