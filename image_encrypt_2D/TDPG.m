
function reconstructed_image = TDPG (Y, Phi1, Phi2, num_rows, num_cols, num_levels)

addpath('C:\Users\Thermos_121380\matlab\github\2DCS-ETC\TDCS-ETC\TDCS-ETC\mywork');                      %    functions


lambda = 20;                          %  the convergence-control factor
max_iterations = 200;                 %  maximum iteration number
TOL = 0.00001;                        %  the given error tolerance
D_prev = 0;                             

pinvPhi1 = pinv(Phi1);
pinvPhi2 = pinv(Phi2');
X = pinvPhi1 * Y * pinvPhi2;         %    Initialization

for i = 1:max_iterations
      
 [X, D] = SPLIterationBivariate (Y,X,Phi1,Phi2,pinvPhi1,pinvPhi2,num_rows,num_cols,lambda,num_levels);                   

 % if the difference between the current solution and the previous solution is smaller 
 % than the given error tolerance, break.
 
  if ((D_prev ~= 0) && (abs(D - D_prev) < TOL))       
          break;      
  end
  D_prev = D;
end


[X,~] = SPLIterationBivariate (Y,X,Phi1,Phi2,pinvPhi1, pinvPhi2,num_rows,num_cols,lambda,num_levels);
reconstructed_image = X;
    


function [X, D] = SPLIterationBivariate(Y,X,Phi1,Phi2,pinvPhi1,pinvPhi2,num_rows,num_cols,lambda,num_levels)
[af, sf] = farras;
L = 12;
r = 0.8;                                    %   the step-size
X1 = X;                                     %    Save the previous solution

% Gradient descent in spatial domain
X_hat = X - r * derivating_of_TV(X, num_rows);     

% Bivariate shrinkage in wavelet domain
X_check = dwt2D(symextend(X_hat, L * 2^(num_levels - 1)), ...
    num_levels, af);      
 
if (nargin == 9)                          
  end_level = 1;
else
  end_level = num_levels - 1;
end

X_check = SPLBivariateShrinkage(X_check, end_level, lambda);
X_bar = idwt2D(X_check, num_levels, sf);   
Irow = (L * 2^(num_levels - 1) + 1):(L * 2^(num_levels - 1) + num_rows);
Icol = (L * 2^(num_levels - 1) + 1):(L * 2^(num_levels - 1) + num_cols);
X_bar = X_bar(Irow, Icol);

% Projection
X = X_bar +  pinvPhi1 * (Y-Phi1 *X_bar * Phi2')* pinvPhi2; 

D = RMS(X, X1);
 
% Bivariate shrinkage in wavelet domain

% x_check = waveletcdf97(X_hat, num_levels);     
% threshold = lambda * sqrt(2 * log(num_rows * num_cols)) * (median(abs(x_check(:)))/0.6745); 
% x_check(abs(x_check)<threshold) = 0;
% X_bar = waveletcdf97(x_check , -num_levels);      
% 
% % Projection
% X = X_bar +  pinvPhi1 * (Y-Phi1 *X_bar * Phi2')* pinvPhi2;
% 
% D = RMS(X, X1);

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


