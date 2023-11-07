clear all;
close all;
clc;
% ///////////// Parameter Setting /////////////////
A1=imread('images.jpg');
A2 = im2double(A1)
A2 = A2(:,:,1)
A = 255 .* A2;
m = size(A,1);
n = size(A,2);
r =8; % Rank of matrix
C =8; % Specify number of samples: C times of dimension.


maxIter =30; % Maximum number of iterations for AltMin method
tol = 0.0001; % Tolerance of difference below which AltMin will stop iteration
show = true; % Iteration log will show if this is true


% ////////////// Running Area ///////////////////
NumOfScale = m * n;
NumOfDimension = m*r + n*r - r^2;
p = C * NumOfDimension / NumOfScale; % Sampling probability
MASK = zeros(m,n);
MASK = RandomPick(MASK,p);
MASK_ = ~MASK
NumOfSample = sum(MASK(:));
AA = A .* MASK;
AA1 = AA ./ 255;
tic;
[Acomplete,difference] = Altmin(AA1,MASK,r,maxIter,tol,show)
toc;
diffALL = norm(MASK_.*(A- Acomplete),'fro')/norm(MASK_.*A,'fro') % for unmasked entries
TM = max(max(Acomplete));
Acomplete = Acomplete./TM;

subplot(1,3,1);imshow(A2);title('Original Image')
subplot(1,3,2);imshow(AA1);title('Masked Image')
subplot(1,3,3);imshow(Acomplete);title('Restored image')