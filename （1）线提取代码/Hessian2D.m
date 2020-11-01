function [Dx,Dy,Dxx,Dxy,Dyy] = Hessian2D(I,Sigma)
%  This function Hessian2 Filters the image with 2nd derivatives of a 
%  Gaussian with parameter Sigma.
% 
% [Dxx,Dxy,Dyy] = Hessian2(I,Sigma);
% 
% inputs,
%   I : The image, class preferable double or single
%   Sigma : The sigma of the gaussian kernel used
%
% outputs,
%   Dxx, Dxy, Dyy: The 2nd derivatives
%
% example,
%   I = im2double(imread('moon.tif'));
%   [Dxx,Dxy,Dyy] = Hessian2(I,2);
%   figure, imshow(Dxx,[]);
%
% Function is written by D.Kroon University of Twente (June 2009)

if nargin < 2, Sigma = 1; end

% Make kernel coordinates
[X,Y]   = ndgrid(-round(3*Sigma):round(3*Sigma)); %高斯卷积核就是这么生成的，加上下面的话，
%也可以用for循环，比如for i=1:5,for j=1:5,这个5是卷积核的size，这里面是用sigma来确定size的，
% Build the gaussian 2nd derivatives filters
DGaussx  = 1/(2*pi*Sigma^4)*(-X).* exp(-(X.^2 + Y.^2)/(2*Sigma^2)); %二维高斯函数一阶对x的偏导
DGaussy  = 1/(2*pi*Sigma^4)*(-Y).* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
DGaussxx = 1/(2*pi*Sigma^4) * (X.^2/Sigma^2 - 1) .* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
DGaussxy = 1/(2*pi*Sigma^6) * (X .* Y)           .* exp(-(X.^2 + Y.^2)/(2*Sigma^2));
DGaussyy = DGaussxx';

Dx  = imfilter(I,DGaussx,'conv');
Dy  = imfilter(I,DGaussy,'conv');
Dxx = imfilter(I,DGaussxx,'conv');
Dxy = imfilter(I,DGaussxy,'conv');
Dyy = imfilter(I,DGaussyy,'conv');
