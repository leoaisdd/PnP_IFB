
clear all
close all
clc
 rng('default')
addpath(genpath('./utilities/'));

%add path to denoisers
addpath(genpath('./denoisers/BM3D/'));
addpath(genpath('./denoisers/TV/'));
addpath(genpath('./denoisers/NLM/'));
addpath(genpath('./denoisers/RF/'));
addpath(genpath('./denoisers/DnCNN/'));
addpath(genpath('./denoisers/gradientl0/'));
addpath(genpath('./denoisers/nuclear_l21/'));
%read test image
z = im2double(imread('./data/House256.png'));
% z = im2double(imread('./data/Cameraman256.png'));
% z = im2double(imread('./data/Lena512.png'));
%  S=50;
%  z=simuyuan(S);
%initialize a blurring filter
h = fspecial('gaussian',[9 9],2);
% h = fspecial('gaussian',10,10);
h=h/sum(h(:));

%reset random number generator
rng(0);

%set noies level
noise_level = 10/255;

M=100;
z11 = M*z/max(z(:));
y11 = imfilter(z11, h, 'circular', 'conv');

 y=poissrnd(y11);
% y=y/max(y(:));
%calculate observed image
% y = imfilter(z,h,'circular')+noise_level*randn(size(z));
% y = proj(y,[0,1]);
tic
%parameters
method = 'nuclear_l21';
switch method
    case 'RF'
        lambda = 0.1;
    case 'NLM'
        lambda = 0.005;
    case 'BM3D'
        lambda = 0.001;
    case 'TV'
        lambda = 0.01;
    case 'DnCNN'
        lambda =100;
    case 'l0'
        lambda=0.005;
%         rho=60*M/lambda;
    case 'nuclear_l21'
        lambda=0.05;
end

%optional parameters
opts.rho     = 1;
opts.gamma   = 1;
opts.max_itr = 20;
opts.print   = true;

%main routine
[out,PSNR] = PlugPlayproposed_deblur(y,h,z11,lambda,method,opts,M);
toc
%display
PSNR_output = psnr(out,z11,M);

% fprintf('\nPSNR = %3.2f dB \n', PSNR_output);

figure(1);
subplot(121);
imshow(y,[0,max(y(:))]);
title('Input');

subplot(122);
imshow(out,[0,max(out(:))]);
tt = sprintf('PSNR = %3.2f dB', PSNR_output);
title(tt);
figure(2)
 plot(PSNR(3:end))
% imshow(y,[0,max(y(:))])