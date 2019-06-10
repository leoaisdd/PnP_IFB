function [out,PSNR] = PlugPlayADMM_deblur(y,h,raw,lambda,method,opts,M)

if nargin<4
    error('not enough input, try again \n');
elseif nargin==4
    opts = [];
end

% Check defaults
if ~isfield(opts,'rho')
    opts.rho = 1;
end
if ~isfield(opts,'max_itr')
    opts.max_itr = 20;
end
if ~isfield(opts,'tol')
    opts.tol = 1e-4;
end
if ~isfield(opts,'gamma')
    opts.gamma=1;
end
if ~isfield(opts,'print')
    opts.print = false;
end


% set parameters
% max_itr   = opts.max_itr;
max_itr   =2000;
tol       = opts.tol;
gamma     = opts.gamma;
rho       = opts.rho;
% M=100;
%initialize variables
dim         = size(y);          
N           = dim(1)*dim(2);    
Hty         = imfilter(y,h,'circular'); 
eigHtH      = abs(fftn(h, dim)).^2;
v           = 0.5*ones(dim);
% x           = v;
u           = zeros(dim);
residual    = inf;
   V = psf2otf(h,[dim(1),dim(2)]);
%set function handle for denoiser
switch method
    case 'BM3D'
        denoise=@wrapper_BM3D;
    case 'TV'
        denoise=@wrapper_TV;
    case 'NLM'
        denoise=@wrapper_NLM;
    case 'RF'
        denoise=@wrapper_RF;
    case 'DnCNN'
        denoise=@wrapper_DnCNN;
    case 'l0'
        denoise=@wrapper_l0;
    case 'nuclear_l21'
        denoise=@wrapper_nuclear_l21
    otherwise
        error('unknown denoiser \n');
end

% main loop

if opts.print==true
    fprintf('Plug-and-Play ADMM --- Deblurring \n');
    fprintf('Denoiser = %s \n\n', method);
%     fprintf('itr \t ||x-xold|| \t ||v-vold|| \t ||u-uold|| \n');
end
itr = 3;
% rho=6;
gamma=1/6;
alpha=0;
%   x=zeros(dim(1),dim(2),max_itr);
x(:,:,1)=v;
x(:,:,2)=v;
t=1;
% gamma=0.95;
while(itr<=max_itr)
 %% sIFB scheme

%      x_old = x(:,:,itr-1);
%      hh=1-y./(imfilter(x(:,:,itr-1), h, 'circular', 'conv')+eps);
%      adjKhh=real(ifft2(conj(V).*fft2(hh)));
%      z=x(:,:,itr-1)-gamma*adjKhh+alpha*(x(:,:,itr-1)-x(:,:,itr-2));
% 
%     sigma  = sqrt(lambda*gamma);
%     x(:,:,itr)   = denoise(z,sigma);
%     x(:,:,itr)  = proj(x(:,:,itr),[0,M]); 
   

%% uIFB scheme
%      alpha=1-3/itr;
      t0=t;
      t=sqrt(4*t0^2+1)/2;
     alpha=(t0-1)/t;
     x_old = x(:,:,itr-1);
     ss=x(:,:,itr-1)+alpha*(x(:,:,itr-1)-x(:,:,itr-2));
     hh=1-y./(imfilter(ss, h, 'circular', 'conv')+eps);
     adjKhh=real(ifft2(conj(V).*fft2(hh)));
     z=ss-gamma*adjKhh;

    sigma  = sqrt(lambda*gamma);
    x(:,:,itr)   = denoise(z,sigma);
    x(:,:,itr)  = proj(x(:,:,itr),[0,M]); 
%%

    PSNR(itr) = psnr(double(x(:,:,itr)),double(raw),M);
     if opts.print==true
         fprintf('%3g \t %3.5e \n', itr, PSNR(itr));
     end

    itr = itr+1;
end
out = x(:,:,itr-1);
end