%% Split Bregman Isotropic Total Variation(TV)denoising

%Ref: Goldstein et al.,"The Split Bregman Method for L1 
% Regularized Problems",SIAM Journal on Imaging Sciences, Vol. 2,no. 2
% pp 323-343, 2009.

% Author: Sovan Mukherjee, May, 2014.

clear all;
close all;
clc;

%% initialize parameters

lamda=0.10;

gamma=1;

%% Shepp-Logan phantom

N = 1024;                                                                  % no. of pixels

im =phantom(N);

sigma= 0.10;                                                               %standard deviation of the noise
% 
f = im+sigma*randn(size(im,1));                                            % add noise

%% ------------------
figure;imshow(im); title('Original Image');colormap('jet');                % show the original image

figure; imshow(f); title('Noisy Image'); colormap('jet');                  % show the noisy image

%% -------------------------------------

n=N^2;

im= im(:);
f=f(:);

TOL =1e-3;                                                                 % tolerance for split-Bregman iteration

%% Split-Bregman TV denoising algorithm
u_init=f;                                                                  %initialize u(denoised image)with f(noisy image)

[del,D1,D2,del_star,delta]=grad_operator(N);                               %calculate gradient matrix

d= zeros(2*n,1);

b= d;

err= 1;
it=1;

%% Run the iteration

tic;

while err >TOL

%% u subproblem

[u_update,~]= cgs(speye(N*N)-delta, ...
    f+gamma*del_star*(d-b),1e-5,100); 

fprintf('iteration# %d\t',it);

%% d subproblem

d= d_subproblem_updated(del,u_update,lamda,gamma,b,n);

%% update b

b= b+del*u_update-d;
% iter_term=norm((u_update-u_init),2)

err=norm((u_update-u_init),2)./...                                         %iteration tolerance
    norm(u_update,2);

u_init=u_update;                                                           %update initial u

fprintf('err=%d\n',err);

it=it+1;

end

t.tv=toc;

fprintf('\nReconstruction time (s)= %f\n',t.tv);

fprintf('\nIteration stopped\n');

fprintf('\nRelative error is= %d\n', err);

%% show denoised image
u_update=reshape(u_update,N,N);
figure;
imshow(u_update);title('denoised image');

