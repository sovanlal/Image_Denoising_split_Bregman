%this function solves for d

% Author: Sovan Mukherjee, April 2014

function d= d_subproblem_updated(del,u,lamda,gamma,b,n)

B= del*u+b;

%% for anisotropic denoising

d= sign(B).*max(abs(B)-(lamda/gamma),0);  

