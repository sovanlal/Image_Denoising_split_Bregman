% This function calculates gradient and DELTA oprator

% Author: Sovan Mukherjee, April, 2014

function [del,del_x,del_y,del_star,delta]=grad_operator(N)
 
% gradient operator
D= spdiags([-ones(N,1) ones(N,1)],[0,1],N,N);

% D(N,:)=0;

del_x = kron(speye(N),D);                                                  % gradient matrix in horizontal direction

del_y= kron(D, speye(N));                                                  % gradient matrix in vertical direction

del= [del_x;del_y];                                                        %gradient                                                       

% transpose of gradient 

del_star= del';

% delta operator

delta= -del_star*del;