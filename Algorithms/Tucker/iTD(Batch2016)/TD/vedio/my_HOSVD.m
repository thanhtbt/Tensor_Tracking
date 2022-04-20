function [G U] = my_HOSVD(X, J)
%
% hosvd - Higher-Order Singular Value Decomposition with truncation.
%           
%
% DESCRIPTION:
% 
%      ([1], 7.118, p. 362)
%
%          X = G x1 U1 x2 U2 ... xN UN
%
%
% INPUTS:
%
%   X - Data tensor, order-N.
%   J - N dimensional vector, determines number of eigenvectors for each 
%       mode.
%   G - Core tensor.
%
%
%
% REFERENCES:
%   [1] Cichocki et al., Nonnegative Matrix and Tensor factorizations, 
%      Wiley, 2009.
%
%
% NOTES:
%   1. Tensor Toolbox required.
%
%
% COPYRIGHT: A.Jukic, ZLAIR, IRB, Zagreb, Croatia
% REVISION:  0.0.1
% DATE:      2012/05/03

% debug mode
DEBUG = 0;

% input
sizeX = size(X);

% default: hosvd without truncation
if(~exist('J','var'))
	J = sizeX;
end


% selected modes
SelMod = find( J~=-1 );
% has to be a row vector because it is used by the for loop (and for loop 
% iterates over rows)
SelMod = SelMod(:)';    

                        
% check sensibility of dimensions of the core
if(any(J(SelMod)>sizeX(SelMod)))
    errtext = ['my_HOSVD: Dimensions in J should be less or equal to' ...
        'dimensions of X.'];
    error(errtext)
end
                        

% prepare
X = tensor(X);
U = cell(length(J),1);

% svd in each mode
for n = SelMod

    if(DEBUG==1)
        fprintf('HOSVD, mode = %d\n',n)
    end
    
	% SVD in mode n
    U{n} = nvecs(X,n,J(n));
    
end

% calculate core
% multiply with inverse (transp) only in modes where factor have been
% calculated - factors in modes with J(n)=-1 are not used
G = ttm(X, U, SelMod, 't');

