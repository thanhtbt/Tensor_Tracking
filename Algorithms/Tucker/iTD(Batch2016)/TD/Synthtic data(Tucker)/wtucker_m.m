
function [X, Result] = tucker_opt_track(s, y, ind, r, OPTS, varargin)
% Optimization-based tensor completion by Tucker factorization
% with missing data. 
%
% USAGE: X = tucker_opt(s, y, ind, varargin)
%
% INPUTS: s - size of the tensor
%
%         y - known entries of the tensor
%
%         ind - indexes of known entries
%
%         r - rank estimates
%
%         optional parameters: xinit - initial approximation of X (default:
%         zeros)
%
%                              tol - if change in X is less than tol, stop
%                              (default: 1e-6)
%
%                              maxiter - maximal number of iterations
%                              (default: 100)
%
%                              verbose - show progress (default: false)
%
% OUTPUTS: X - completed tensor
%
%          r - estimated ranks
%
%          errors - objective function through iterations
%
%          C, A - core tensor and factors (A is a cell array, C is a 
%                 tensor)
%
% M. Filipovic, 22.-23.11.2012.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Result =[];
Utrue  = OPTS.Utrue;
Xtrue  = OPTS.Xtrue;



% Tensor order
l = length(s);
maxiter = 100;

%%{
gd = false;
% Number of gradient descent iterations
maxiter_gd = 200;

% Tolerance for gradient descent
gradtol = 1e-6;

% Line search parameters
% (Initial) step size
mi0 = 1;
maxiter_ls = 50;
beta_ls = 0.5;
alpha_ls = 0.05;
%}



tol = 1e-6;
verbose = false;

% Initialization
X = zeros(s);
% X = randn(s);

% Optional parameters
% Check for optional parameters
if(mod(length(varargin),2) == 1)
    error('Optional arguments come in pairs.');
end
for i = 1:2:length(varargin)
    switch lower(varargin{i})                
        case 'xinit'
            X = varargin{i+1};         
        case 'tol'
            tol = varargin{i+1};
        case 'gd'
            gd = varargin{i+1};
%         case 'mi'
%             mi0 = varargin{i+1};    
        case 'maxiter'
            maxiter = varargin{i+1}; 
        case 'verbose'
            verbose = varargin{i+1};
    end
end

X(ind) = y;


RE  = zeros(1,maxiter); 
SEP = zeros(1,maxiter); 



% Initialize factors
if ~exist('C', 'var')
    if ~exist('A', 'var') 
        warning('off', 'MATLAB:eigs:TooManyRequestedEigsForRealSym')
%         hosvd_start = tic;
        [C, A] = my_HOSVD(X, r);
%         hosvd_time = toc(hosvd_start)
    else
        % Initialize C
    end
    
end % if ~exist C

% Weight matrix
W = false(s);
W(ind) = true;
Xw = zeros(s);
Xw(ind) = y;


% warning('off', 'optim:fminunc:SwitchingMethod')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errors = zeros(1, maxiter);
% fit = zeros(1, maxiter);

iter = 0;
done = false;

if verbose
    %fprintf('Tucker-opt (%d iters)', maxiter);
elseif verbose == 2
    %hh = waitbar(0, 'Tensor completion in progress...');
end

while iter < maxiter && ~done  %% edited by LT Thanh 08/05/2020
    
    iter = iter+1;
    
    if verbose
        if mod(iter-1, 50) == 0
           % fprintf('\n.')
        else
          %  fprintf('.')
        end
    end
    
    Xsave = X;
    
    for i = 1:l
        
   
        % Use poblano toolbox
        %%{       
        ncg_display = 'off';
        ncg_stoptol = 1e-5;
        ncg_RelFuncTol = 1e-5;
        ncg_maxiter = 100;
        fun = @(a) fminunc_fun(a, W, X, C, A, Xw, i);
        
        opt_alg = 'ncg';
        if strcmp(opt_alg, 'ncg')
            ncg_update = 'PR';
            ncg_out = ncg(fun, reshape(A{i}, numel(A{i}), 1), ...
                'Update', ncg_update, 'Display', ncg_display, ...
                'StopTol', ncg_stoptol, 'RelFuncTol', ncg_RelFuncTol, ...
                'MaxIters', ncg_maxiter);
        elseif strcmp(opt_alg, 'lbfgs')
            ncg_out = lbfgs(fun, reshape(A{i}, numel(A{i}), 1), ...
                'Display', ncg_display, ...
                'StopTol', ncg_stoptol, 'RelFuncTol', ncg_RelFuncTol, ...
                'MaxIters', ncg_maxiter);
        end
        
        A{i} = reshape(ncg_out.X, size(A{i}));
        
        %}
        
        
        if gd
            %%{
            temp = ttm(tensor(C), A, -i);
            temp_mat = double(tenmat(temp, i));        

            for znj = 1:maxiter_gd

                X_old = X;%double(ttm(temp, A{i}, i));
                objf_ten = W .* X_old - Xw;

                % Gradient of the objective with respect to i-th factor            
                grad = 2 * double(tenmat(tensor(objf_ten), i)) * temp_mat';

                grad_norm = norm(grad, 'fro')^2;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Line search
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                mi = mi0;
                Ain = A{i} - mi * grad;            
                % f_old = sum((y - X_old(ind)).^2);
                f_old = norm(objf_ten(:))^2;
                X = double(ttm(temp, Ain, i)); 
                f_current = sum((y - X(ind)).^2);

                iter_ls = 0;        
                while (f_current - f_old) > ...
                        -(alpha_ls * mi * grad_norm) && iter_ls < maxiter_ls

                    mi = beta_ls * mi;
                    Ain = A{i} - mi * grad;

                    X = double(ttm(temp, Ain, i));
                    f_current = sum((y - X(ind)).^2);

                    iter_ls = iter_ls + 1;
                end

                Aold = A{i};

                if iter_ls < maxiter_ls                
                    A{i} = Ain;            
                end

%                 if norm(A{i} - Aold, 'fro') < gradtol
%                     break
%                 end
            SEP(1,znj) =  sub_est_per(A{1,1},Utrue,'SEP');    
            end % for znj
            %}
        end % if gd
        
    end
    
  
    % Use poblano toolbox
    %%{        
    ncg_display = 'off';
    ncg_stoptol = 1e-5;
    ncg_RelFuncTol = 1e-5;
    ncg_maxiter = 100;
    fun = @(cc) fminunc_core(cc, W, C, A, Xw);        

    opt_alg = 'ncg';
    if strcmp(opt_alg, 'ncg')
        ncg_update = 'PR';
        ncg_out = ncg(fun, reshape(double(C), numel(double(C)), 1), ...
            'Update', ncg_update, 'Display', ncg_display, ...
            'StopTol', ncg_stoptol, 'RelFuncTol', ncg_RelFuncTol, ...
                'MaxIters', ncg_maxiter);
    elseif strcmp(opt_alg, 'lbfgs')
        ncg_out = lbfgs(fun, reshape(double(C), ...
            numel(double(C)), 1), 'Display', ncg_display, ...
            'StopTol', ncg_stoptol, 'RelFuncTol', ncg_RelFuncTol, ...
                'MaxIters', ncg_maxiter);
    end

    C = reshape(ncg_out.X, size(C));

    X = double(ttensor(tensor(C), A));

    err = sum((y - X(ind)).^2);
    %}
        
    if gd
        %%{
        for znj = 1:maxiter_gd                                

            grad = 2 * double(ttm(tensor(W .* X - Xw), A, 't'));

            grad_norm = norm(grad(:))^2;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Line search
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            mi = mi0;
            Cn = C - mi * grad;

            f_old = f_current;

            X = double(ttensor(tensor(Cn), A));   

            err = sum((y - X(ind)).^2);

            iter_ls = 0;        
            while (err - f_old) > ...
                    -(alpha_ls * mi * grad_norm) && iter_ls < maxiter_ls

                mi = beta_ls * mi;
                Cn = C - mi * grad;

                X = double(ttensor(tensor(Cn), A));

                err = sum((y - X(ind)).^2);

                iter_ls = iter_ls + 1;
            end

            Cold = C;

            if iter_ls < maxiter_ls
                C = Cn;        
            end

%             if norm(C(:) - Cold(:)) < gradtol
%                 break
%             end

        end    

        %}
    end % if gd
    
    errors(iter) = err;
        
    if (norm(X(:) - Xsave(:)) < tol)
        done = true;
    end
    
    if verbose == 2
        waitbar(iter / maxiter, hh)
    end
    
    
   
    ER = X(ind) - Xtrue(ind);
    RE(1,iter)  = norm(ER(:)) / norm(Xtrue(ind));
end % while

if verbose
elseif verbose == 2
    close(hh)
end

errors = errors(1:iter);

Result.SEP   = SEP;
Result.RE    = RE;

Result.Factor = A;
end % return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
% Factors
function [f, grad] = fminunc_fun(a, W, X, C, A, Xw, i)

N = numel(X);

temp = ttm(tensor(C), A, -i);

temp_mat = double(tenmat(temp, i));

Amat = reshape(a, size(A{i}));
Xa = ttm(temp, Amat, i);

objf_ten = W .* double(Xa) - Xw;
f = norm(reshape(objf_ten, N, 1))^2;

grad = 2 * reshape(double(tenmat(tensor(objf_ten), i)) * temp_mat', ...
    size(a));

end% return

% Core
function [f, grad] = fminunc_core(cc, W, C, A, Xw)

X = double(ttensor(tensor(reshape(cc, size(C))), A));

objf_ten = tensor(W .* X - Xw);

f = norm(objf_ten)^2;

grad = 2 * reshape(double(ttm(objf_ten, A, 't')), size(cc));

end %return
%}


