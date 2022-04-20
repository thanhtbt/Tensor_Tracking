function [A,C,Out] = iHOOI(data,known,Nway,coreNway,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% incomplete Higher-order orthogonality iteration for HOSVD from missing values
% Input: 
%       data:       observed values of the underlying tensor
%       known:      indices of observed values
%       Nway:       size of the underlying tensor
%       coreNway:   initial size of the core tensor (can be adaptively
%                   increased)
%       opts.
%           maxit:  maximum number of iterations (500)
%           tol:    stopping tolerance (1e-4)
%           maxT:   maximum running time (inf)
%           rank_max: maximal rank estimates
%           rank_adj: whether adaptively increase rank (zeros(1,N): no)
%           rank_inc: how many to increase rank each time (ones(1,N))
%           td:     whether increase coreNway to rank_max (0: no)
%
% Output:
%       A: factor matrices
%       C: core tensor
%       Out.
%           iter: number of actual iterations
%           hist_obj: objective values
%           hist_rel: relative errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Yangyang Xu, 11/18/2014
%
% Reference:
% Yangyang Xu. On higher-order singular value decomposition from incomplete
% data. arXiv preprint. arXiv:1411.4324

N = length(Nway); % number of modes / ways of the underlying tensor
%% Parameters and defaults
if isfield(opts,'maxit')      maxit = opts.maxit;        else maxit = 500;   end
if isfield(opts,'tol')        tol = opts.tol;            else tol = 1e-4;    end
if isfield(opts,'maxT')       maxT = opts.maxT;          else maxT = inf;    end
if isfield(opts,'rank_inc') rank_inc = opts.rank_inc; else rank_inc = ones(1,N);    end
if isfield(opts,'rank_max') rank_max = opts.rank_max; else rank_max = 50*ones(1,N); end
% td:
%   0: coreNway may not increase to rank_max
%   1: coreNway finally equals rank_max
if isfield(opts,'td')         td = opts.td;              else td = 0;        end

% rank_adj:
%    0: fix rank
%    1: increase rank
if isfield(opts,'rank_adj') rank_adj = opts.rank_adj; else rank_adj = zeros(1,N);   end

%% Data preprocessing and initialization
[known,id] = sort(known); data = data(id);

nrmb = norm(data); obj0 = 0.5*nrmb^2;

if isfield(opts,'A0')
    A0 = opts.A0;
else
    A0 = cell(1,N);
    for n = 1:N
        % randomly generate each factor matrix
        [A0{n},~] = qr(randn(Nway(n),coreNway(n)),0);
    end
end

obj = obj0;

A = A0;
if (isfield(opts,'A0') && isfield(opts,'C0')) || isfield(opts,'M0')
    if isfield(opts,'M0')
        X = tensor(opts.M0);
        X(known) = data;
    else
        X = ttm(opts.C0,opts.A0);
        X(known) = data;
    end
else
    X = tensor(zeros(Nway),Nway); 
    X(known) = data;
    X = ttm(ttm(X,A,'t'),A);
end
Y = ttm(ttm(X,A,'t'),A);

% residual change used to detect how much progress is made
% if reschg < reschg_tol and rank_adj(n) == 1, adaptively increase the rank
reschg = 1; 
reschg_tol = max(1e-2,10*tol); % 1e-2 is a heuristic and be changed to a smaller numberrank_inc_set = find(rank_adj==1 & coreNway < rank_max);

% set of modes whose ranks can be increased
rank_inc_set = find(rank_adj==1 & coreNway < rank_max);
rank_inc_num = length(rank_inc_set);

midN = ceil(N/2);

start_time = tic;

%fprintf('Iteration:     ');

%% main iterations
for k = 1:maxit
    %fprintf('\b\b\b\b\b%5i',k);
    
    % update X
    X = Y;
    X(known) = data;
    
    % update factor matrices
    Z = ttm(X,A,[midN+1:N],'t');
    getZ = true;
    for n = 1:N
        if n > midN && getZ
            Z = ttm(X,A,[1:midN],'t');
            getZ = false;
        end
        
        if n <= midN
            Bn = tenmat(ttm(Z,A,[1:n-1,n+1:midN],'t'),n);
        else 
            if midN == N-1
                Bn = tenmat(Z,n);
            else
                Bn = tenmat(ttm(Z,A,[midN+1:n-1,n+1:N],'t'),n);
            end
        end
        
        Bn = Bn.data;
        
        if size(Bn,1) <= size(Bn,2) || coreNway(n) > size(Bn,2)
            [U,S] = eig(Bn*Bn');
        else
            [V,S] = eig(Bn'*Bn);
            V = V(:,end-coreNway(n)+1:end);
            s = diag(S);
            s = s(end-coreNway(n)+1:end);
            U = Bn*V*diag(s.^(-0.5));
        end

        A{n} = U(:,end-coreNway(n)+1:end);
    end
    
    
    
    % --- diagnostics, reporting, stopping checks ---
    Y = ttm(ttm(X,A,'t'),A);
    obj = .5*norm(Y-X)^2;
    ratio = sqrt(2*obj)/sqrt(2*obj0);   reschg = abs(1-ratio);
    Out.reschg(k) = reschg;
    Out.obj(k) = obj;
    res = abs(obj-obj0)/(obj0+1);
    
    if rank_inc_num > 0 && reschg < reschg_tol
        rank_inc_adaptive();
    end
    
    crit = res < tol;
    if crit; nstall = nstall+1; else nstall = 0; end
    
    switch td
        case 0
            % not necessary to increase rank to rank_max
            if nstall>=3 || sqrt(2*obj)/nrmb < tol 
                break; 
            end
        case 1
            % stop only when rank increased to rank_max
            if rank_inc_num==0 && (nstall>=3 || sqrt(2*obj)/nrmb < tol)
                break;
            end
    end   
    
    if toc(start_time)>maxT; break; end;
    
    A0 = A; obj0 = obj;
end % end of main iterations
C = ttm(X,A,'t');
%fprintf('\n'); 
Out.iter = k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rank_inc_adaptive()
        % increase the estimated rank
        [~,adj_mode] = max(rank_max-coreNway);
        
        if coreNway(adj_mode)<rank_max(adj_mode)
            inc_num = min(rank_inc(adj_mode),rank_max(adj_mode)-coreNway(adj_mode));
            for ii = 1:inc_num
                rdnx = randn(Nway(adj_mode),1);
                rdnx = rdnx-A{adj_mode}*(A{adj_mode}'*rdnx);
                rdnx = rdnx/norm(rdnx);
                A{adj_mode} = [A{adj_mode},rdnx];
            end
            A0{adj_mode} = A{adj_mode};
            coreNway(adj_mode) = coreNway(adj_mode)+inc_num;
            if coreNway(adj_mode) == rank_max(adj_mode)
                rank_inc_num = rank_inc_num-1;
            end
            
            if rank_inc_num == 0
                nstall = 0;
            end          
        end        
    end
end