function [A,C,Out] = ALSaS(data,known,Nway,coreNway,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- alternating least squares method for HOSVD with missing values --%
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
if isfield(opts,'tol')        tol = opts.tol;            else tol   = 1e-4;    end
if isfield(opts,'maxT')       maxT = opts.maxT;          else maxT  = inf;    end
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

M = tensor(zeros(Nway));
M(known) = data;
nrmb = norm(data); obj0 = 0.5*nrmb^2;

estMnrm = sqrt(nrmb^2*(prod(Nway)/length(known)));

Nnum = Nway.*coreNway;
totalNum = prod(coreNway)+sum(Nnum);

pinvA = cell(1,N);

if isfield(opts,'A0')
    A0 = opts.A0;
    for n = 1:N
        [Q,R] = qr(A0{n},0);
        pinvA{n} = R\Q';
    end
else
    A0 = cell(1,N);
    for n = 1:N
        % randomly generate each factor matrix
        A0{n} = randn(Nway(n),coreNway(n));
        A0{n} = A0{n}/norm(A0{n},'fro')*estMnrm^(Nnum(n)/totalNum);
        [Q,R] = qr(A0{n},0);
        pinvA{n} = R\Q';
    end
end

if isfield(opts,'C0')
    C0 = opts.C0;
else
    % randomly generate core tensor
    C0 = tensor(randn(coreNway),coreNway);
    C0 = tensor(C0/norm(C0)*estMnrm^(prod(coreNway)/totalNum));
end

obj = obj0;

pinvA0 = pinvA;

if isfield(opts,'M0')
    M0 = tensor(opts.M0);
else
    M0 = ttm(C0,A0);
end

res0 = norm(M0(known)-data);

A = A0; C = C0; M = M0;
M(known) = data;

% residual change used to detect how much progress is made
% if reschg < reschg_tol and rank_adj(n) == 1, adaptively increase the rank
reschg = 1; 
reschg_tol = max(1e-2,10*tol); % 1e-2 is a heuristic and be changed to a smaller number

% set of modes whose ranks can be increased
rank_inc_set = find(rank_adj==1 & coreNway < rank_max);
rank_inc_num = length(rank_inc_set);

midN = ceil(N/2);

sol_C = 1;
start_time = tic;
%fprintf('Iteration:     ');

%% main iterations
for k = 1:maxit
    %fprintf('\b\b\b\b\b%5i',k);
    
    % update core tensor C
    if sol_C
        C = ttm(M,pinvA);
    end
    
    sol_C = 1;
    
    % update factor matrices A
    Z = ttm(C,A,[midN+1:N]);
    getZ = true;
    
    for n = 1:N
        if n > midN && getZ
            Z = ttm(C,A,[1:midN]);
            getZ = false;
        end
        
        if n <= midN
            Bn = tenmat(ttm(Z,A,[1:n-1,n+1:midN]),n);
        else 
            if midN == N-1
                Bn = tenmat(Z,n);
            else
                Bn = tenmat(ttm(Z,A,[midN+1:n-1,n+1:N]),n);
            end
        end
        
        A{n} = tenmat(M,n)*Bn';
        Bnsq = Bn.data*Bn.data';
        A{n} = A{n}.data*pinv(Bnsq);
        pinvA{n} = pinv(A{n});
    end
    
    % update M and impute unknown entries
    M = ttm(C,A);
    obj0 = obj;    res = norm(M(known)-data);
    obj = 0.5*res^2;    
    M(known) = data;
    
    % check progress and adaptively increase ranks
    ratio = res/res0;   reschg = abs(1-ratio);
    Out.reschg(k) = reschg;    
    
    if rank_inc_num > 0 && reschg < reschg_tol
        rank_inc_adaptive();
    end
        
    % --- diagnostics, reporting, stopping checks ---
    relerr1 = abs(obj-obj0)/(obj0+1);    relerr2 = res/nrmb;
    
    % reporting
    Out.hist_obj(k) = obj;
    Out.hist_rel(1,k) = relerr1;
    Out.hist_rel(2,k) = relerr2;
    
    % check stopping criterion
    crit = relerr1 < tol;
    if crit; nstall = nstall+1; else nstall = 0; end
    switch td
        case 0
            % not necessary to increase rank to rank_max
            if nstall>=3 || relerr2<tol C = ttm(M,pinvA); break; end
        case 1
            % stop only when rank increased to rank_max
            if rank_inc_num==0 && (nstall>=3 || relerr2<tol)
                C = ttm(M,pinvA);
                break;
            end
    end
    if toc(start_time)>maxT; break; end;
    
    A0 = A; C0 = C; M0 = M; pinvA0 = pinvA;
    res0 = res;
end % end of main iterations
%fprintf('\n');
Out.iter = k;
% normalize A
for n = 1:N
    [Q,R] = qr(A{n},0);
    A{n} = Q;
    C = ttm(C,R,n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rank_inc_adaptive()
        % increase the estimated rank
        [~,adj_mode] = max(rank_max-coreNway);
        
        if coreNway(adj_mode)<rank_max(adj_mode)
            inc_num = min(rank_inc(adj_mode),rank_max(adj_mode)-coreNway(adj_mode));
            [Q,R] = qr(A{adj_mode},0);
            for ii = 1:inc_num
                rdnx = randn(Nway(adj_mode),1);
                rdnx = rdnx-Q*(Q'*rdnx);
                rdnx = rdnx/norm(rdnx);
                Q = [Q,rdnx];
            end
            A{adj_mode} = Q; A0{adj_mode} = Q;
            pinvA{adj_mode} = Q'; pinvA0{adj_mode} = Q';
            coreNway(adj_mode) = coreNway(adj_mode)+inc_num;
            if coreNway(adj_mode) == rank_max(adj_mode)
                rank_inc_num = rank_inc_num-1;
            end
            
            if rank_inc_num == 0
                nstall = 0;
            end    
            
            C = ttm(M,pinvA); C0 = C;
            sol_C = 0;
        end
        if min(coreNway) == 1
            rank_inc_adaptive();
        end
    end
end
