%% Generate synthetic 3-order tensor
clear; close all; clc;
Nway = [10,15,16]; % dimension of tensor
coreNway = [5,3,5];

% randomly generate core tensor
G = tensor(randn(coreNway));
A = cell(1,ndims(G));
% randomly generate factor matrices
for i = 1:ndims(G)
    A{i} = randn(Nway(i),coreNway(i));
end
% generate tensor
Mt = full(ttensor(G,A));
M = Mt + 1e-4*randn(Nway);
N = ndims(M);

sr = 0.3; % take 30% samples
p = round(sr*prod(Nway));
known = randsample(prod(Nway),p); 
data = M(known);

opts = [];
opts.maxit = 500; opts.tol = 1e-5;
opts.rank_adj = ones(1,N); % use rank-increasing strategy
opts.rank_max = 1*coreNway;

% initial estimate of multilinear rank
EstCoreNway = [1 1 1];

%% ALSaS
t0 = tic;
[A,C,Out] = ALSaS(data,known,Nway,EstCoreNway,opts);
time = toc(t0);

Mrec = full(ttensor(C,A));

% Reporting
relerr = norm(Mrec(:)-Mt(:))/norm(Mt(:));
fprintf('ALSaS: time = %4.2e, relerr = %4.2e\n\n',time, relerr);

%% iHOOI
t0 = tic;
[A,C,Out] = iHOOI(data,known,Nway,EstCoreNway,opts);
time = toc(t0);

Mrec = full(ttensor(C,A));

% Reporting
relerr = norm(Mrec(:)-M(:))/norm(M(:));
fprintf('iHOOI: time = %4.2e, relerr = %4.2e\n\n',time, relerr);