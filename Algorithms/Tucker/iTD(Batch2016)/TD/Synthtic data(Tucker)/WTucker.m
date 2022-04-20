function [Factor_Es,PER] = WTucker(X,Omega,tucker_rank,OPTS);



X = double(X);
ind = find(Omega);
tensor_dim  = size(X);


Xtrue = OPTS.Xtrue;

PER  = [];
Factor_Es = [];


rest = 1*tucker_rank(1:end);
initial = randn(tensor_dim);
topt_maxiter = 2;
poblano_maxiter = tensor_dim(end);
topt_verbose = true;
[X_re,Result] = ...
    wtucker_m(size(X), X(ind), ind, rest, OPTS, ...
    'Xinit', initial, 'maxiter', topt_maxiter, ...
    'poblanomax', poblano_maxiter, 'verbose', topt_verbose);


er = Xtrue - X_re;
PER.X = norm(er(:)) / norm(Xtrue(:));
PER.U = Result.SEP;

end 