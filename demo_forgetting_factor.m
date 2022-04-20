%% This demo illustrates the effect of the forgetting factor 

clc; clear all; close all; 
run_path;

n_exp      = 3;                 % Number of independent runs
fac_noise  = 1e-3;              % Noise factor 
sampling   = 0.9;               % (1-sampling) = missing density
tensor_dim = [20 20 500];       % Tensor dimension
R          = 5;                 % CP rank

% 
Lambda = [0.1:0.1:1];           % Forgetting factor
MAG    = [1e-5 1e-4 1e-3 1e-2 1e-1]; % Time-varing factor

N = length(tensor_dim);
T = tensor_dim(end);

PER_U = cell(length(Lambda),length(MAG));
PER_X = cell(length(Lambda),length(MAG));

for ii = 1:length(Lambda)
    for jj = 1 : length(MAG)
        PER_U{ii,jj} = zeros(N-1,T);
        PER_X{ii,jj} = zeros(1,T);
    end
end


for n = 1 : n_exp
for ii = 1 : length(Lambda)

    OPTS.lambda  = Lambda(ii);
    for jj = 1 : length(MAG)
        fprintf('\n // Run %d:   [Lambda = %f  | time-varying =  %f] \n',n, Lambda(ii),MAG(jj))
        fprintf('        ')
        mag = MAG(jj);
        epsilon   = mag*ones(1,tensor_dim(end));
        %% Generate True Tensor

        Size_Fixed = tensor_dim(1:end-1);
        Num_Slides = tensor_dim(end);
        X_stream   = cp_online_generator(Size_Fixed,R,Num_Slides,epsilon);
        OPTS.Factor_True = X_stream(:,2);
        OPTS.Slide_True  = X_stream(:,1);

        %% Add Missing + Gaussian Noise
        X_cell = cell(tensor_dim(end),1);
        Omega  = cell(tensor_dim(end),1);
        for t  = 1 : tensor_dim(end)
            Omega_t = rand(tensor_dim(1:end-1));
            Omega_t = 1 .* (Omega_t <= sampling);
            Omega{t,1}  = Omega_t;
            X_cell{t,1} = OPTS.Slide_True{t,1}  + fac_noise * randn(tensor_dim(1:end-1));
            X_cell{t,1} = Omega_t .* X_cell{t,1};
        end
        t_start = tic;
        [~,PER] = ACP(X_cell,Omega,R,OPTS);
        toc(t_start);
        PER_U{ii,jj} = PER_U{ii,jj} + PER(1:N-1,1:end);
        PER_X{ii,jj} = PER_X{ii,jj} + PER(N,1:end);
    end

end

end


%% 
t_train = 100;

PER_SUM = zeros(length(Lambda),length(MAG));
for jj = 1:length(MAG)
    c = MAG(jj);
    for ii = 1:length(Lambda)
        PER_PLOT = PER_X{ii,jj}/n_exp;
        PER_SUM(ii,jj) = sum(PER_PLOT(t_train:end))/length(PER_PLOT(t_train:end));
    end

end

figure; surf(PER_SUM(:,2:end));
set(gca,'ColorScale','log')
set(gca, 'ZScale', 'log')
caxis([1e-2,1]); %set value range of colorbar

