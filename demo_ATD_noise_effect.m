%% This demo illustrates the effect of noise on performance of ATD

clear; clc; close all;
run_path;

n_exp        = 3;                   % Number of independent runs
tensor_dim   = [20 20 20 1000];     % Dimension of streaming tensors
tucker_rank  = [2 2 2 2];           % Tucker rank
sampling     = [0.9];               % (1-Sampling) = Missing density
Noise_factor = [0 1e-1 1e-2 1e-3];  % Noise level
time_varying = 0;                   % Time-varying factor
epsilon      = time_varying*ones(1,tensor_dim(4));
epsilon(500) = 1;                    % Create an abrupt change at t = 500

%% Main Program
PER_ATD = cell(1,length(Noise_factor));

for m = 1 : length(Noise_factor)
 
    fac_noise   = Noise_factor(m);
       
    fprintf('\n // Noise level %f : \n', fac_noise)
    [Xtrue_cell,Factor_True] = tucker_generate_tensor(tensor_dim,tucker_rank,0,epsilon);
    OPTS.Factor_True = Factor_True;
    OPTS.Slide_True  = Xtrue_cell;
    OPTS.lambda      = 0.8;
 
    PER_Our = zeros(length(tensor_dim)+1,tensor_dim(4));
    for k = 1 : n_exp
        fprintf('   + Run [%d/%d]: ',k,n_exp)
        % Missing + Gaussian Noise
        X_cell = cell(tensor_dim(4),1);
        Omega  = cell(tensor_dim(4),1);
        for t  = 1 : tensor_dim(4)
            Omega_t          = rand(tensor_dim(1:end-1));
            Omega_t          = 1 .* (Omega_t <= sampling);
            Omega{t,1}       = Omega_t;
            X_cell{t,1}      = Xtrue_cell{t,1} + fac_noise*randn(tensor_dim(1:end-1));
            X_cell_Miss{t,1} = Omega_t .* X_cell{t,1};
        end
        
        t_start = tic;
        [~,PER] = ATD(X_cell_Miss,Omega,tucker_rank,OPTS);
        toc(t_start)
        PER_Our = PER_Our + PER;
        
    end
    
    PER_ATD{1,m} = PER_Our/n_exp;
    clear  PER_Our
end


X1  = PER_ATD{1,1};
X2  = PER_ATD{1,2};
X3  = PER_ATD{1,3};
X4  = PER_ATD{1,4};


%% PLOT RESULTS
makerSize = 14;
numbMarkers = 50;
LineWidth = 2;
set(0, 'defaultTextInterpreter', 'latex');
color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
gree_o = [0, 0.5, 0];
black_o = [0.25, 0.25, 0.25];

blue_n  = color(1,:);
oran_n  = color(2,:);
yell_n  = color(3,:);
viol_n  = color(4,:);
gree_n  = color(5,:);
lblu_n  = color(6,:);
brow_n  = color(7,:);
lbrow_n = [0.5350    0.580    0.2840];

%%
fig = figure;

k = 1;
K = tensor_dim(4);
hold on;

d1 = semilogy(1:k:K,X1(1,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d11 = plot(1:100:K,X1(1,1:100:end),...
    'marker','h','markersize',makerSize,...
    'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,X1(1,1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d2 = semilogy(1:k:K,X2(1,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,X2(1,1:100:end),...
    'marker','o','markersize',makerSize,...
    'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,X2(1,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d3 = semilogy(1:k:K,X3(1,1:k:end),...
    'linestyle','-','color','g','LineWidth',LineWidth);
d31 = plot(1:100:K,X3(1,1:100:end),...
    'marker','s','markersize',makerSize,...
    'linestyle','none','color','g','LineWidth',LineWidth);
d32 = semilogy(1:1,X3(1,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color','g','LineWidth',LineWidth);


d4 = semilogy(1:k:K,X4(1,1:k:end),...
    'linestyle','-','color','k','LineWidth',LineWidth);

d41 = plot(1:100:K,X4(1,1:100:end),...
    'marker','d','markersize',makerSize,...
    'linestyle','none','color','k','LineWidth',LineWidth);
d42 = semilogy(1:1,X4(1,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color','k','LineWidth',LineWidth);



lgd = legend([d12 d22 d32 d42],'$\sigma = 0$','$\sigma = 10^{-1}$','$\sigma = 10^{-2}$','$\sigma = 10^{-3}$');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);

set(gca, 'YScale', 'log')
xlabel('Time Index - $t$','interpreter','latex','FontName','Times New Roman');
ylabel('SEP$(\mathbf{U}_{tr}, \mathbf{U}_{es})$','interpreter','latex','FontName','Times New Roman');
set(gca,'FontSize', 20);
grid on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);





