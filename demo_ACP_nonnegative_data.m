%% This demo illustrates the tracking ability of Nonnegative ACP (NACP)

clc;
clear all;
run_path;

n_exp            = 3;                     % Number of independent runs
noise_factor     = 1e-3;                  % Noise factor
sampling         = 0.9;                   % (1-sampling) = missing density 
tensor_dim       = [20 20 20 1000];       % Tensor Dimension
R                = 5;                     % CP rank
time_varying     = 1e-3;                  % Time-varying factor
epsilon          = time_varying*ones(1,tensor_dim(end));
epsilon(500)     = 1;                     % Create an abrupt change

%% Main Program
PER_NACP = zeros(length(tensor_dim),tensor_dim(end));

for ii = 1 : n_exp
    
    fprintf('\n Run [%d/%d]: ', ii,n_exp)
    %% Generate True Tensor
    N = length(tensor_dim);
    Factor = cell(1,N);
            
    % True Nonnegative Tensor Data
    Size_Fixed = tensor_dim(1:end-1);
    Num_Slides = tensor_dim(end);
    X_stream   = cp_ntf_generator_sparse(Size_Fixed,R,Num_Slides,epsilon);
    
    OPTS.Factor_True = X_stream(:,2);
    OPTS.Slide_True  = X_stream(:,1);   
   
    %%Add Missing + Gaussian Noise
    X_cell = cell(tensor_dim(end),1);
    Outliers_cell = cell(tensor_dim(end),1);

    Omega  = cell(tensor_dim(end),1);
    sampling = 1;
    for t  = 1 : tensor_dim(end)
        Omega_t     = rand(tensor_dim(1:end-1));
        Omega_t     = 1 .* (Omega_t <= sampling);
        Omega{t,1}  = Omega_t;
        
        X_cell{t,1} = OPTS.Slide_True{t,1}  + noise_factor * randn(tensor_dim(1:end-1));
        X_cell{t,1} = Omega_t .* X_cell{t,1};
           
    end
    
    %% Algorithm 
    t_start = tic;
    [Factor_Es,PER_Our] = NACP(X_cell,Omega,R,OPTS);
    toc(t_start);
    PER_NACP = PER_NACP + PER_Our;
        
        
end

PER_NACP = PER_NACP/n_exp;


%% PLOT RESULTS
makerSize = 14;
numbMarkers = 50;
LineWidth = 2;
set(0, 'defaultTextInterpreter', 'latex');
color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
gree_o  = [0, 0.5, 0];
black_o = [0.25, 0.25, 0.25];
blue_n  = color(1,:);
oran_n  = color(2,:);
yell_n  = color(3,:);
viol_n  = color(4,:);
green_n  = color(5,:);
lblu_n  = color(6,:);
brow_n  = color(7,:);
lbrow_n = [0.5350    0.580    0.2840];

%% 
fig = figure;
subplot(211);
hold on;
k = 1;
K = tensor_dim(end);
ell = 100;
d1 = semilogy(1:k:K,PER_NACP(1,1:k:K),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
d11 = plot(1:ell:K,PER_NACP(1,1:ell:K),...
 'marker','o','markersize',makerSize,...
   'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,PER_NACP(1,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);



lgd = legend([d12],'\texttt{NACP}','Orientation','horizontal');
lgd.FontSize = 16;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);

ylabel(' RE$ (\mathbf{U}_{tr}, \mathbf{U}_{es})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
yticks([1e-4 1e-1 1e2]);
yticklabels({'10^{-4}','10^{-1}','10^2'});

h=gca;
set(gca, 'YScale', 'log')
set(h,'Xtick',0:ell:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 16);
% axis([0 K 5e-5 1e2]);
grid on;
box on;


% fig2 = figure;
subplot(212);
hold on;

d1 = semilogy(1:k:K,PER_NACP(N,1:k:K),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
d11 = plot(1:ell:K,PER_NACP(N,1:ell:K),...
 'marker','o','markersize',makerSize,...
   'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,PER_NACP(N,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);


xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel(' RE($ \mathcal{X}_{tr},  \mathcal{X}_{es}$)','interpreter','latex','FontSize',13,'FontName','Times New Roman');
yticks([1e-4 1e-2 1e0 1e2]);
yticklabels({'10^{-4}','10^{-2}','10^{0}','10^2'});

h=gca;
set(gca, 'YScale', 'log')
set(h,'Xtick',0:ell:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 16);
% axis([0 K 1e-4 1e2]);
grid on;
box on;




