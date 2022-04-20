clear; clc; close all;

run_path;
n_exp = 1;

tensor_dim  = [20 20 20 500];
tucker_rank = [2 2 2 2];

fac_noise   = 1e-2;
sampling    = 0.7;
mag         = 0;  
epsilon     = mag*ones(1,tensor_dim(end));

lambda      = 0.7;

%% metrices
PER_ATD     = zeros(size(tensor_dim,2)+1,tensor_dim(end));
PER_ALSaS   = zeros(size(tensor_dim,2)+1,tensor_dim(end));
PER_iHOOI   = zeros(size(tensor_dim,2)+1,tensor_dim(end));
PER_WTucker = zeros(size(tensor_dim,2)+1,tensor_dim(end));

fprintf('Processing ... ')
for k = 1 : n_exp
    %% Generate Tensor
    [Xtrue_cell,Factor_True] = tucker_generate_tensor(tensor_dim,tucker_rank,0,epsilon);
    OPTS.Factor_True = Factor_True;
    OPTS.Slide_True  = Xtrue_cell;
    OPTS.lambda      = lambda;
   
    % Missing + Gaussian Noise
    X_cell = cell(tensor_dim(end),1);
    Omega  = cell(tensor_dim(end),1);
    for t  = 1 : tensor_dim(end)
        Omega_t = rand(tensor_dim(1:end-1));
        Omega_t = 1 .* (Omega_t <= sampling);
        Omega{t,1}       = Omega_t;
        X_cell{t,1}      = Xtrue_cell{t,1} + fac_noise*randn(tensor_dim(1:end-1));
        X_cell_Miss{t,1} = Omega_t .* X_cell{t,1};
    end
    %%  Online Algorithm
    X_miss     = zeros(tensor_dim);
    X_true     = zeros(tensor_dim);
    Omega_full = zeros(tensor_dim);

    for ii = 1:tensor_dim(end)
        X_miss(:,:,:,ii)     = X_cell_Miss{ii,1};
        X_true(:,:,:,ii)     = X_cell{ii,1};
        Omega_full(:,:,:,ii) = Omega{ii,1};
    end
    OPTS2.Xtrue = X_true;
    
    
    t_start = tic;
    [~,PER] = ATD(X_cell_Miss,Omega,tucker_rank,OPTS);
    t_end   = toc(t_start);
    PER_ATD =  PER;
    fprintf(' + ATD (Proposed): %f(s) \n',t_end)
    

    
    
  %% Batch Algorithm
    t_start = tic;
    [~,PER] = ALSaS_Batch(X_cell_Miss,tucker_rank,OPTS);
    t_end   = toc(t_start);
    fprintf(' + ALSaS: %f(s) \n',t_end)    
    PER_ALSaS = PER_ALSaS + PER;
    
    t_start = tic;
    [~,PER] = iHOOI_Batch(X_cell_Miss,tucker_rank,OPTS);
    t_end   = toc(t_start);
    fprintf(' + iHOOI: %f(s) \n',t_end)    
    PER_iHOOI = PER_iHOOI + PER;
    
    t_start = tic;
    [~,PER] = WTucker_Batch(X_cell,tucker_rank,OPTS);
    t_end = toc(t_start);
    PER_WTucker = PER_WTucker + PER;
    fprintf(' + WTucker: %f(s) \n',t_end)  
     
end


PER_ATD = PER_ATD/n_exp;
PER_WTucker = PER_WTucker/n_exp;
PER_ALSaS = PER_ALSaS/n_exp;
PER_iHOOI = PER_iHOOI/n_exp;

% PLOT RESULTS
makerSize = 14;
numbMarkers = 50;
LineWidth = 2;
set(0, 'defaultTextInterpreter', 'latex');
color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
gree_o  = 'g'; %[0, 0.5, 0];
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
k   = 5;
K   = tensor_dim(end);
KL  = 50;
subplot(311);
hold on;

d2 = semilogy(1:k:K,PER_iHOOI(4,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,PER_iHOOI(4,1:100:end),...
 'marker','o','markersize',makerSize,...
   'linestyle','none','color',blue_o,'LineWidth',2);
d22 = semilogy(1:1,PER_iHOOI(4,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',2);


d3 = semilogy(1:k:K,PER_ALSaS(4,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
d31 = plot(1:100:K,PER_ALSaS(4,1:100:end),...
 'marker','s','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',2);
d32 = semilogy(1:1,PER_ALSaS(4,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',2);

d4 = semilogy(1:k:K,PER_WTucker(4,1:k:end),...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);
d41 = plot(1:100:K,PER_WTucker(4,1:100:end),...
 'marker','h','markersize',makerSize,...
   'linestyle','none','color',viol_n,'LineWidth',2);
d42 = semilogy(1:1,PER_WTucker(4,1:1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color',viol_n,'LineWidth',2);



d1 = semilogy(1:k:K,PER_ATD(4,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d11 = plot(1:100:K,PER_ATD(4,1:100:end),...
 'marker','d','markersize',makerSize,...
   'linestyle','none','color',red_o,'LineWidth',2);
d12 = semilogy(1:1,PER_ATD(4,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',2);


lgd = legend([d22 d32 d42 d12],'\texttt{iHOOI}','\texttt{ALSaS}','\texttt{WTucker}','\texttt{ATD(Proposed)}');
lgd.FontSize = 16;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);

ylabel('RE $(\mathcal{X}_{tr}, \mathcal{X}_{es})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

% 
h=gca;
set(gca, 'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:KL:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 16);
grid on;
box on;




subplot(312);
hold on;

d2 = semilogy(1:k:K,PER_iHOOI(1,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,PER_iHOOI(1,1:100:end),...
 'marker','o','markersize',makerSize,...
   'linestyle','none','color',blue_o,'LineWidth',2);
d22 = semilogy(1:1,PER_iHOOI(1,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',2);


d3 = semilogy(1:k:K,PER_ALSaS(1,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
d31 = plot(1:100:K,PER_ALSaS(1,1:100:end),...
 'marker','s','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',2);
d32 = semilogy(1:1,PER_ALSaS(1,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',2);

d4 = semilogy(1:k:K,PER_WTucker(1,1:k:end),...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);
d41 = plot(1:100:K,PER_WTucker(1,1:100:end),...
 'marker','h','markersize',makerSize,...
   'linestyle','none','color',viol_n,'LineWidth',2);
d42 = semilogy(1:1,PER_WTucker(1,1:1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color',viol_n,'LineWidth',2);


d1 = semilogy(1:k:K,PER_ATD(1,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d11 = plot(1:100:K,PER_ATD(1,1:100:end),...
 'marker','d','markersize',makerSize,...
   'linestyle','none','color',red_o,'LineWidth',2);
d12 = semilogy(1:1,PER_ATD(1,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',2);

ylabel('SEP $(\mathbf{U}_{tr}, \mathbf{U}_{es})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

% 
h=gca;
set(gca, 'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:KL:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 16);
grid on;
box on;


subplot(313)
hold on;

d2 = semilogy(1:k:K,PER_iHOOI(5,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,PER_iHOOI(5,1:100:end),...
 'marker','o','markersize',makerSize,...
   'linestyle','none','color',blue_o,'LineWidth',2);
d22 = semilogy(1:1,PER_iHOOI(5,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',2);


d3 = semilogy(1:k:K,PER_ALSaS(5,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
d31 = plot(1:100:K,PER_ALSaS(5,1:100:end),...
 'marker','s','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',2);
d32 = semilogy(1:1,PER_ALSaS(5,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',2);

d4 = semilogy(1:k:K,PER_WTucker(5,1:k:end),...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);
d41 = plot(1:100:K,PER_WTucker(5,1:100:end),...
 'marker','h','markersize',makerSize,...
   'linestyle','none','color',viol_n,'LineWidth',2);
d42 = semilogy(1:1,PER_WTucker(5,1:1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color',viol_n,'LineWidth',2);



d1 = semilogy(1:k:K,PER_ATD(5,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
d11 = plot(1:100:K,PER_ATD(5,1:100:end),...
 'marker','d','markersize',makerSize,...
   'linestyle','none','color',red_o,'LineWidth',2);
d12 = semilogy(1:1,PER_ATD(5,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',2);

xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('Time(s)','interpreter','latex','FontSize',13,'FontName','Times New Roman');

h=gca;
set(gca, 'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:KL:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 16);
grid on;
box on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);

