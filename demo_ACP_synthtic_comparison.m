%% 1. This demo illustrates a performance comparsion of the-state-of-the-art
%%    streaming CP decomposition methods
%% 2. The experiment setup for 3-way tensors only (since compared algorithm 
%%    TeCPSDG, OLSTEC, CP-PETRELS are designed for handling 3-way tensors only)

%% 3-way streaming tensor
clc; clear all; run_path;
n_exp      = 2;               % Number of independent runs
tensor_dim = [100 100 1000];  % Dimension of streaming tensors
R          = 5;               % CP Rank
sampling   = 0.9;             % (1 - sampling) = Missing density

fac_noise    = 1e-3;          % noise factor
time_varying = 1e-3;          % time-varying factor
epsilon      = time_varying*ones(1,tensor_dim(end));
epsilon(600) = 1;             % aimed to create an abrupt change at t = 600

PER_ACP        = zeros(length(tensor_dim),tensor_dim(end));
PER_TeSGD      = zeros(length(tensor_dim),tensor_dim(end));
PER_OLSTEC     = zeros(length(tensor_dim),tensor_dim(end));
PER_CP_PETRELS = zeros(length(tensor_dim),tensor_dim(end));

%% Main Program
for ii = 1 : n_exp
    
    fprintf('\n Run [%d/%d]: \n', ii,n_exp)
    %% Generate True Tensor
    N = length(tensor_dim);
    Factor = cell(1,N);
    for n = 1 : N
        Factor{1,n} = randn(tensor_dim(n),R);
    end
        
    % True Data
    Size_Fixed       = tensor_dim(1:end-1);
    Num_Slides       = tensor_dim(end);
    X_stream         = cp_online_generator(Size_Fixed,R,Num_Slides,epsilon);
    OPTS.Factor_True = X_stream(:,2);
    OPTS.Slide_True  = X_stream(:,1);
    OPTS.lambda      = 0.5;
    
    % Add Missing + Gaussian Noise
    X_cell = cell(tensor_dim(end),1);
    Omega  = cell(tensor_dim(end),1);
    for t  = 1 : tensor_dim(end)
        Omega_t     = rand(tensor_dim(1:end-1));
        Omega_t     = 1 .* (Omega_t <= sampling); 
        Omega{t,1}  = Omega_t; % Obsevation Mask
        X_cell{t,1} = OPTS.Slide_True{t,1}  + fac_noise * randn(tensor_dim(1:end-1));
        X_cell{t,1} = Omega_t .* X_cell{t,1};
    end
    
    %% 
    X_Data  = zeros(tensor_dim);
    Omega_c = zeros(tensor_dim);
    for ii  = 1 : tensor_dim(end)
        X_Data(:,:,ii)  = X_cell{ii,1};
        Omega_c(:,:,ii) = Omega{ii,1};     % Mask for streaming CP algorithms
    end
    
    %% TeSGD  
    maxepochs               = 1;
    verbose                 = 2;
    tolcost                 = 1e-8;
    permute_on              = false;
    clear options;
    options.maxepochs       = maxepochs;
    options.tolcost         = tolcost;
    options.permute_on      = permute_on;
    options.stepsize        = 0.1;
    options.lambda          = 0.001;  % Forgetting paramter
    options.mu              = 0.01;   % Regualization paramter
    options.tw_flag         = 0;      % 0:Exponential Window, 1:Truncated Window (TW)
    options.tw_len          = 10;     % Window length for Truncated Window (TW) algorithm
    options.store_subinfo   = true;
    options.store_matrix    = false;
    options.verbose         = verbose;
   
    t_start = tic;
    [PER,~] = TeCPSGD_M(X_Data,Omega_c,R,options,OPTS);
    t_end   = toc(t_start);
    PER_TeSGD(1,:) = PER_TeSGD(1,:) + PER.A;
    PER_TeSGD(2,:) = PER_TeSGD(2,:) + PER.B;
    PER_TeSGD(3,:) = PER_TeSGD(3,:) + PER.X;
    
    fprintf('+ TeSGD (2015): %f(s)\n',t_end)

    
    %%  OLSTEC
    clear options;
    options.maxepochs       = maxepochs;
    options.tolcost         = tolcost;
    options.permute_on      = permute_on;    
    options.lambda          = 0.7;   % Forgetting paramter
    options.mu              = 0.05;  % Regualization paramter
    options.tw_flag         = 0;     % 0:Exponential Window, 1:Truncated Window (TW)
    options.tw_len          = 10;    % Window length for Truncated Window (TW) algorithm
    options.store_subinfo   = true;     
    options.store_matrix    = false; 
    options.verbose         = verbose; 

    OPTS_OLSTEC = [];
    OPTS_OLSTEC.Factor_True = OPTS.Factor_True;
    OPTS_OLSTEC.Slide_True  = OPTS.Slide_True;
    OPTS_OLSTEC.tensor_dims = tensor_dim;
    Xinit.A = randn(tensor_dim(1), R);
    Xinit.B = randn(tensor_dim(2), R);
    Xinit.C = randn(tensor_dim(3), R);
    OPTS_OLSTEC.Xinit = Xinit;

    t_start   = tic;
    [PER,~,~] = OLSTEC_M(X_Data,Omega_c,R,options,OPTS_OLSTEC);
    t_end     = toc(t_start);
    PER_OLSTEC(1,:) = PER_OLSTEC(1,:)  + PER.A;
    PER_OLSTEC(2,:) = PER_OLSTEC(2,:)  + PER.B;
    PER_OLSTEC(3,:) = PER_OLSTEC(3,:)  + PER.X;

    fprintf('+ OLSTEC (2019): %f(s)\n',t_end)
    
    %% CP-PETRELS
    clear options
    t_start = tic;
    [PER,~] = PETRELS_SOAP(X_Data,Omega_c,R,OPTS);
    t_end = toc(t_start);
    PER_CP_PETRELS(1,:) = PER_CP_PETRELS(1,:) + PER.A;
    PER_CP_PETRELS(2,:) = PER_CP_PETRELS(2,:) + PER.B;
    PER_CP_PETRELS(3,:) = PER_CP_PETRELS(3,:) + PER.X;
    
    fprintf('+ CP_PETRELS (2016): %f(s)\n',t_end)
    
    
    
    %% Our Method: ACP
    t_start = tic;
    [~,PER] = ACP(X_cell,Omega,R,OPTS);
    t_end   = toc(t_start);
    PER_ACP = PER_ACP + PER;
    fprintf('+ Our Method (2022): %f(s)\n',t_end)


end

PER_ACP         = PER_ACP/n_exp;
PER_OLSTEC      = PER_OLSTEC/n_exp;
PER_TeSGD       = PER_TeSGD/n_exp;
PER_CP_PETRELS  = PER_CP_PETRELS/n_exp;

%% PLOT RESULTS
makerSize   = 14;
numbMarkers = 50;
LineWidth   = 2;
set(0, 'defaultTextInterpreter', 'latex');
color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
gree_o  = 'g';%[0, 0.5, 0];
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
k = 5;
K = tensor_dim(end);
KL = 200;
subplot(211);
hold on;

d2 = semilogy(1:k:K,PER_TeSGD(3,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,PER_TeSGD(3,1:100:end),...
 'marker','o','markersize',makerSize,...
   'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,PER_TeSGD(3,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d3 = semilogy(1:k:K,PER_OLSTEC(3,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
d31 = plot(1:100:K,PER_OLSTEC(3,1:100:end),...
 'marker','s','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',LineWidth);
d32 = semilogy(1:1,PER_OLSTEC(3,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);

d4 = semilogy(1:k:K,PER_CP_PETRELS(3,1:k:end),...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);
d41 = plot(1:100:K,PER_CP_PETRELS(3,1:100:end),...
 'marker','h','markersize',makerSize,...
   'linestyle','none','color',viol_n,'LineWidth',LineWidth);
d42 = semilogy(1:1,PER_CP_PETRELS(3,1:1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);



d1 = semilogy(1:k:K,PER_ACP(3,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
d11 = plot(1:100:K,PER_ACP(3,1:100:end),...
 'marker','d','markersize',makerSize,...
   'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,PER_ACP(3,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);



lgd = legend([d22 d32 d42 d12],'\texttt{TeCPSGD}','\texttt{OLSTEC}','\texttt{CP-PETRELS}','\texttt{ACP(Proposed)}');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);
ylabel('RE $(\mathcal{X}_{tr}, \mathcal{X}_{es})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

% 
h=gca;
set(gca, 'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:KL:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 20);
% axis([0 K 5e-3 1e0]);
grid on;
box on;

subplot(212);
hold on;
d2 = semilogy(1:k:K,PER_TeSGD(1,1:k:end),...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);
d21 = plot(1:100:K,PER_TeSGD(1,1:100:end),...
 'marker','o','markersize',makerSize,...
   'linestyle','none','color',blue_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,PER_TeSGD(1,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color',blue_o,'LineWidth',LineWidth);


d3 = semilogy(1:k:K,PER_OLSTEC(1,1:k:end),...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
d31 = plot(1:100:K,PER_OLSTEC(1,1:100:end),...
 'marker','s','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',LineWidth);
d32 = semilogy(1:1,PER_OLSTEC(1,1:1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);

d4 = semilogy(1:k:K,PER_CP_PETRELS(2,1:k:end),...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);
d41 = plot(1:100:K,PER_CP_PETRELS(2,1:100:end),...
 'marker','h','markersize',makerSize,...
   'linestyle','none','color',viol_n,'LineWidth',LineWidth);
d42 = semilogy(1:1,PER_CP_PETRELS(2,1:1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color',viol_n,'LineWidth',LineWidth);


d1 = semilogy(1:k:K,PER_ACP(1,1:k:end),...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

d11 = plot(1:100:K,PER_ACP(1,1:100:end),...
 'marker','d','markersize',makerSize,...
   'linestyle','none','color',red_o,'LineWidth',LineWidth);
d12 = semilogy(1:1,PER_ACP(1,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('RE $(\mathbf{U}_{tr}, \mathbf{U}_{es})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

h=gca;
set(gca, 'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:KL:K,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 20);
% axis([0 K 5e-3 10e0]);
grid on;
box on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);


