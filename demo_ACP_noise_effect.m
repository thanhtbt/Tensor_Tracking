%% This demo illustrates the effect of the noise level on the performance of ACP
clc; clear all; close all; 

run_path;
n_exp       = 2;
Noise_level = [0 1e-2 1e-4 1e-6];
tensor_dim  = [20 20 20 500]; R = 2;
epsilon     = zeros(1,tensor_dim(end));  

PER_ACP     = cell(1,length(Noise_level));

for m = 1 : length(Noise_level)
    fac_noise  = Noise_level(m);
    fprintf('\n Noise factor = %f\n', fac_noise);  
    
    PER_Our = zeros(length(tensor_dim),tensor_dim(end));
    for ii = 1: n_exp
        
        fprintf('\n + Run [%d/%d]: ', ii,n_exp)
        %% Generate True Tensor
        N = length(tensor_dim);
        Factor = cell(1,N);
        for n = 1 : N
            Factor{1,n} = randn(tensor_dim(n),R);
        end
        
        % True Data
        Size_Fixed = tensor_dim(1:end-1);
        Num_Slides = tensor_dim(end);
        X_stream   = cp_online_generator(Size_Fixed,R,Num_Slides,epsilon);
        
        OPTS.Factor_True = X_stream(:,2);
        OPTS.Slide_True  = X_stream(:,1);
        OPTS.lambda      = 0.5;
        
        % Add Missing + Gaussian Noise
        X_cell = cell(tensor_dim(end),1);
        Omega  = cell(tensor_dim(end),1);
        for t  = 1 : tensor_dim(end)
            Omega_t     = rand(tensor_dim(1:end-1));
            Omega_t     = 1 .* (Omega_t <= 1);
            Omega{t,1}  = Omega_t;
            X_cell{t,1} = OPTS.Slide_True{t,1}  + fac_noise * randn(tensor_dim(1:end-1));
            X_cell{t,1} = Omega_t .* X_cell{t,1};
        end
        
        %% Algorithms
        t_start = tic;
        [~,PER] = ACP(X_cell,Omega,R,OPTS);
        toc(t_start);
        PER_Our = PER_Our + PER;
        
    end
    PER_ACP{1,m} = PER_Our/n_exp;
    
end



X11 = PER_ACP{1,1};
X12 = PER_ACP{1,2};
X13 = PER_ACP{1,3};
X14 = PER_ACP{1,4};

fig = figure;
hold on;
plot(X11(end,:),'r','LineWidth',1.5)
plot(X12(end,:),'b','LineWidth',1.5)
plot(X13(end,:),'g','LineWidth',1.5)
plot(X14(end,:),'k','LineWidth',1.5)
set(gca, 'YScale', 'log')
lgd = legend('$\sigma = 0$','$\sigma = 10^{-2}$','$\sigma = 10^{-4}$','$\sigma = 10^{-6}$');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);
xlabel('Time Index - $t$','interpreter','latex','FontName','Times New Roman');
ylabel('RE $(\mathbf{U}_{tr}, \mathbf{U}_{es})$','interpreter','latex','FontName','Times New Roman');
set(gca,'FontSize', 20);
grid on;
set(fig, 'units', 'inches', 'position', [0.5 0.5 8 7]);

