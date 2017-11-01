% script name: "test_MRA_EM_NU"

clear;
L = 25;
x = rand(L, 1);

% Generate N observations with noise variance sigma^2
N = 2000;
sigma = norm(x); %2.5;
base_noise = randn(L,N)*sigma;
s_values = 3:.5:8;
number_of_trials = 10;
parallel_pool = 0;

% initializing
Err_EM    = zeros(number_of_trials,numel(s_values));
Err_EM_NU = zeros(number_of_trials,numel(s_values));
Err_LS    = zeros(number_of_trials,numel(s_values));
Err_Spec  = zeros(number_of_trials,numel(s_values));

% accelerate EM
if parallel_pool
    parallel_nodes = 2;
    if isempty(gcp('nocreate'))
        parpool(parallel_nodes, 'IdleTimeout', 240);
    end
end

%main loop
tic

for ind = 1:numel(s_values)
    curr_s = s_values(ind);
    rho = exp(-0.5*(1:L).^2/curr_s^2)'; rho = rho/sum(rho); rho = rho(:);
    
    X_new = generate_observations(x, N, 0, rho);
    X = X_new + base_noise;
    
    for trialnum=1:number_of_trials
        disp(['Trial no. ' num2str(trialnum),' out of ',num2str(number_of_trials)])
        x_est_em = MRA_EM(X, sigma);
        x_est_em = align_to_reference(x_est_em, x);
        Err_EM(trialnum,ind) = norm(x_est_em-x)/norm(x);
        
        [x_est_NU, dist_est_NU] = MRA_EM_NU(X, sigma);
        x_est_NU = align_to_reference(x_est_NU, x);
        Err_EM_NU(trialnum,ind) = norm(x_est_NU-x)/norm(x);
        
        [x_est_ls, ls_est_dist] = MRA_LS(X, sigma);
        x_est_ls = align_to_reference(x_est_ls, x);
        Err_LS(trialnum,ind) = norm(x_est_ls-x)/norm(x);
        
        warning('off');
        [x_est, est_dist] = spectral_method(X, sigma);
        warning('on');
        x_est = align_to_reference(x_est, x);
        Err_Spec(trialnum,ind) = norm(x_est-x)/norm(x);
    end
end
htime = toc;
fprintf('\t\t It took about: %d mins \n', floor(htime/60));

if number_of_trials>1
    Err_EM = mean(Err_EM);
    Err_EM_NU = mean(Err_EM_NU);
    Err_LS = mean(Err_LS);
    Err_Spec = mean(Err_Spec);
end

Err_Spec(Err_Spec>1) = inf;

figure; hold on;
plot(s_values,Err_EM,'LineWidth',3);
plot(s_values,Err_EM_NU,'LineWidth',3);
plot(s_values,Err_LS,'LineWidth',3);

legend('EM','EM NU','LS');%,'Spec');
xlabel('S');
ylabel('Relative error');

poolobj = gcp('nocreate');
delete(poolobj);
