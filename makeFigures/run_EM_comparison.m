% script name: "run_EM_comparison"
%
% This is a comparison betweem a standard EM, estimating just the signal,
% and an adpated EM,estimating both signal and distribution.
%
% NS, Nov, 17

clear; close all;
L = 25;
x = rand(L, 1); x= x/norm(x);

% Generate N observations with noise variance sigma^2
N = 2000;
sigma = 1; % = norm(x);
base_noise = randn(L,N)*sigma;
s_values = 3:13; 
number_of_trials = 1; % paper run with = 100;
parallel_pool = 0;

% initializing
Err_EM    = zeros(number_of_trials,numel(s_values));
Err_EM_NU = zeros(number_of_trials,numel(s_values));

% accelerate EM
if parallel_pool
    parallel_nodes = 4;
    if isempty(gcp('nocreate'))
        parpool(parallel_nodes, 'IdleTimeout', 240);
    end
end

%main loop
tic
for trialnum=1:number_of_trials
    disp(['Iteration no. ',num2str(trialnum)])
    for ind = 1:numel(s_values)
        curr_s = s_values(ind);
        rho = exp(-0.5*(1:L).^2/curr_s^2)'; rho = rho/sum(rho); rho = rho(:);
        
        X_new = generate_observations(x, N, 0, rho);
        X = X_new + base_noise;
        
        x_est_em = MRA_EM(X, sigma);
        x_est_em = align_to_reference(x_est_em, x);
        Err_EM(trialnum,ind) = norm(x_est_em-x)/norm(x);
        
       [x_est_NU, dist_est_NU] = MRA_EM_NU(X, sigma);
        x_est_NU = align_to_reference(x_est_NU, x);
        Err_EM_NU(trialnum,ind) = norm(x_est_NU-x)/norm(x);
        
    end
end
htime = toc;
fprintf('\t\t It took about: %d mins \n', floor(htime/60));

if number_of_trials>1
    Err_EM = mean(Err_EM);
    Err_EM_NU = mean(Err_EM_NU);
end

figure; hold on;
plot(s_values,Err_EM,'+-b','LineWidth',3);
plot(s_values,Err_EM_NU,'*-r','LineWidth',3);
legleg = legend('EM','adapted EM','location','SouthEast');
set(legleg,'Interpreter','latex');
xlabel('s','Interpreter','latex');
ylabel('Relative error','Interpreter','latex');
xlim([min(s_values),max(s_values)])
set(gca,'FontSize',22)

folder_name = ['EM_comarisons',date];
name_it = 'compare_em';
mkdir(folder_name);
cd(folder_name);
saveas(gcf, name_it, 'fig');
print('-depsc2', name_it);
print(gcf, '-dpdf', name_it);


% saving the second plot
ind1 = 1;
ind2 = floor(numel(s_values)/2);
ind3 = floor(numel(s_values));

rho1 = exp(-0.5*(1:L).^2/floor(s_values(ind1))^2)'; rho1 = rho1/sum(rho1); rho1 = rho1(:);
rho2 = exp(-0.5*(1:L).^2/floor(s_values(ind2))^2)'; rho2 = rho2/sum(rho2); rho2 = rho2(:);
rho3 = exp(-0.5*(1:L).^2/ceil(s_values(ind3))^2)'; rho3 = rho3/sum(rho3); rho3 = rho3(:);

figure; hold on;
plot(1:L,rho3,'LineWidth',3);
plot(1:L,rho2,'LineWidth',3);
plot(1:L,rho1,'LineWidth',3);
xlim([1,L])
set(gca,'FontSize',22)
leg1 = legend(['s=',num2str(ceil(s_values(ind3)))],['s=',num2str(floor(s_values(ind2)))],['s=',num2str(floor(s_values(ind1)))],'location','NorthEast');%,'Interpreter','latex');%,'Spec');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',20);

name_it = 'distributions';
saveas(gcf, name_it, 'fig');
saveas(gcf, name_it, 'png');
print('-depsc2', name_it);

save('data_comparison_EM');
cd '../'

poolobj = gcp('nocreate');
delete(poolobj);
