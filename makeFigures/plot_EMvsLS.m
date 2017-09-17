% script name: "plot_EMvsLS"
%
%
% May, 19
% ------------------------------------
close all
clear

% main parameters
repeat_trial = 10;   % how many trials for each scenario
N = 25;             % signal length
parl = 1;
save_it = 1;

% parpool for accelerating EM
if parl
    parallel_nodes = 2;
    if isempty(gcp('nocreate'))
        parpool(parallel_nodes, 'IdleTimeout', 240);
    end
end

% main parameters
tic
M = 2000;
sigma_values = 2.5:.15:5.2; % 1.5:.2:4.5;
s = 3:.25:7;

% parameters grid
[s_mesh, sigma_mesh ] = meshgrid(s,sigma_values);
numS = numel(sigma_mesh);

err_grd = zeros(repeat_trial,numS);
err_EM = zeros(repeat_trial,numS);

for trial_num=1:repeat_trial    
    % make signal
    x = randn(N,1);
    
    % runs
    for t = 1:numS
        % create data
        curr_s = s_mesh(t);
        sigma = sigma_mesh(t);
        rho = exp(-0.5*(1:N).^2/curr_s^2)'; rho = rho/sum(rho);
        X = generate_observations(x, M, sigma, rho );
        % IP
        [x_est, ~] = MRA_LS(X, sigma);
        x_est = align_to_reference(x_est, x);
        err_grd(trial_num,t) = norm(x_est-x)/norm(x);
        % EM
        x_est_em = MRA_EM(X, sigma);
        x_est_em = align_to_reference(x_est_em,x);
        err_EM(trial_num,t) = norm(x_est_em-x)/norm(x);
    end
end
save('err_grd','err_grd'); save('err_EM','err_EM');

% summary
if repeat_trial>1
    err_grd = mean(err_grd);
    err_EM  = mean(err_EM);
end
totaltime = toc;
disp(['Total time is: ', int2str(floor(totaltime/60)) ,' mins'])

diff_Arr = (err_grd-err_EM); diff_Arr = diff_Arr/max(max(abs(diff_Arr)));
diff_Arr = reshape(diff_Arr,size(s_mesh));

% ======== binaric map of the differences ===========
figure;
surf(s_mesh, sigma_mesh ,(diff_Arr<0)*1); view(2)
xlabel('Uniformity (s)');
ylabel('Noise (\sigma)');
colorbar;
xlim([min(s) max(s)]);
ylim([min(sigma_values) max(sigma_values)]);

% ======== green-red map of the differences ===========
figure;
surf(s_mesh, sigma_mesh ,diff_Arr*1); view(2)
xlabel('Uniformity (s)');
ylabel('Noise (\sigma)');
set(gca,'FontSize',22); 
xlim([min(s) max(s)]);
ylim([min(sigma_values) max(sigma_values)]);

a = abs(min(diff_Arr(:))); b = max(diff_Arr(:));
t1 = ceil(256*(a/(a+b)));
redColorMap = [linspace(1, 0, t1).^(.15), zeros(1, (256-t1))];
greenColorMap = [zeros(1, 100), linspace(0, 1, 156).^(1.5)];  %red effect
colorMap = [greenColorMap; redColorMap; zeros(1, 256)]';
colormap(colorMap);
colorbar

% diff_Arr = 5*(err_grd-err_EM);
% diff_Arr = reshape(diff_Arr,size(s_mesh));
% figure;
% surf(s_mesh, sigma_mesh ,diff_Arr*1); view(2)
% %title('Region where IP better than EM');
% xlabel('Uniformity (s)');
% ylabel('Noise (sigma)');
% set(gca,'FontSize',22); colorbar;

if save_it
    folder_name = ['CompareIPandEM_',date];
    mkdir(folder_name);
    cd(folder_name);
    saveas(gcf,'IPvsEM','fig');
    print('-depsc2', 'IPvsEM');
end

figure;
surf(s_mesh, sigma_mesh ,reshape(err_grd-err_EM,size(s_mesh))); 
title('The difference: errorIP-errorEM');
xlabel('Uniformity (s)');
ylabel('Noise (sigma)');
set(gca,'FontSize',22); colorbar;
if save_it
    saveas(gcf,'IPvsEM_diff','fig');        
end

figure;
rho1(:,1) = exp(-0.5*(1:N).^2/s(1)^2)'; rho1(:,1) = rho1(:,1)/sum(rho1(:,1));
rho1(:,2) = exp(-0.5*(1:N).^2/s(floor(length(s)/2))^2)'; rho1(:,2) = rho1(:,2)/sum(rho1(:,2));
rho1(:,3) = exp(-0.5*(1:N).^2/s(end)^2)'; rho1(:,3) = rho1(:,3)/sum(rho1(:,3));

plot(1:N, rho1,'LineWidth',3); 
legend(['Most concentrated (s=',num2str(s(1)),')'],...
    ['Mid-range concentrated (s=',num2str(s(floor(length(s)/2))),')'],['Most uniform (s=',num2str(s(end)),')'])
set(gca,'FontSize',22); 
xlim([1,N]);
if save_it
    saveas(gcf,'different_dists','fig'); 
    print('-depsc2','different_dists');
end

% ====================================
if parl
    poolobj = gcp('nocreate');
    delete(poolobj);
end
if save_it
    save('data')
    cd '../'
end
