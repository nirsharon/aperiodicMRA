% script name: "vary_sigma_comparison"
%
% The script compares between the following methods for 3 different
% scenarios. The methods are:
% 1 -- LS
% 2 -- (adapted) EM
% 3 -- Spectral method (over Fourier domain)
%
% NS, Nov. 2017
% ------------------------------------

clear;
% main parameters
repeat_trial = 20;   % how many trials to average for each scenario
L = 25;             % signal length
tic

save_it = 1;
parallel_pool = 1;  % accelerating EM
is_scenario1  = 0;  % fixed sigma, changing M
is_scenario2  = 1;  % changing sigma, fixed M

if parallel_pool
    parallel_nodes = 2;
    if isempty(gcp('nocreate'))
        parpool(parallel_nodes, 'IdleTimeout', 240);
    end
end

% ====================================
% scenario 1 - fixed sigma, changing M
if is_scenario1
    tic
    sigma = 0.5;
    N = sort([floor(logspace(2,5,5)), 5000, 10000]); %10*floor(logspace(2,5,6));
    numS = numel(N);
    
    err_grd = zeros(repeat_trial,numS);
    % err_Jen = zeros(repeat_trial,numS);
    err_EM = zeros(repeat_trial,numS);
    err_SDP = zeros(repeat_trial,numS);
    err_spec = zeros(repeat_trial,numS);
    
    for trial_num=1:repeat_trial
        
        % make data
        x = randn(L, 1);
        xfft = fft(x); xfft = xfft./abs(xfft);
        xfft(1) = 1; xfft(2) = 1;  xfft(end) = 1;
        x = real(ifft(xfft));
        
        rho = rand(L,1); rho = rho/sum(rho);
        
        % runs, different m
        for m = 1:numS
            X = generate_observations(x, N(m), sigma, rho );
            
            % gradient
            [x_est, ~] = MRA_LS(X, sigma);
            x_est = align_to_reference(x_est,x);
            err_grd(trial_num,m) = norm(x_est-x)/norm(x);
            
            % Jennrich
            % [x_est_jen, ~] = spectral_jennrich_like(X, sigma);
            % x_est_jen = align_to_reference(x_est_jen,x);
            % err_Jen(trial_num,m) = norm(x_est_jen-x)/norm(x);
            
            % SPD
            [x_est_sdp, ~] = SDP_solver(X, sigma);
            x_est_sdp = align_to_reference(x_est_sdp, x);
            err_SDP(trial_num,m) = norm(x_est_sdp-x)/norm(x);
            
            % EM
            x_est_em = MRA_EM_NU(X, sigma);
            x_est_em = align_to_reference(x_est_em,x);
            err_EM(trial_num,m) = norm(x_est_em-x)/norm(x);
            
            % direct spectral
            [x_est, ~] = spectral_method(X, sigma);
            x_est = align_to_reference(x_est,x);
            err_spec(trial_num,m) = norm(x_est-x)/norm(x);
        end
        
    end
    % summary
    if repeat_trial>1
        err_grd = mean(err_grd);
        % err_Jen = mean(err_Jen);
        err_EM  = mean(err_EM);
        err_SDP = mean(err_SDP);
        err_spec  = mean(err_spec);
    end
    
    raw_err_SDP = err_SDP;
    raw_err_spec = err_spec;
    err_SDP(err_SDP>1.5)=inf;
    err_spec(err_spec>1.5)=inf;
    
    % ploting the error
    figure;
    loglog(N, err_grd,'LineWidth',3); hold on;
    loglog(N, err_EM,'r','LineWidth',3.6);
   % loglog(N, err_SDP,'-.','LineWidth',3.3);
    loglog(N(raw_err_SDP<1.1), raw_err_SDP(raw_err_SDP<1.1),'-.','LineWidth',3.3);
    loglog(N, err_spec,'--k','LineWidth',3.5);
    legend('LS','EM','SDP','Spectral','Location','SouthWest');
    
    xlabel('N');
    ylabel('Relative error');
    xlim([min(N),max(N)])
    xticks([1000 10000 100000])
    % totalNum = [err_grd(:); err_EM(:); err_SDP(:); err_spec(:)];
    % ylim([min(totalNum) inf])
    set(gca,'FontSize',22)
    
    if save_it
        folder_name = ['NU_ComparisonTest_',date];
        name_it = 'comparison_scenario1';
        mkdir(folder_name);
        cd(folder_name);
        saveas(gcf, name_it, 'fig');
        print('-depsc2', name_it);
        save('data_scenario1');
        cd '../'
    end
end

% ====================================
% scenario 2 - changing sigma, fixed M
if is_scenario2
    sigma = logspace(-3,0,6)*2;
    N = 10000;
    numS = numel(sigma);
    
    err_grd = zeros(repeat_trial,numS);
    err_EM = zeros(repeat_trial,numS);
    err_spec = zeros(repeat_trial,numS);
    
    for trial_num=1:repeat_trial
        
        % make data
        x = randn(L, 1); x = x/norm(x);
%         xfft = fft(x); xfft = xfft./abs(xfft);
%         xfft(1) = 1; xfft(2) = 1;  xfft(end) = 1;
%         
%         x = real(ifft(xfft));
        rho = rand(L,1); rho = rho/sum(rho);
        
        % runs, different m
        for m = 1:numS
            X = generate_observations(x, N, sigma(m), rho );
            
            % gradient LS
            [x_est, ~] = MRA_LS(X, sigma(m));
            x_est = align_to_reference(x_est,x);
            err_grd(trial_num,m) = norm(x_est-x)/norm(x);
            
            % EM
            x_est_em = MRA_EM_NU(X, sigma(m));
            x_est_em = align_to_reference(x_est_em,x);
            err_EM(trial_num,m) = norm(x_est_em-x)/norm(x);
            
            % direct spectral
            [x_est, ~] = spectral_method(X, sigma(m));
            x_est = align_to_reference(x_est,x);
            err_spec(trial_num,m) = norm(x_est-x)/norm(x);
            
        end
    end
    % summary
    if repeat_trial>1
        err_grd = mean(err_grd);
        err_EM  = mean(err_EM);
        err_spec  = mean(err_spec);
    end

    err_grd(err_grd>1)=1;
    err_spec(err_spec>1)  = 1;
    err_spec(err_spec>1)=1;

    % ploting
    figure;
    loglog(sigma,err_grd,'LineWidth',3); hold on;
    loglog(sigma,err_EM,'r','LineWidth',3.6);
    loglog(sigma,err_spec,'--k','LineWidth',3.5);
    leg1 = legend('LS','Adapted EM','Spectral','Location','SouthEast');
    set(leg1,'Interpreter','latex')
    
    xlabel('\sigma');%,'Interpreter','latex');
    ylabel('Relative error','Interpreter','latex');
    xlim([min(sigma),max(sigma)])
    totalNum = [err_grd(:); err_EM(:); err_spec(:)];
    ylim([min(totalNum) inf])
    set(gca,'FontSize',22)
    
    if save_it
        folder_name = ['All_methods_comparison_',date];
        if exist(folder_name,'dir')==7
            cd(folder_name);
        else
            mkdir(folder_name);
            cd(folder_name);
        end
        name_it = 'compare_methods';
        saveas(gcf, name_it, 'fig');
        print('-depsc2', name_it);
        print('-dpdf', name_it);
        save('data_scenario')
        cd '../'
    end
end

poolobj = gcp('nocreate');
delete(poolobj);

totaltime = toc;
disp(['Total time is: ', int2str(floor(totaltime/60)) ,' mins'])

