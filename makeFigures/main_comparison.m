% script name: "main_comparison"
%
% The script compares between the following methods for 3 different
% scenarios. The methods are:
% 1 -- SDP
% 2 -- LS fitting based on gradient optimization (interior point)
% 3 -- Spectral method (over Fourier domain)
% 4 -- Jennrich-like spectral method (spectral over Real domain).
% => Currently this method is not included as it is inferior numerically
% 5 -- EM algorithm
%
% NS, Sep. 2017
% ------------------------------------

clear;
% main parameters
repeat_trial = 10;   % how many trials to average for each scenario
L = 25;             % signal length
tic

save_it = 1;
parallel_pool = 0;  % accelerating EM
is_scenario1  = 1;  % fixed sigma, changing M
is_scenario2  = 0;  % changing sigma, fixed M
is_scenario3  = 0;  % non uniformity of distibution varies

if parallel_pool                
    parallel_nodes = 2;
    if isempty(gcp('nocreate'))
        parpool(parallel_nodes, 'IdleTimeout', 240);
    end
end

% ====================================
% scenario 1 - fixed sigma, changing M
if is_scenario1
    sigma = 0.4;
    %logspace(2,5,7)
    M = [500:5000:10600, 25000, 40000, 80000, 100000]; %10101; %floor(logspace(2,5,7)); %
    numS = numel(M);
    
    err_grd = zeros(repeat_trial,numS);
    err_Jen = zeros(repeat_trial,numS);
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
            X = generate_observations(x, M(m), sigma, rho );
            
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
            if norm(x_est_sdp-x)/norm(x)>10
                disp('here')
            end
            err_SDP(trial_num,m) = norm(x_est_sdp-x)/norm(x);
           
            % EM
            x_est_em = MRA_EM(X, sigma);
            x_est_em = align_to_reference(x_est_em,x);
            err_EM(trial_num,m) = norm(x_est_em-x)/norm(x);
                       
            % direct spectral
            [x_est, ~] = spectral_method(X, sigma);
            x_est = align_to_reference(x_est,x);
            if (norm(x_est-x)/norm(x))>1
                disp('here');  % direct_spectral_v2(X, sigma, x, rho);
            end
            err_spec(trial_num,m) = norm(x_est-x)/norm(x);
        end
        
    end
    % summary
    if repeat_trial>1
        err_grd = mean(err_grd);
        err_Jen = mean(err_Jen);
        err_EM  = mean(err_EM);
        err_SDP = mean(err_SDP);
        err_spec  = mean(err_spec);
    end
    
    % ploting the error
    figure; hold on;
%     plot(M,err_grd,'LineWidth',3); hold on;
%     %plot(M,err_Jen,'--r','LineWidth',3.8);
%     plot(M,err_EM,'r','LineWidth',3.6);
%     plot(M,err_SDP,'-.','LineWidth',3.3);
%     plot(M,err_spec,'--k','LineWidth',3.5);
%     
   semilogx(M,err_grd,'LineWidth',3);
   semilogx(M,err_EM,'r','LineWidth',3.6);
   semilogx(M,err_SDP,'-.','LineWidth',3.3);
   semilogx(M,err_spec,'--k','LineWidth',3.5); 
   legend('LS','EM','SDP','Spectral','Location','NorthWest');

    xlabel('N');
    ylabel('Relative error');
    set(gca,'FontSize',24)
    
    if save_it
      folder_name = ['ComparisonTest_',date];
      name_it = 'comparison_scenario1';
      mkdir(folder_name);   
      cd(folder_name);
      saveas(gcf, name_it, 'fig');        
      print('-depsc2', name_it);
      save('data')
      cd '../'  
    end
    
end
% ====================================
% scenario 2 - changing sigma, fixed M
if is_scenario2
    sigma = logspace(-3,0,5);   %0.1:.3:1.5;
    M = 10000;
    numS = numel(sigma);
    
    err_grd = zeros(repeat_trial,numS);
    err_Jen = zeros(repeat_trial,numS);
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
            X = generate_observations(x, M, sigma(m), rho );
                        
            % gradient LS
            [x_est, ~] = MRA_LS(X, sigma(m));
            x_est = align_to_reference(x_est,x);
            err_grd(trial_num,m) = norm(x_est-x)/norm(x);
            
%             % Jennrich
%             [x_est_jen, ~] = spectral_jennrich_like(X, sigma(m));
%             x_est_jen = align_to_reference(x_est_jen,x);
%             err_Jen(trial_num,m) = norm(x_est_jen-x)/norm(x);
            
            % SPD
            [x_est_sdp, ~] = SDP_solver(X, sigma(m));
            x_est_sdp = align_to_reference(x_est_sdp, x);
            err_SDP(trial_num,m) = norm(x_est_sdp-x)/norm(x);

            % EM
            x_est_em = MRA_EM(X, sigma(m));
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
        err_Jen = mean(err_Jen);
        err_EM  = mean(err_EM);
        err_SDP = mean(err_SDP);
        err_spec  = mean(err_spec);
    end
   
    
    % ploting the error
  %  figure; 
%     plot(sigma,err_grd,'LineWidth',3);
%     plot(sigma,err_EM,'r','LineWidth',3.6);
%     plot(sigma,err_SDP,'-.','LineWidth',3.3);
%     plot(sigma,err_spec,'--k','LineWidth',3.5);
%    semilogx(sigma,err_grd,'LineWidth',3); hold on;
%    semilogx(sigma,err_EM,'r','LineWidth',3.6);
%    semilogx(sigma,err_SDP,'-.','LineWidth',3.3);
%    semilogx(sigma,err_spec,'--k','LineWidth',3.5);
   loglog(sigma,err_grd,'LineWidth',3); hold on;
   loglog(sigma,err_EM,'r','LineWidth',3.6);
   loglog(sigma,err_SDP,'-.','LineWidth',3.3);
   loglog(sigma,err_spec,'--k','LineWidth',3.5);
   legend('LS','EM','SDP','Spectral','Location','NorthWest');
    
   %plot(sigma,err_Jen,'--r','LineWidth',3.8);
    %legend('LS','Jennrich-like','EM','SDP','Spectral');
    xlabel('\sigma');
    ylabel('Relative error');
    xlim([min(sigma),max(sigma)])
    totalNum = [err_grd(:); err_EM(:); err_SDP(:); err_spec(:)];
    ylim([min(totalNum) inf])
    set(gca,'FontSize',24)
    
    if save_it
      folder_name = ['comparison_',date];
      if exist(folder_name,'dir')==7
        cd(folder_name);
      else
        mkdir(folder_name);     
        cd(folder_name);
      end
      name_it = 'comparison_scenario2';
      saveas(gcf, name_it, 'fig');        
      print('-depsc2', name_it);
      save('data')
      cd '../'  
    end
end
% ====================================
% scenario 3 - non uniformity varies
if is_scenario3
    sigma = 1.5;
    s = .5:2:8;
    M = 5000;
    numS = numel(s);
    
    err_grd = zeros(repeat_trial,numS);
    err_Jen = zeros(repeat_trial,numS);
    err_EM = zeros(repeat_trial,numS);
    err_SDP = zeros(repeat_trial,numS);
    err_spec = zeros(repeat_trial,numS);
    
    for trial_num=1:repeat_trial
        
        % make data
        x = randn(L, 1);
        %rho = rand(N,1); rho = rho/sum(rho);
        
        % runs, different m
        for m = 1:numS
            % update the distribution
            rho = exp(-0.5*(1:L).^2/s(m)^2)'; rho = rho/sum(rho);
            X = generate_observations(x, M, sigma, rho );
            
            % gradient
            [x_est, ~] = MRA_LS(X, sigma);
            x_est = align_to_reference(x_est,x);
            err_grd(trial_num,m) = norm(x_est-x)/norm(x);
            
%             % Jennrich
%             [x_est_jen, ~] = spectral_jennrich_like(X, sigma);
%             x_est_jen = align_to_reference(x_est_jen,x);
%             err_Jen(trial_num,m) = norm(x_est_jen-x)/norm(x);
            
            % EM
            x_est_em = MRA_EM(X, sigma);
            x_est_em = align_to_reference(x_est_em,x);
            err_EM(trial_num,m) = norm(x_est_em-x)/norm(x);
            
            % SPD
            [x_est, ~] = SDP_solver(X, sigma);
            x_est = align_to_reference(x_est,x);
            err_SDP(trial_num,m) = norm(x_est-x)/norm(x);
            
            % direct spectral
            [x_est, ~] = spectral_method(X, sigma);
            x_est = align_to_reference(x_est,x);
            err_spec(trial_num,m) = norm(x_est-x)/norm(x);
        end
    end
    % summary
    if repeat_trial>1
        err_grd = mean(err_grd);
        err_Jen = mean(err_Jen);
        err_EM  = mean(err_EM);
        err_SDP = mean(err_SDP);
        err_spec  = mean(err_spec);
    end
    
    % ploting the error
    figure; hold on;
        
    plot(s,err_SDP,'-.','LineWidth',3.3);
    plot(s,err_grd,'LineWidth',3);
    %plot(s,err_Jen,'--r','LineWidth',3.8);
    plot(s,err_EM,'r','LineWidth',3.6);
    plot(s,err_spec,'--k','LineWidth',3.5);
   legend('LS','EM','SDP','Spectral','Location','best');

    %legend('LS','Jennrich-like','EM','SDP','Spectral');
    xlabel('Tendency to uniform (s)');
    ylabel('Error');
    xlim([min(s),max(s)])
    set(gca,'FontSize',24)
    
    if save_it
      folder_name = ['comparison_',date];
      if exist(folder_name,'dir')==7
        cd(folder_name);
      else
        mkdir(folder_name);     
        cd(folder_name);
      end
      name_it = 'comparison_scenario2';
      saveas(gcf, name_it, 'fig');        
      print('-depsc2', name_it);
      save('data')
      cd '../'  
    end
    
end

poolobj = gcp('nocreate');
delete(poolobj);

totaltime = toc;
disp(['Total time is: ', int2str(floor(totaltime/60)) ,' mins'])

