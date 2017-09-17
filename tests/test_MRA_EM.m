% script name: "test_MRA_EM"

L = 11;
x = rand(L, 1);  
rho = rand(L,1); rho = rho/sum(rho);

% Generate N observations with noise variance sigma^2
N = 5000;
sigma = 0.00001;
X = generate_observations(x, N, sigma, rho);

x_est = MRA_EM(X, sigma);
x_est = align_to_reference(x_est, x);
norm(x_est-x)
figure; plot(1:L,x); hold on; plot(1:L,x_est); 
title('EM'); legend('original','estimated');

%=======================
[x_est, est_dist] = spectral_method(X, sigma);
x_est = align_to_reference(x_est, x);
norm(x_est-x)
figure; plot(1:L,x); hold on; plot(1:L,x_est); legend('original','estimated');
title('Spectral'); legend('original','estimated');

% est_dist = align_to_reference(est_dist, rho);
% norm(est_dist-rho)
% figure; plot(1:L,rho); hold on; plot(1:L,est_dist);
% title('Spectral -- distribution'); legend('original','estimated');
