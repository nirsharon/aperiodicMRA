% script name: "test_spectral_method"

clear;
L = 25;

x = rand(L,1);
rho = rand(L,1); rho = rho/sum(rho);

% Generate N observations with noise variance sigma^2
N = 1000;
sigma = 0.005;
X = generate_observations(x, N, sigma, rho);

% older version: direct_spectral_v2(X, sigma, x, rho);
[x_est, est_dist] = spectral_method(X, sigma);

x_est = align_to_reference(x_est, x);
norm(x_est-x)
est_dist = align_to_reference(est_dist, rho);

figure; plot(1:L,x); hold on; plot(1:L,x_est); legend('original','estimated');
figure; plot(1:L,rho); hold on; plot(1:L,est_dist); legend('original','estimated');

