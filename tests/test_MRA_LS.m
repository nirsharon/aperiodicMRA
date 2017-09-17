% script name: "test_MRA_LS"

clear;
L = 25;

x = rand(L, 1);  
rho = rand(L,1); rho = rho/sum(rho);

% Generate N observations with noise variance sigma^2
N = 1000;
sigma = .2;
X = generate_observations(x, N, sigma, rho);

% % check moments accuracy
% circ = @(v) toeplitz([v(1); v(end:-1:2)], v)';
% norm(mean(X,2) - cconv(x,rho,L))
% norm(((1/N)*X*(X')-sigma^2*eye(L))-circ(x)*diag(rho)*circ(x)')
% 

% [x_est2, est_dist2] = InteriorPt_test_vestion(X, sigma);
% x_est2 = align_to_reference(x_est2, x);
% norm(x_est2-x)/norm(x)

% MRA_LS(X, sigma, 0, x, rho);
[x_est, ~] = MRA_LS(X, sigma);
x_est = align_to_reference(x_est, x);
err_ls = norm(x_est-x)/norm(x);

%est_dist = align_to_reference(est_dist, rho);
figure; plot(1:L,x); hold on; plot(1:L,x_est); 
title('LS'); legend('original','estimated');

% %=======================
% [x_est, est_dist] = InteriorPt_test_vestion(X, sigma);
% x_est = align_to_reference(x_est, x);
% norm(x_est-x)
% est_dist = align_to_reference(est_dist, rho);
% figure; plot(1:L,x); hold on; plot(1:L,x_est); legend('original','estimated');
% title('Fixed coefficient'); legend('original','estimated');

%=======================
[x_est, est_dist] = spectral_method(X, sigma); %TrustRegion(X,sigma, x, rho);% 
x_est = align_to_reference(x_est, x);
err_spec = norm(x_est-x)/norm(x);
est_dist = align_to_reference(est_dist, rho);
figure; plot(1:L,x); hold on; plot(1:L,x_est); legend('original','estimated');
title('Spectral'); legend('original','estimated');

disp(['The LS error is: ',num2str(err_ls)]);
disp(['The spectral error is: ',num2str(err_spec)]);
