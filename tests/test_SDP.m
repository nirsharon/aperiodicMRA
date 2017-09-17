% script name: "test_SDP"

L = 11;

x0 = randn(L, 1);   xfft = fft(x0); xfft = xfft./abs(xfft); 
xfft(1) = 1; xfft(2) = 1;  xfft(end) = 1;
x = real(ifft(xfft));
rho = rand(L,1); rho = rho/sum(rho);

% Generate N observations with noise variance sigma^2
N = 5000;
sigma = 0.5; %.05;
X = generate_observations(x, N, sigma, rho);

% SDP_solver(X, sigma, x, rho)
[x_est, est_dist] = SDP_solver(X, sigma);
x_est = align_to_reference(x_est, x);
norm(x_est-x)
est_dist = align_to_reference(est_dist, rho);
figure; plot(1:L,x); hold on; plot(1:L,x_est); 
title('SDP'); legend('original','estimated');

%=======================
[x_est, est_dist] = spectral_method(X, sigma);
x_est = align_to_reference(x_est, x);
norm(x_est-x)
est_dist = align_to_reference(est_dist, rho);
figure; plot(1:L,x); hold on; plot(1:L,x_est); legend('original','estimated');
title('Direct spectral'); legend('original','estimated');
