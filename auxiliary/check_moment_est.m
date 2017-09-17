
% check_moment_est
clear;
L = 15;
x = rand(L, 1);
rho = ones(L,1); %rand(L,1); 
rho = rho/sum(rho);

N = 1000:1000:40000;
sigma = 1.5;

% check moments accuracy
circ = @(v) toeplitz([v(1); v(end:-1:2)], v)';
for j=1:numel(N)
    X = generate_observations(x, N(j), sigma, rho);
    m1_est(j) = norm(mean(X,2) - cconv(rho,x,L));
    m2_est(j) = norm(((1/N(j))*X*(X')-sigma^2*eye(L))-circ(x)*diag(rho)*circ(x)','fro');
end

figure; plot(N',[m1_est;m2_est]'); legend('\mu','M')
figure; plot(N', (m1_est./m2_est).^2); hold on; plot(N',(1/(L+L*sigma^2))*ones(numel(N),1))
