function [xest, dist_est] = SDP_solver(X, sigma, true_x, true_rho)
% solving MRA by SDP:
% We use CVX and convex relaxation to form a solution. 
% We assume the signal is of the following structure:
% 1 -- Unit first Fourier entry (DC component)
% 2 -- Unit second Fourier entry (and thus the last one too, real signal)
% These constraints are part of the limits of the convex relaxation, and
% spare us from using the first moment for scaling the signal.
%
% NS, Sep 17

% basic parameters
[L,N] = size(X);
if ~exist('sigma', 'var') || isempty(sigma)
    sigma = std(sum(X, 1))/sqrt(L);      % sample sigma
end

% forming the second moments
M = (1/N)*X*(X')-sigma^2*eye(L);
% debuggin: small norm(circ(true_x)*diag(true_rho([L,1:(L-1)]))*circ(true_x)'- M)

% the power spectrum (modulus of DFT of x)
P_est = sqrt(max((1/N)*sum(abs(fft(X)).^2,2)-L*sigma^2,0));

% tranform the second moment matrix M
F             = fft(eye(L));
D_normalize_x = diag(P_est.^(-1));
M_hat         = D_normalize_x*(F*M*F')*D_normalize_x;
M_hat = (M_hat+M_hat')/2;

% solving SDP with CVX
warning('off');
cvx_precision high;
if nargin<4
    cvx_begin sdp quiet
else
    cvx_begin sdp
end
variable Z(L,L) hermitian
variable z(L) complex

minimize (norm(M_hat.*conj(Z) - circulant(z), 'fro'));
subject to
Z >= 0;      % positive semidefiniteness constraint
diag(Z) == 1;
Z(1,2) == 1;
z(1) == 1; %/N;
z(2:end) == conj(z(end:-1:2));

cvx_end
warning('on');

% extract the eigenvalue
Z = conj(Z);
dist_est = flipud(ifft(z));
[ev, eg] = eigs(Z,2);

% verify the (somehow artificial) constraints on the signal
xf_est = conj(ev(:,1)*conj(ev(1))/(abs(ev(1))^2));
xf_est(3:(end-1)) = (conj(xf_est((end-1):-1:3))+ (xf_est(3:(end-1))))/2;
xf_est(2) = 1; xf_est(end) = 1;

% the estimation
xest = real(ifft(xf_est.*P_est));

if nargin==4  % debugging mode
    xest = align_to_reference(xest, true_x);
    err = norm(xest - true_x)/norm(true_x);
    distribution_err = norm(dist_est - true_rho)/norm(true_rho);
    fprintf('the error is = %.4f\n',err);
    fprintf('the distribution error is = %.4f\n',err);
    fprintf('the spectral gap is = %.4f\n',eg(1)/eg(2));
end

end
