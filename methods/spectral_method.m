function [est_x, est_dist] = spectral_method(X, sigma, true_x, true_rho)
% Solving MRA by a direct spectral solver, using the diagonal form of
% 
% \[     (D_{\abs{x}^{-1}F) M2 (D_{\abs{x}^{-1}F)' ,\]
%
% where F is the (unnormalized) Fourier matrix, \abs{x} is the power
% spectrum of $x$, and M2 is our estimation for the second moment.
%
% Last two inputs are merely for debugging
%
% NS, September 17.

[L,N] = size(X);

if ~exist('sigma', 'var') || isempty(sigma)
    sigma = std(sum(X, 1))/sqrt(L);      % sample sigma
end

circ = @(v) toeplitz([v(1); v(end:-1:2)], v)';

% estimating the two (sample) moments
mu = mean(X,2);
M = (1/N)*X*(X')-sigma^2*eye(L);

% the power spectrum (modulus of DFT of x)
P_est = (1/N)*sum(abs(fft(X)).^2,2)-L*sigma^2;
if ~isequal(1*(P_est<1e-10),zeros(L,1))
    warning('Not enough samples for power spectrum positive estimation. Projection applied.');
    P_est(P_est<1e-10) = 1e-10;
end
P_est = P_est.^(.5);

% tranform the second moment matrix M and extract eigenvectors
F             = fft(eye(L));
D_normalize_x = diag(P_est.^(-1));
conj_mat      = (D_normalize_x*F);
M_hat         = (conj_mat*M*conj_mat');
try
 [V, D] = eig(M_hat);
catch
  M_hat = (M_hat+M_hat')/2;  % impose symmetry if needed
  [V, D] = eig(M_hat);
end

% choosing distinct entry if exists 
raw_rho = abs(diag(D));
tol = 1e-4;
[ii,jj,kk] = uniquetol(raw_rho, tol);
f = histc(kk,1:numel(jj));         % frequencies corresponding to ii
idx = f<2;                         % "single" frequencies
idx1 = ismember(raw_rho,ii(idx));  % ind of the "singles"
if (1*idx1)==zeros(size(raw_rho))
    warning('spectral_method: No distinct distribution entry. Result may be inaccurate');
    % choosing eigenvector that leads to least imaginary solution
    Z = fft(V);
    [~, ind] = min(sum(imag(Z).^2));
else 
    % choosing eigenvector that leads to least imaginary solution  
    % ------- other alternative: a cost function like ----------
    % \[ norm(M_hat-circ(xhat)*diag(est_dist)*circ(xhat)') ,\]
    % where xhat = (sign(V(1,j))/sign(mean(X(:))))*V(:,j)
    % ----------------------------------------------------------
    [~, ind_in_idx1] = min(sum(imag(fft(V(:,idx1))).^2));
    val = raw_rho(idx1);    
    ind = find(~(raw_rho-val(ind_in_idx1)));
end

% estimating the signal
xj_hat = sqrt(L)*(P_est.*V(:,ind));
est_x  = real(ifft(xj_hat));
est_x  = est_x*(mean(mu)/mean(est_x)); % scale correction

% estimating the distribution
est_dist = (circ(est_x))\mu;

if nargin==4 % debugging mode
    % if true values exist, output the aligned signal and rho
    est_x = align_to_reference(est_x, true_x);
    est_dist = align_to_reference(est_dist, true_rho);
    figure; plot(est_x); hold on; plot(true_x);
    norm(est_x-true_x)
    figure; plot(est_rho); hold on; plot(true_rho);
end

end
