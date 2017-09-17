function [X, shifts,rho_emp] = generate_observations(x, M, sigma, rho, noisetype)
% Given a signal x of length N, generates a matrix X of size N x M such
% that each column of X is a randomly, circularly shifted version of x with
% iid Gaussian noise of variance sigma^2 added on top. If x is complex,
% the noise is also complex.

    x = x(:);
    N = length(x);
    
    X = zeros(N, M);
    % shifts = randi(N, M, 1); uniform distribution
    rho_vec = cumsum(rho);
    shifts_ind = rand(M, 1);
    shifts = discretize(shifts_ind, [0; rho_vec]);
    rho_emp = hist(shifts,N); rho_emp = (rho_emp/sum(rho_emp))';
  % verifying the distribution 
%   figure; 
%   subplot(211); hist(shifts,N); title('empirical distribution')
%   subplot(212); stem(rho); title('underlying distribution')
%

    for m = 1 : M
        X(:, m) = circshift(x, shifts(m));
    end
    
    if ~exist('noisetype', 'var') || isempty(noisetype)
        noisetype = 'Gaussian';
    end
    
    switch lower(noisetype)
        case 'gaussian'
            if isreal(x)
                X = X + sigma*randn(N, M);
            else
                X = X + sigma*(randn(N, M) + 1i*randn(N, M))/sqrt(2);
            end
        case 'uniform'
            if isreal(x)
                X = X + sigma*(rand(N, M)-.5)*sqrt(12);
            else
                error('Uniform complex noise not supported yet.');
            end
        otherwise
            error('Noise type can be ''Gaussian'' or ''uniform''.');
    end

end
