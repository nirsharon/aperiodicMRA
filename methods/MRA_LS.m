function [x_est, est_dist] = MRA_LS(X, sigma, initial_guess, true_x, true_rho)
% solving MRA with LS fitting. 
% We are using both first and second moments with a weight on the latter 
% which depends on the noise level.
% The "initial_guess" variable is for using the spectral method as initial
% guess. The defualt is a random guess.
%
% NS, Sep 17.

if nargin==5
    dub_mode = 1;
else
    dub_mode = 0;
end

if nargin<3
    initial_guess = 0;
end

% basic parameters
[L,N] = size(X);
if ~exist('sigma', 'var') || isempty(sigma)
    sigma = std(sum(X, 1))/sqrt(L);      % sample sigma
end

% the moments
mu = mean(X,2);
P = (1/N)*X*(X'); %-sigma^2*eye(L); 

% the cost function weight 
%lambda = 1/(L*(1+3*sigma^2));
lambda = 1/(L*(1+sigma^2));

% cost function
fun = @(x) costFunctionV2(x(1:L), x((L+1):2*L), mu, P, lambda, sigma);

% initial guess
if initial_guess
    [x0, est_dist ] = spectral_method(X, sigma);
else
    est_dist = rand(L,1); est_dist = est_dist/sum(est_dist);
    x0 = rand(L,1); x0 = x0*(2*mean(mu)/mean(x0)); %twice the scale
end
x0 = [x0 ;est_dist(:)];

% constraints
A = [zeros(L), -eye(L)];
b = zeros(L,1);
Aeq = [zeros(1,L),ones(1,L)];
beq = 1;

% solver options
options = optimoptions(@fmincon,'Algorithm','interior-point');
options.MaxIterations          = 2000;
options.MaxFunctionEvaluations = 10000;
options.FunctionTolerance = 1e-25;
options.ConstraintTolerance =  1e-8;
options.StepTolerance =  1e-16; 
options.Display = 'off';
options.SpecifyObjectiveGradient = true(1);
options.OptimalityTolerance = 1e-20; %NIR
options.FiniteDifferenceStepSize = 1e-25;
if dub_mode
    options.Display = 'iter'; %'off'; %'final';
end

% solving
[x, ~, ~, O] = fmincon(fun,x0,A,b,Aeq,beq,[],[],[],options);
%[x, fval2, exitflag2, output2] = fmincon(fun,x,A,b,[],[],[],[],[],options);

% concluding
x_est    = x(1:L);
est_dist = x((L+1):2*L); 

if dub_mode
    x_est = align_to_reference(x_est, true_x);
    [v_stop,g_stop] = costFunctionV2(x(1:L), x((L+1):2*L), mu, P, lambda, sigma);
    [v_real,g_real] = costFunctionV2(true_x, true_rho, mu, P, lambda, sigma);
end

end
