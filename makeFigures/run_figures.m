% script name: "run_figures"
%
% This script generates the first four figures of the paper
% 
% To modify the setting, naming, and saving, open each script separately 
%
% NS, September -- November  2017

% the counter example of uniquness with periodic distribution, i.e., Proposition III.6, including figures
make_figure_counter_example
% the comparison between : spectral, LS, and adapted-EM
vary_sigma_comparison
% in the following, make sure the number of trials is large enough. Add parallel pool if needed
run_EM_comparison


% older comparisons: main_comparison, plot_EMvsLS



