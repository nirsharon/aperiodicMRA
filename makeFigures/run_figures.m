% script name: "run_figures"
%
% This script generates the first four figures of the paper
% 
% To modify the setting, naming, and saving, open each script separately 
%
% NS, September -- November  2017

% Figure I.1
visual_fig_data
% the counter example of uniquness with periodic distribution, i.e., Proposition III.6, including Figure III.1
make_figure_counter_example
% the basis figures and calculations for Figure VII.1
visual_demo_LS
% in the following, make sure the number of trials is large enough. Add parallel pool if needed (Figure VII.2)
run_EM_comparison
% the comparison between : spectral, LS, and adapted-EM (Figure VII.3)
vary_sigma_comparison


% older comparisons: main_comparison, plot_EMvsLS



