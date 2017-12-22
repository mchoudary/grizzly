%% Show plots of template attack results
% 
% Author: Omar Choudary

%% Reset environment
% close all;
clear;
set(0, 'DefaulttextInterpreter', 'latex'); % Use TeX interpretation
tic

%% Setup paths and parameters
rpath = 'figures/';
font_size = 12;
options = 'gp';
yrange = [0, 6];
style = 'cmapjet';
nr_traces_vec = [1:10, 20:10:100, 200, 500, 1000];
len_na_vec = length(nr_traces_vec);


%% Load template results and related data
fdata = 'results/a2_bat_fb_stochastic_f9_all_dlinear_n200r_slr_g1000_r10.mat';
data = load(fdata, 'results');
results = data.results;
np = size(results{1}.idx_profile, 1);
rand_iter = results{1}.rand_iter;
nr_exp_ab = length(results);
L = cell(nr_exp_ab, 1);
G = zeros(nr_exp_ab, len_na_vec);
slines_ab = cell(nr_exp_ab, 1);
style = 'fancy';
for k=1:nr_exp_ab
    if strcmp(results{k}.atype, 'templatelda')
        L{k} = 'A2 BAT FB, DLINEAR, S-LDA';
    elseif strcmp(results{k}.atype, 'templatepca')
        L{k} = 'A2 BAT FB, DLINEAR, S-PCA';            
    else
        L{k} = 'A2 BAT FB, DLINEAR, SELECTION';
    end
    g = get_ge_from_success_info(results{k}.success_info, nr_traces_vec);
    G(k,:) = g.joint;
    slines_ab{k} = get_line_properties_templates(k, style);
end

%% Plot results
title_results = sprintf('A2 BAT FB, np=%d', np);
rprefix = sprintf('a2_bat_fb_stoc_dlinear_n%dr_ls_r%d_', np, rand_iter);
make_figures_ge(G, nr_traces_vec, ...
                     rpath, rprefix, ...
                     title_results, L, font_size, ...
                     slines_ab, options, yrange);
                 
%% Exit the matlab interpreter, needed when run in script mode
% exit

