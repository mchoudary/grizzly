function [results] = run_template_attack(...
    m_data_profile, metadata_profile, idx_profile, ...
    m_data_attack, metadata_attack, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, ...
    rand_iter, nr_traces_vec)
%RUN_TEMPLATE_ATTACK Runs a template attack
%   [results] = RUN_TEMPLATE_ATTACK(...
%       m_data_profile, metadata_profile, idx_profile, ...
%       m_data_attack, metadata_attack, idx_attack, ...
%       bytes, atype, cmethod, cparams, discriminant, ...
%       rand_iter, nr_traces_vec)
%   runs a template attack with the given parameters and returns a results
%   structure that is defined below. This method is intended for any memory
%   mapped data that has a similar structure as the data from the E2
%   experiment. See get_mmap for more details.
%
%   m_data_profile and metadata_profile should be the memory mapped object
%   and associated metadata info for the profiling data. Use get_mmap 
%   on the selected data to obtain these objects.
%
%   idx_profile should be a vector of indices specifying which traces from
%   each group should be used for the profile data.
%
%   m_data_attack and metadata_attack should be the memory mapped object
%   and associated metadata info for the attack data. Use get_mmap 
%   on the selected data to obtain these objects. You can pass the same
%   objects for profiling and attack, which should make the attack faster.
%
%   idx_attack should be a vector of indices specifying which traces from
%   each group should be used for the attack data.
%   
%   bytes should be a vector of indices specifying which bytes (starting
%   from 0) will be used for the attack. This might be useful in order to
%   restrict the attack only to the bytes 0-15 for example (i.e. using 4
%   bits).
%
%   atype should be a string specifying the type of template attack to be
%   used. Currently supported are:
%   - 'mvn': which relies on the multivariate normal probability density
%   function to compute templates and probabilities.
%
%   cmethod should be a string specifying the compression method. Currently
%   supported methods are: 'sample', 'PCA' and 'LDA'.
%
%   cparams should be a structure of params specific to the compression
%   method. For each compression method the params are as follows:
%   - 'sample':
%       -> cparams.curve is a string specifying the signal strength curve to
%       be used.
%       -> cparams.sel is a string specifying the class of selection.
%       -> cparams.p1 is a parameter for the class of selection.
%   - 'PCA':
%       -> cparams.pca_threshold
%       -> cparams.pca_alternate
%       -> cparams.pca_dimensions
%   - 'LDA':
%       -> cparams.lda_dimensions
%
%   discriminant should be a string specifying the type of discriminant to
%   be used. The possible options are:
%   - 'linear': uses a pooled common covariance matrix with a linear
%      discriminant.
%   - 'linearnocov': does not use a covariance matrix. Might be useful in
%      particular with LDA, where the covariance should be the
%      identity if the eigenvectors are chosen carefully.
%   - 'log': uses individual covariances and log-determinants to compute
%      the group specific templates (mean and covariance).
%
%   rand_iter should be a positive integer specifying the number of
%   iterations to run the evaluation (guessing_entropy) computation. The
%   returned results may contain either the individual or the average
%   results. Check below for details.
%
%   nr_traces_vec is a vector containing the number of attack traces to be
%   used for each element.
%
%   The 'results' structure contains the following:
%   -> results.metadata_profile: the metadata structure for profile.
%   -> results.idx_profile: the idx_profile vector.
%   -> results.metadata_attack: the metadata structure for attack.
%   -> results.idx_attack: the idx_attack vector.
%   -> results.bytes: the bytes vector.
%   -> results.atype: the atype string.
%   -> results.cmethod: the compression method string.
%   -> results.cparams: the cparams structure.
%   -> results.discriminant: the discriminant string.
%   -> results.rand_iter: the number of iterations.
%   -> results.nr_traces_vec: the vector with number of attack traces.
%   -> results.M: the matrix of group means.
%   -> results.B: the between-groups matrix.
%   -> results.W: the matrix of variances and covariances across all data.
%   -> results.x_profile: the profiling data, after compression.
%   -> results.x_attack: the attack data, after compression.
%   -> results.success_info: guessing entropy information, as returned by
%   the get_success_info_like method.
%
%   See the paper "Efficient Template Attacks" by Omar Choudary and Markus
%   Kuhn, presented at CARDIS 2013.

%% Check and initialise parameters
np = length(idx_profile);
nr_groups = length(bytes);
results.metadata_profile = metadata_profile;
results.idx_profile = idx_profile;
results.metadata_attack = metadata_attack;
results.idx_attack = idx_attack;
results.bytes = bytes;
results.atype = atype;
results.cmethod = cmethod;
results.cparams = cparams;
results.discriminant = discriminant;
results.rand_iter = rand_iter;
results.nr_traces_vec = nr_traces_vec;

%% Get the sums of squares and cross products
fprintf('Obtaining sums of squares and cross products...\n');
[M, B, W] = compute_ssp_e2_mmap(m_data_profile, metadata_profile, ...
                                idx_profile, bytes);
results.M = M;
results.B = B;
results.W = W;
xmm = mean(M, 1);
toc

%% Get compression parameters
if strcmp(cmethod, 'sample')
    fprintf('Computing selection curves and selected samples ...\n');
    [curves] = get_signal_strength_ssp(M, B, W, np);
    
    interest_points = get_selection(curves.(cparams.curve), cparams.sel, cparams.p1);
    
    handle_prepare = @prepare_data_template;
    pp1 = interest_points;
    pp2 = [];
    pp3 = [];
    pp4 = [];
    pp5 = [];
elseif strcmp(cmethod, 'PCA')
    fprintf('Computing PCA parameters...\n');
    [U, D, xmm, K] = compute_params_pca(M, cparams.pca_threshold, ...
                                        cparams.pca_alternate);
    if cparams.pca_dimensions > 0
        U = U(:,1:cparams.pca_dimensions);
    else
        U = U(:,1:K);
    end
    
    handle_prepare = @prepare_data_template_pca_v2;
    pp1 = U;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
elseif strcmp(cmethod, 'LDA')
    fprintf('Computing Fishers LDA parameters...\n');
    Spool = W / (nr_groups*(np-1));
    [A,D] = compute_params_lda(B, Spool);
    FW = A(:,1:cparams.lda_dimensions);

    handle_prepare = @prepare_data_template_pca_v2;
    pp1 = FW;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
else
    error(sprintf('Unknown compression method: %s', cmethod));    
end

%% Load raw leakage data for profile
% Using evalc to avoid console output
fprintf('Computing profiling data...\n');
[~, x_profile] = evalc(['compute_features_e2_mmap(', ...
                        'm_data_profile, metadata_profile, idx_profile,', ...
                        'handle_prepare, pp1, pp2, pp3, pp4, pp5, bytes)']);
results.x_profile = x_profile;
toc

%% Load raw leakage data for attack
fprintf('Computing attack data...\n');
[~, x_attack] = evalc(['compute_features_e2_mmap(', ...
                        'm_data_attack, metadata_attack, idx_attack,', ...
                        'handle_prepare, pp1, pp2, pp3, pp4, pp5, bytes)']);
results.x_attack = x_attack;
toc

%% Compute templates
if strcmp(atype, 'mvn')
    fprintf('Computing mvn template and evaluation parameters...\n');
    [tmiu, tsigma] = compute_template(x_profile);
    
    handle_eval = @evaluate_discriminant;
    if strcmp(discriminant, 'linear')
        c0 = mean(tsigma, 3);
        ic0 = inv(c0);
        pe3 = tmiu;
        pe4 = ic0;
        pe5 = [];
        pe6 = [];
    elseif strcmp(discriminant, 'linearnocov')
        pe3 = tmiu;
        pe4 = [];
        pe5 = [];
        pe6 = [];
    elseif strcmp(discriminant, 'log')
        n = size(tsigma, 3);
        tsinv = zeros(size(tsigma));
        tlogdet = zeros(n, 1);
        for k = 1:n
            tsinv(:,:,k) = inv(tsigma(:,:,k));
            tlogdet(k) = logdet(tsigma(:,:,k), 'chol');
        end
        pe3 = tmiu;
        pe4 = tsinv;
        pe5 = tlogdet;
        pe6 = [];
    else
        error('discriminant not supported: %s', discriminant);
    end 
else
    error('template attack type not supported: %s', atype);
end

%% Compute the success information
fprintf('Computing success info...\n');
[results.success_info] = get_success_info_like(x_attack, rand_iter, ...
                                       nr_traces_vec, ...
                                       handle_eval, pe3, pe4, pe5, pe6);
toc


