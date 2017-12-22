function [results] = run_template_attack_adapt(...
    s_profile, s_helper, ...
    m_data_attack, metadata_attack, idx_attack, ...
    bytes, atype, cmethod, cparams, discriminant, ...
    rand_iter, nr_traces_vec, eparams)
%RUN_TEMPLATE_ATTACK_ADAPT Runs a template attack
%   [results] = RUN_TEMPLATE_ATTACK_ADAPT(...
%       s_profile, s_helper, ...
%       m_data_attack, metadata_attack, idx_attack, ...
%       bytes, atype, cmethod, cparams, discriminant, ...
%       rand_iter, nr_traces_vec, eparams)
%   runs a template attack with the given parameters and returns a results
%   structure that is defined below.
%
%   This method is somewhat similar to run_template_attack but it adapts the
%   profiling data with the aim of matching better the attack data.
%   This can be used to evaluate a real scenario where the attacker is
%   forced to use different devices and can only acquire a limited subset
%   of traces from the target device. The aim of the attacker therefore is
%   to use the available attack traces and as much training data as possible
%   (possibly from a different device) to adapt the profiling data in order
%   to match as well as possible the data in the attack traces.
%
%   s_profile is a structure with data for profiling. This structure should
%   contain traces of the same length and generally over the same
%   operations as the attack traces. This structure should have the
%   following fields:
%   - .nr_sets: the number of profiling data sets available
%   - .mmap_data: a cell of length nr_sets, containing the memory mapped
%   data sets.
%   - .metadata: a cell of length nr_sets, containing the metadata related
%   to the memory mapped object. Use get_mmap on the selected data to
%   obtain these objects.
%   - .idx: a cell of length nr_sets, containing the vectors of
%   indeces specifying which traces from each group should be used for a
%   particular training set.
%
%   s_helper is a structure with helper data for profiling. This might be
%   for example, traces from the target device over instructions that can
%   be manipulated by the attacker and that may aid in the improvement
%   of templates. Its use depends on the attack type 'atype'. If in doubt,
%   pass [] (empty).
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
%   - 'mvn': this is just the standard multivariate attack, no adaptation.
%   - 'famvn': standard multivariate attack, but using factor analysis for
%   covariance.
%   - 'mvn_offset_median': using the median from attack traces to compensate
%   offset in profiling data for a template attack based on mvn.
%   - 'multi': combining traces from multiple boards during profile, then
%   running a normal mvn template attack on another dataset (board).
%   - 'multi_offset_median': combining traces from multiple boards during
%   profile and using the median from attack traces to compensate
%   offset in profiling data for a template attack based on mvn.
%   - 'roffset': adding a random offset to each trace before profiling in
%   order to train to avoid such offset. Values given in eparams. Attack
%   using mvn.
%   - 'roffset_median': similar to 'roffset' but in addition it also
%   compensates for the offset as in 'mvn_offset_median'.
%   - 'boffset': adds random offset to the M and B matrices, allowing PCA
%   to be more efficient across different campaigns.
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
%       -> cparams.p2 is an additional optional parameter.
%   - 'PCA':
%       -> cparams.pca_threshold
%       -> cparams.pca_alternate
%       -> cparams.pca_dimensions
%   - 'LDA':
%       -> cparams.lda_threshold
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
%   used for each evaluation step.
%
%   eparams is a structure of extra parameters that may be needed for some
%   options. Some of these extra parameters include:
%   - 'nr_factors': number of factors to use for factor analysis (which is
%   enabled using the 'famvn' atype).
%   - 'roffset': cell of length nr_sets, having a matrix of size
%   np_set x nr_points per set defining the exact random offset to be added
%   to each trace when using the atype 'roffset'.
%   - 'save_xdata': if this field exists and is non zero then the x_profile
%   and x_attack data will be saved. Otherwise, these will not be saved as
%   they take a considerable amount of space.
%   - 'save_eval': if this field exists and is non zero then the
%   handle_eval pointer and the data for evaluation (generally the
%   templates) will be saved. Otherwise, these will not be saved as
%   they take a considerable amount of space.
%   - 'use_elv': for both PCA and LDA, use the ELV method for determining
%   the eigenvectors to be used instead of the default cummulative
%   percentage method.use_eig
%
%   The 'results' structure contains the following:
%   -> results.s_profile: the structure for profiling data.
%   -> results.s_helper: the structure for helper data.
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
%   -> results.L: matrix of factor loadings (only for 'famvn')
%   -> results.P: matrix of individual factors (only for 'famvn')
%   -> results.x_profile: the profiling data, after compression. (optional)
%   -> results.x_attack: the attack data, after compression. (optional)
%   -> results.handle_prepare: function used to extract features for templates.
%   -> results.pp1 ... results.pp5: parameters of results.handle_prepare.
%   -> results.handle_eval: function used to evaluate templates. (optional)
%   -> results.pe3 ... results.pe6: parameters for results.handle_eval. (optional)
%   -> results.success_info: guessing entropy information, as returned by
%   the get_success_info_like method.
%
%   See the paper "Efficient Template Attacks" by Omar Choudary and Markus
%   Kuhn, presented at CARDIS 2013.

%% Check and initialise parameters
nr_groups = length(bytes);
results.s_profile = s_profile;
results.s_helper = s_helper;
results.metadata_attack = metadata_attack;
results.idx_attack = idx_attack;
results.bytes = bytes;
results.atype = atype;
results.cmethod = cmethod;
results.cparams = cparams;
results.discriminant = discriminant;
results.rand_iter = rand_iter;
results.nr_traces_vec = nr_traces_vec;
if nargin < 13
    eparams = [];
end
results.eparams = eparams;
fprintf('Running run_template_attack_adapt() ...\n');

%% Get the sums of squares and cross products
fprintf('Obtaining sums of squares and cross products for all sets ...\n');
if strcmp(atype, 'roffset') || strcmp(atype, 'roffset_median')
    [M, B, W, np] = compute_ssp_e2_mmap_multi(s_profile, bytes, [], eparams.roffset);
elseif strcmp(atype, 'boffset')
    [M, ~, W, np] = compute_ssp_e2_mmap_multi(s_profile, bytes);
    [m,n] = size(M);
    xm = mean(M, 1);
    R = mean(xm)*rand(m,1)*ones(1,n);
    M = M + R;
    xm2 = mean(M, 1);
    XB = M - ones(m,1)*xm2;
    B = XB'*XB;
else
    [M, B, W, np] = compute_ssp_e2_mmap_multi(s_profile, bytes);
end
xmm = mean(M, 1);
results.M = M;
results.B = B;
results.W = W;
toc
    
%% Estimate correlation via factor analysis if 'famvn' specified
if strcmp(atype, 'famvn')
    if ~isfield(eparams, 'nr_factors')
        error('Need eparams.nr_factors for famvn');
    else
        nr_factors = eparams.nr_factors;
    end
    C = W / (nr_groups*(np - 1));
    dvec = sqrt(diag(C));
    Dinv = diag(1./dvec);
    R = Dinv*C*Dinv;
    [U, S, ~] = svd(R);
    d = diag(S);
    L = U(:,1:nr_factors)*diag(sqrt(d(1:nr_factors)));
    results.L = L;
    P = R - L*L';
    P = diag(P);
    results.P = P;
    RE = L*L' + diag(P); % Estimated correlation matrix from factor analysis
    CE = diag(dvec) * RE * diag(dvec); % estimated covariance from RE
    W = CE * (nr_groups*(np - 1)); % Replace W with data from CE
end

%% Get compression parameters
if strcmp(cmethod, 'sample')
    fprintf('Computing selection curves and selected samples ...\n');
    [curves] = get_signal_strength_ssp(M, B, W, np);
    
    if ~isfield(cparams, 'p2')
        cparams.p2 = [];
    end
    interest_points = get_selection(curves.(cparams.curve), cparams.sel, ...
        cparams.p1, cparams.p2);
    
    handle_prepare = @prepare_data_template_xmm;
    pp1 = interest_points;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
elseif strcmp(cmethod, 'PCA')
    fprintf('Computing PCA parameters...\n');
    [U, D, xmm, K] = compute_params_pca(M, cparams.pca_threshold, ...
                                        cparams.pca_alternate);
    if ~isempty(eparams) && isfield(eparams, 'use_elv')
        params = [];
        params.method = 'maximal';
        params.max_elvs = round(size(U,1) / 100);
        [idx] = get_elv_order(U, D, params);
        U = U(:,idx);
    end
    if isfield(cparams, 'pca_dimensions') && (cparams.pca_dimensions > 0)
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
    [A ,D, K] = compute_params_lda(B, Spool, nr_groups, cparams.lda_threshold);
    if ~isempty(eparams) && isfield(eparams, 'use_elv')
        params = [];
        params.method = 'maximal';
        params.max_elvs = round(size(A,1) / 100);
        [idx] = get_elv_order(A, D, params);
        A = A(:,idx);
    end
    if isfield(cparams, 'lda_dimensions') && (cparams.lda_dimensions > 0)
        FW = A(:,1:cparams.lda_dimensions);
    else
        FW = A(:,1:K);
    end

    handle_prepare = @prepare_data_template_pca_v2;
    pp1 = FW;
    pp2 = xmm;
    pp3 = [];
    pp4 = [];
    pp5 = [];
else
    error('Unknown compression method: %s', cmethod);
end

%% Store handle_prepare data
results.handle_prepare = handle_prepare;
results.pp1 = pp1;
results.pp2 = pp2;
results.pp3 = pp3;
results.pp4 = pp4;
results.pp5 = pp5;

%% Load raw leakage data for profile
fprintf('Computing profiling data from all sets...\n');
if strcmp(atype, 'roffset') || strcmp(atype, 'roffset_median')
    [x_profile] = compute_features_e2_mmap_multi(s_profile, ...
                                    handle_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5, bytes, ...
                                    eparams.roffset);
else
    [x_profile] = compute_features_e2_mmap_multi(s_profile, ...
                                    handle_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5, bytes);
end
if isfield(eparams, 'save_xdata') && (eparams.save_xdata ~= 0)
    results.x_profile = x_profile;
end
toc

%% Load raw leakage data for attack
fprintf('Computing attack data...\n');
if strcmp(atype, 'mvn_offset_median') || strcmp(atype, 'multi_offset_median') ...
        || strcmp(atype, 'roffset_median')
    s_adapt = [];
    s_adapt.type = 'offset_median';
    s_adapt.xmm = xmm;
    [~, x_attack] = evalc(['compute_features_e2_mmap_adapt(', ...
                            'm_data_attack, metadata_attack, idx_attack,', ...
                            's_adapt,', ...
                            'handle_prepare, pp1, pp2, pp3, pp4, pp5, bytes)']);
else
    [~, x_attack] = evalc(['compute_features_e2_mmap(', ...
                        'm_data_attack, metadata_attack, idx_attack,', ...
                        'handle_prepare, pp1, pp2, pp3, pp4, pp5, bytes)']);
end
if isfield(eparams, 'save_xdata') && (eparams.save_xdata ~= 0)
    results.x_attack = x_attack;
end
toc

%% Compute templates
if strcmp(atype, 'mvn') || strcmp(atype, 'mvn_offset_median') ...
        || strcmp(atype, 'multi') || strcmp(atype, 'multi_offset_median') ...
        || strcmp(atype, 'roffset') || strcmp(atype, 'roffset_median') ...
        || strcmp(atype, 'boffset')
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
        error('discriminant not supported for mvn: %s', discriminant);
    end
elseif strcmp(atype, 'famvn')
    fprintf('Computing mvn template and evaluation parameters...\n');
    [tmiu, tsigma] = compute_template(x_profile);
    
    handle_eval = @evaluate_discriminant;
    if strcmp(discriminant, 'linear')
        c0 = mean(tsigma, 3);
        [U, S, ~] = svd(c0);
        d = diag(S);
        L = U(:,1:nr_factors)*diag(sqrt(d(1:nr_factors)));
        P = c0 - L*L';
        P = diag(P);
        c0_f = L*L' + diag(P);
        ic0 = inv(c0_f);
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
            [U, S, ~] = svd(tsigma(:,:,k));
            d = diag(S);
            L = U(:,1:nr_factors)*diag(sqrt(d(1:nr_factors)));
            P = tsigma(:,:,k) - L*L';
            P = diag(P);
            sk_f = L*L' + diag(P);
            tsinv(:,:,k) = inv(sk_f);
            tlogdet(k) = logdet(sk_f, 'chol');
        end
        pe3 = tmiu;
        pe4 = tsinv;
        pe5 = tlogdet;
        pe6 = [];    
    else
        error('discriminant not supported for famvn: %s', discriminant);
    end
else
    error('template attack type not supported: %s', atype);
end

%% Store evaluation data if requested
if isfield(eparams, 'save_eval') && (eparams.save_eval ~= 0)
    results.handle_eval = handle_eval;
    results.pe3 = pe3;
    results.pe4 = pe4;
    results.pe5 = pe5;
    results.pe6 = pe6;
end

%% Compute the success information
fprintf('Computing success info...\n');
[results.success_info] = get_success_info_like(x_attack, rand_iter, ...
                                       nr_traces_vec, ...
                                       handle_eval, pe3, pe4, pe5, pe6);
toc

end

