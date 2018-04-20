function [success_info] = get_success_info_generic(...
                                        X_attack, V_attack, V_discriminant, ...
                                        nr_iter, nr_traces_vec, ...
                                        func_discriminant, varargin)
%GET_SUCCESS_INFO_GENERIC Returns success information from generic data
%   [success_info] = GET_SUCCESS_INFO_GENERIC(...
%                                   X_attack, V_attack, V_discriminant, ...
%                                   nr_iter, nr_traces_vec, ...
%                                   func_discriminant, varargin)
%   computes several measures of success rate based on a discriminant
%   function func_discriminant (e.g. evaluate_pdf) and generic test data.
%
%   X_attack should be a matrix of size nr_traces x nr_samples x nr_groups
%   having the test data which has been already
%   preprocessed (e.g. via prepare_data_template) to the same format as the
%   varargin parameters passed for func_discriminant (see below).
%
%   V_attack should be a vector of length nr_groups providing the data
%   values corresponding to the traces in each group of X_attack. That is,
%   the same value in V_attack will apply to all the nr_traces in each
%   group of X_attack.
%
%   V_discriminant should be a vector of length nr_disc_values specifying the
%   values corresponding to the parameters (e.g. means) used in the
%   evaluation of func_discriminant. This is needed to estimate measures
%   such as the guessing entropy and depth matrix.
%
%   nr_iter specifies how many times to run each test on
%   randomly picked samples to produce the desired results.
%   Note that the execution time increases linearly with this parameter.
%
%   nr_traces_vec: the groups of traces to test. nr_traces_vec
%   should be a vector of integers, where each element represents the
%   size of one test group. 
%   E.g. if nr_traces_vec = [5, 10] then 2 * rand_iter tests will be
%   performed using averaged traces, 2 * rand_iter tests using independent
%   traces, etc. The first test will perform rand_iter iterations using
%   5 random samples for all cases (independent traces, averaged traces,
%   etc.) and the second test will use 10 random samples. The same random
%   samples at each rand_iter iteration are used for all the cases, i.e.
%   for averaged traces, for independent traces, etc.
%
%   func_discriminant should be a method that returns a discriminant of a
%   sample vector or matrix of samples (one per row), given the sample
%   vector/matrix on which to compute the discriminant (p1) and some
%   parameters (varargin),
%
%   Regardless of the number of target traces to be used in the evaluation
%   of func_discriminant, it is possible to pass the other
%   parameters (p2, ...) for multiple classes. The result will be computed
%   for the given samples for all the classes given in the parameters.
%   If the given function needs less parameters just pass empty for rest.
%
%   Note that func_discriminant only needs to
%   return values such that we can sort the results in decreasing order
%   of their "likelihood", that is the func_discriminant does not have to be a
%   strict likelihood/pdf/pmf. All we need is a function that provides
%   relative order, so even a log-likelihood works fine.
%
%   An example of func_discriminant function is compute_discriminant.
%
%   success_info is a structure containing the results of this
%   method. success_info has these substructures:
%   - depth: matrices of size nr_attack_groups x nr_iter, containing the index
%   of the correct group when sorting the joint likelihoods for each iteration.
%   This data can be used to compute the guessing entropy and success rate
%   but also the partial guessing entropy that I defined in my work.
%   Note that the actual data is found under the 'group(x)' substructures,
%   as explained below.
%   - rindex: matrices of size nr_traces x nr_iter, containing the indeces
%   of the random selection for each group size.
%   
%   All structures (success_info.depth.avg, success_info.rindex, etc.) have
%   substructures corresponding to the nr_traces_vec
%   parameters. These structures are labeled group1, group2, etc...
%   For example, if nr_avg_traces = [1 10 50] then leakage_info.depth.avg
%   would have the following substructures:
%       - group1: results for tests with nr_avg_traces(1) = 1 avg traces
%       - group2: results for tests with nr_avg_traces(2) = 10 avg traces
%       - group3: results for tests with nr_avg_traces(3) = 50 avg traces
%
%   These last structures in turn contain the actual data.
%
%   Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)
%
%   See also compute_discriminant.

%% Initialize and check parameters
nr_traces = size(X_attack, 1);
nr_groups = size(X_attack, 3);
nr_test_groups = length(nr_traces_vec);

%% Compute the success info for each test group size
for i=1:nr_test_groups
    fprintf('Computing success info for group size %d\n', ...
             nr_traces_vec(i));
    
    %% Set up the depth vectors
    success_info.depth.(['group' num2str(i)]) = zeros(nr_groups, nr_iter);    
    
    %% Select random traces to use for all the tests and iterations
    rindex = zeros(nr_traces_vec(i), nr_iter);
    for k=1:nr_iter
        rindex(:,k) = randperm(nr_traces, nr_traces_vec(i));
    end   
    success_info.rindex.(['group' num2str(i)]) = rindex;
    
    %% Compute success info for each group
    for group=1:nr_groups        
        %% Perform the tests for nr_iter
        for count=1:nr_iter  
            %% Select data
            data = X_attack(rindex(:,count), :, group);

            %% First compute likelihood values
            l_joint = func_discriminant(data, varargin{:});

            %% Then compute depth vectors
            % The main assumption here is that the group corresponding to
            % the highest likelihood is the correct group. Therefore make
            % sure that the function you have passed as func_discriminant
            % has this property.
            [~, si] = sort(l_joint(:), 1, 'descend');
            success_info.depth.(['group' num2str(i)])(group, count) = ...
                find(V_discriminant(si) == V_attack(group));
        end
    end
    toc
end

end
    
    