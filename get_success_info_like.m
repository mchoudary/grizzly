function [success_info] = get_success_info_like(test_data, ...
                                                nr_iter, ...
                                                nr_traces_vec, ...
                                                func_discriminant, ...
                                                varargin)
%GET_SUCCESS_INFO_FUNC Returns success rate information from some data
%   [success_info] = GET_SUCCESS_INFO_LIKE(test_data, ...
%                                          nr_iter, ...
%                                          nr_traces_vec, ...
%                                          func_discriminant, ...
%                                          varargin)
%   computes several measures of success rate based on a discriminant
%   function func_discriminant (e.g. evaluate_pdf) and some test data.
%
%   This function is similar to get_success_info_like but uses a slightly
%   different function (func_discriminant), to evaluate the discriminant score
%   of a trace or set of traces jointly rather than independently.
%   
%   Another important difference is that in this method only the depth
%   values over each iteration are kept, discarding the individual
%   values. Kepping the depth values over each iteration allows
%   a better estimation of the standard error of the guessing entropy.
%
%   test_data should contain the test data which has been already
%   preprocessed by a suitable function, such as prepare_data_template.
%   The test_data must be a matrix of size
%   nr_traces x nr_interest_points x nr_groups.
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
%   sample vector or matrix of samples (one per row), given some parameters
%   (varargin), the sample vector/matrix on which to compute the
%   discriminant (p1) and a particular class or classes (one per row)
%   for which to compute the result.
%
%   The main difference of func_discriminant, compared to func_eval in
%   get_success_info_func, is that for multiple traces the func_discriminant
%   should compute the joint likelihood of the given parameters for that
%   combination of target traces.
%
%   Regardless of the number of target traces to be used in the evaluation
%   of func_discriminant, it is possible to pass the other
%   parameters (p2, ...) for multiple classes. In that case, the result
%   will be computed for the given samples for all the classes.
%   If the given function needs less parameters just pass empty for rest.
%
%   Note that func_discriminant only needs to
%   return values such that we can sort the results in decreasing order
%   of their "likelihood", that is the func_discriminant does not have to be a
%   strict likelihood/pdf/pmf. All we need is a function that provides
%   relative order, so even a log-likelihood works fine.
%
%   An example of func_discriminant function is evaluate_mvn_like_log.
%
%   success_info is a structure containing the results of this
%   method. success_info has these substructures:
%   - depth: matrices of size nr_groups x nr_iter, containing the index
%   of the correct group when sorting the likelihoods for each iteration.
%   This data can be used to compute the guessing entropy and success rate.
%   Note that the actual data is found under the 'group(x)' substructures,
%   as explained below.
%   - rindex: matrices of size nr_traces x nr_iter, containing the indeces
%   of the random selection for each group size.
%
%   The structure depth (see above) has the following substructures:
%   - avg: results from average of selected attack traces taken as one trace
%   (i.e removing noise as much as possible)
%   - joint: results from computing the joint likelihood on all selected
%   attack traces.
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

%% Initialize and check parameters
nr_traces = size(test_data, 1);
nr_groups = size(test_data, 3);
nr_test_groups = length(nr_traces_vec);
good_class = 1:nr_groups;

%% Compute the success info for each test group size
for i=1:nr_test_groups
    fprintf('Computing success info for group size %d\n', ...
             nr_traces_vec(i));
    
    %% Set up the depth vectors
    success_info.depth.avg.(['group' num2str(i)]) = zeros(nr_groups, nr_iter);
    success_info.depth.joint.(['group' num2str(i)]) = zeros(nr_groups, nr_iter);
    
    %% Select random traces to use for all the tests and iterations
    rindex = randi([1, nr_traces], nr_traces_vec(i), nr_iter);
    success_info.rindex.(['group' num2str(i)]) = rindex;
    
    %% Compute success info for each group
    for group=1:nr_groups        
        %% Perform the tests for nr_iter
        for count=1:nr_iter  
            %% Select data
            data = test_data(rindex(:,count), :, group);
            data_avg = (1/nr_traces_vec(i))*(ones(1,nr_traces_vec(i))*data);            

            %% First compute likelihood values
            l_avg = func_discriminant(data_avg, (1:nr_groups)', varargin{:});
            l_joint = func_discriminant(data, (1:nr_groups)', varargin{:});

            %% Then compute depth vectors
            % The main assumption here is that the group corresponding to
            % the highest likelihood is the correct group. Therefore make
            % sure that the function you have passed as func_discriminant
            % has this property.
            [~, si] = sort(l_avg(:), 1, 'descend');
            success_info.depth.avg.(['group' num2str(i)])(group, count) = ...
                find(si == good_class(group));
            [~, si] = sort(l_joint(:), 1, 'descend');
            success_info.depth.joint.(['group' num2str(i)])(group, count) = ...
                find(si == good_class(group));
        end
    end
end

end
    
    