function [xdata] = compute_features_generic_mmap(map, D_all, ...
                                    V_sel, idx_group, ...
                                    func_prepare, ...
                                    pp1, pp2, pp3, pp4, pp5)
%COMPUTE_FEATURES_GENERIC_MMAP Computes features from sample vectors
%   [xdata] = COMPUTE_FEATURES_GENERIC_MMAP(map, bindex, ...
%                              func_prepare, ...
%                              pp1, pp2, pp3, pp4, pp5, bytes)
%   computes some features from sample vectors. This method is targeted
%   to the E2-related experiments, where a series of bytes are loaded to
%   some registers, and generally the value of a particular byte is being
%   analysed while the values of the other bytes are either fixed or
%   random.
%
%   This method will apply the func_prepare method (details below) on the
%   data represented by the memmapfile object map, in order to compute some
%   features from each sample (trace) in the input data. In addition, this
%   method will split the data into the different group categories
%   available, such that the output of this method will be a matrix xdata
%   of size nr_rows x nr_features x nr_groups.
%
%   map should be a memmapfile object with 1 data entry containing the data
%   matrix X of size nr_points x nr_trials.
%
%   D_all should be a vector of length nr_trials containing the values
%   corresponding to each trace in X. This method assumes there will be at
%   least length(idx_group) values for each V_sel in D_all.
%
%   V_sel should be a vector of length nr_values that specifies the values
%   (out of those in D_all) for which the features specified by
%   func_prepare will be extracted. For each value in V_sel this method
%   will first find the index of all traces that correspond to that value
%   and then will select only those traces specified by idx_group.
%
%   idx_group should be a vector of length nr_traces specifying the indeces
%   of the traces to use, from the set of all traces per each value/group.
%
%   The func_prepare method should transform the data (e.g. compress it)
%   for the further processing. The output data can then be used to test
%   leakage, classification, etc. An exmaple of func_prepare method is
%   prepare_data_pca which uses some PCA parameters to retain only some
%   principal directions of the data. The prototype of func_prepare is:
%   data_out = func_prepare(data_in, pp1, pp2, pp3, ...) where data_in is
%   the raw input trace and pp1, pp2, pp3 are parameters needed by the
%   specific function.
%
%   The output xdata is a matrix of size nr_traces x nr_features x nr_groups,
%   where nr_rows is the number of traces per group,
%   nr_features is the number of features remaining per trace after
%   applying the func_prepare method and nr_groups is the number of groups
%   specified in the metadata.
%
%   See also compute_features, get_mmap.

%% Check and initialize data
nr_values = length(V_sel);
nr_traces = length(idx_group);
D_all = D_all(:);

%% Extract data for each value
[D_sorted, si] = sort(D_all, 1, 'ascend');
for k=1:nr_values
    k_first = find(D_sorted == V_sel(k), 1, 'first');
    lindex = si(k_first + idx_group - 1);
    
    % Below note transpose and conversion in case we had integer class
    L = double(map.data(1).X(:,lindex)');
    
    data_out = func_prepare(L, pp1, pp2, pp3, pp4, pp5);    
    if k==1        
        nr_features = size(data_out, 2);
        xdata = zeros(nr_traces, nr_features, nr_values);
    end    
    xdata(:,:,k) = data_out;
end
                              
end
    
