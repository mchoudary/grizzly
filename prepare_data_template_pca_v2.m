function [data_out] = prepare_data_template_pca_v2(data_in, W, xmm, varargin)
%PREPARE_DATA_TEMPLATE_PCA_V2 Prepares data for leakage analysis
%   [data_out] = PREPARE_DATA_TEMPLATE_PCA_V2(data_in, W, xmm, varargin)
%   transforms the input trace(s) into the required data for leakage
%   analysis by means of the given parameters.
%
%   This function projects the data into a PCA subspace by means of the
%   eigenvectors W and the average mean xmm.
%
%   data_in should be a matrix of size nr_traces x nr_points.
%
%   W should be a matrix of size nr_points x K, where K is the number
%   of dimensions that were retained during the computation of the PCA
%   parameters.
%
%   xmm is a vector of length nr_points containing the average of the
%   mean traces over all the data groups. This will be used to normalize
%   the data (by substracting xmm from each trace) before projecting it
%   unto the PCA space.
%
%   The output is a matrix of size nr_traces x K.
%
%   See also compute_template_pca_v2.

%% Check and initialize stuff
if isempty(W)
    error('W is empty');
end
xmm = xmm(:)';
m = size(data_in, 1);

%% Process the data
data_out = (data_in - ones(m,1)*xmm) * W;

end

