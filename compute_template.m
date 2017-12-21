function [tmiu, tsigma] = compute_template(X, mle)
%COMPUTE_TEMPLATE Compute template parameters
%   [tmiu, tsigma] = COMPUTE_TEMPLATE(X, mle)
%   computes the templates (mean, covariance matrix) for the given data.
%
%   X should be a 3-dimensional matrix of size
%   nr_samples x nr_interest_points x nr_groups containing the data for
%   which templates should be computed. The templates (tmiu, tsigma) are
%   computed for each group. That is, for group k the templates will be
%   computed from X(:,:,k).
%
%   mle is an optional integer parameter that specifies if the templates
%   should use the maximum likelihood estimate (pass mle non-zero) or the
%   unbiased covariance estimator (ignore mle or pass zero). The maximum
%   likelihood estimate for the covariance matrix divides by nr_samples
%   while the unbiased estimator divides by (nr_samples-1).
%
%   This method assumes that any pre-processing has been done and all the
%   given data in X will be used to compute the templates. Use a function
%   such as compute_features to extract only a subset of interest points
%   for each sample trace.
%
%   This method returns the following:
%   - tmiu: a matrix of size nr_groups x nr_interest_points, containing the
%     mean vector of each group.
%   - tsigma: a matrix of size nr_interest_points x nr_interest_points x nr_groups,
%     having the covariance matrix of each group.


%% Initialize and check parameters
nr_samples = size(X,1);
nr_interest_points = size(X, 2);
nr_groups = size(X, 3);
tmiu = zeros(nr_groups, nr_interest_points);
tsigma = zeros(nr_interest_points, nr_interest_points, nr_groups);
mct = 1/nr_samples;
if nargin < 2 || isempty(mle) || mle == 0
    sct = (1/(nr_samples-1));
else
    sct = mct;
end

%% Compute the templates for each group
for k=1:nr_groups
    x = X(:,:,k);
    tmiu(k,:) = mct*(ones(1,nr_samples)*x);
    xm = x - ones(nr_samples,1)*tmiu(k,:);
    tsigma(:,:,k) = sct*(xm'*xm);
end

end