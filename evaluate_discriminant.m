function [d] = evaluate_discriminant(data, groups, ...
                                     miu, sinv, slogdet, prior, varargin)
%EVALUATE_DISCRIMINANT Evaluates a discrimnant score
%   [d] = EVALUATE_DISCRIMINANT(data, groups, ...
%                               miu, sinv, slogdet, prior, varargin)
%   evaluates a discriminant score for the given samples and parameters.
%
%   The idea behind using a discriminant score is to compare the
%   discriminants of some given data for many different groups
%   for which miu and sinv are known, and then decide to which group
%   belongs X by checking for which group the discriminant score is
%   highest.
%
%   data should be a matrix of size nr_samples x nr_points containing
%   the test samples for which to compute the discriminant scores.
%   All the samples in data will be used to compute a discriminant score
%   derived from the joint likelihood of the samples. This discriminant
%   score will be computed for each group requested, as explained below.
%
%   groups should be a vector of length nr_test_groups.
%   Regardless of the number of samples (nr_samples) in the data matrix,
%   groups should specify the indexes for miu for which the discriminant
%   will be computed.
%
%   miu should be a matrix of size nr_groups x nr_points
%   containing the mean vectors (one per row) for all groups.
%   
%   sinv should be a matrix of nr_points x nr_points x nr_groups having the
%   inverse of the unbiased covariance matrix for each group, or a matrix
%   of size nr_points x nr_points if the same covariance should be used
%   for all groups (e.g. when using the pooled covariance). Note that if in
%   your application you can use a pooled covariance estimator this will
%   improve much both the results and the computation time. That is because
%   in that case this method will use a linear discriminant (depends only
%   linearly on x) rather than a quadratic one.
%   Pass [] (empty) if only the mean miu should be used as template. This
%   is generally not recommended but is very useful if the data has been
%   preprocessed to have an identity covariance matrix, such as in the case
%   of using Fisher's LDA. In this case the computation will be much faster
%   and more accurate.
%
%   When sinv is given as a different matrix for each group you must
%   also provide slogdet as a vector of length nr_groups having the
%   log-determinant of the covariance matrices for which the inverse was
%   given in sinv. In all other cases pass empty, as this will be ignored.
%
%   miu (and sinv/logdet if not using a common covariance) should be
%   provided for all the classes that are specified by the parameter groups
%   (see above).
%
%   prior is a vector of length nr_groups having the prior probability of
%   each group. This should be such that sum(prior) = 1.
%   Pass [] (empty) to use equal probability.
%
%   This method returns a vector d of size nr_test_groups x 1, with the
%   discriminant scores of the data for each group specified by the
%   groups vector.
%
%   Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)

%% Check and initilise parameters
groups = groups(:);
ns = size(sinv, 3);

%% Compute the output
if ns > 1
    d = compute_discriminant(data, miu(groups,:), ...
                             sinv(:,:,groups), slogdet(groups), prior);
else
    d = compute_discriminant(data, miu(groups,:), ...
                             sinv, slogdet, prior);
end

end