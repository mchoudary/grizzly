function d = compute_discriminant(X, miu, sinv, slogdet, prior)
%COMPUTE_DISCRIMINANT Computes a discriminant score
%   [d] = COMPUTE_DISCRIMINANT(X, miu, sinv)
%   computes the discriminant score of X based on miu and sinv.
%
%   X should be a matrix of size nr_samples x nr_points containing
%   the test samples for which to compute the discriminant score.
%   All the samples in X will be used to compute a discriminant score
%   derived from the joint likelihood of the samples. This discriminant
%   score will be computed for all groups for which miu and sinv are given.
%
%   The idea behind using a discriminant score is to compare the
%   discriminants of some given data (such as X) for many different groups
%   for which miu and sinv are known, and then decide to which group
%   belongs X by checking for which group the discriminant score is
%   highest.
%   
%   X should be a matrix of size N x D, having N D-dimensional observations
%   for which the discriminant score will be computed.
%
%   miu must be a matrix of size nr_groups x D having the mean vector of
%   the nr_groups groups for which the discriminants should be computed.
%
%   sinv should be a matrix of nr_points x nr_points x D having the
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
%   prior is a vector of length nr_groups having the prior probability of
%   each group. This should be such that sum(prior) = 1.
%   Pass [] (empty) to use equal probability.
%
%   See the "Applied Multivariate Statistical Analysis" book for details.
%
%   Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)


%% Check and initialise parameters
if nargin < 5
    error('Insufficient arguments');
end
N = size(X, 1);
D = size(X, 2);
nr_groups = size(miu, 1);
if size(miu, 2) ~= D
    error('Incompatible size of miu with X');
end
if ~isempty(sinv) && (size(sinv, 1) ~= D || size(sinv, 2) ~= D)
    error('Bad dimension sinv');
end
if ~isempty(sinv) && size(sinv,3) > 1 && size(sinv, 3) ~= nr_groups
    error('Incompatible third dimension for sinv');
end
if ~isempty(slogdet)
    slogdet = slogdet(:);
end
if ~isempty(sinv) && size(sinv,3) > 1 && ...
        (isempty(slogdet) || length(slogdet) ~= nr_groups)
    error('Incorrect vector slogdet');
end
if ~isempty(prior)
    prior = prior(:);
    if length(prior) ~= nr_groups
        error('Incorrect length of prior');
    end
end
d = zeros(nr_groups, 1);

%% Compute discriminant from given parameters
if ~isempty(sinv) && size(sinv,3) > 1
    %% Compute quadratic discriminant when sinv is given for each group
    % This is the slowest case
    ct1 = -N/2;
    ct2 = -1/2;
    for k = 1:nr_groups
        dsum = 0;
        for j = 1:N
            x = X(j,:) - miu(k,:);
            dsum = dsum + x*sinv(:,:,k)*x';
        end
        d(k) = ct1*slogdet(k) + ct2*dsum;
    end
elseif ~isempty(sinv) && size(sinv,3) == 1
    %% Compute linear discriminant using common covariance
    % Much faster than previous method since due to linearity we can first
    % add the test traces and then make only one multiplication
    ct = -N/2;
    xs = ones(1,N) * X;
    for k = 1:nr_groups
        d(k) = miu(k,:)*sinv*xs' + ct*(miu(k,:)*sinv*miu(k,:)');        
    end
elseif isempty(sinv)
    %% Compute linear discriminant using with identity covariance
    % Fastest method, use only when actually having an identity covariance
    ct = -N/2;
    xs = ones(1,N) * X;
    for k = 1:nr_groups
        d(k) = miu(k,:)*xs' + ct*(miu(k,:)*miu(k,:)');        
    end
end

%% Add priors to each discriminant if given
if ~isempty(prior)
    d = d + N*log(prior);
end

end



