function [curves] = get_signal_strength_ssp(M, B, W, nr_traces)
%GET_SIGNAL_STRENGTH_SSP Returns signal strength estimates
%   [curves] = GET_SIGNAL_STRENGTH_SSP(M, B, W, nr_traces)
%   returns the signal strength estimates (DOM, SOSD, etc.) for side
%   channel analysis, based on the matrices M, B and W that contain the
%   group means, treatment sum of squares and crossproducts matrix and
%   residual sum of squares and crossproducts matrix respectively.
%
%   M is the group means matrix, of size nr_groups x nr_points.
%
%   B is the treatment sum of squares and crossproducts matrix of size
%   nr_points x nr_points. You can pass [], in which case the 'ftest'
%   vector will not be computed.
%
%   W is the residual sum of squares and crossproducts matrix of size
%   nr_points x nr_points. W = C * (nr_groups*(nr_traces-1)), where C is
%   the common (pooled) covariance of all groups.
%
%   nr_traces is the number of traces per group that were used to compute
%   M, B and W.
%
%   You can obtain M, B and W from compute_ssp_e2_mmap for example.
%   
%   This method returns a structure, with the following data:
%   - 'dom': the signal vector based on Difference of Means (DOM). See the
%   "Template Attacks" paper by Chari et al.
%   - 'sosd': the signal vector based on Sum of Squared Differences (SOSD).
%   See the paper "Templates vs Stochastic Methods".
%   - 'snr': the signal vector based on the signal to noise ratio. See the
%   book "Power Analysis Attacks: Revealing the Secrets of Smartcards".
%   - 'sost': the sost vector based on the SOST method. See "Templates vs
%   Stochastic Methods" by Gierlichs et al.
%   - 'ftest': a vector based on the F-test from ANOVA. See the "Applied
%   Multivariate Statistical Analysis" book by Johnson and Wichern.
%   Computed only if B is not empty.
%
%   See also compute_ssp_e2_mmap, get_mmap, read_metadata.


%% Check and initialize stuff
nr_groups = size(M, 1);
nr_points = size(M, 2);
if ~isempty(B) && (size(B,1) ~= nr_points || size(B,2) ~= nr_points)
    error('Incorrect size of B');
end
if size(W,1) ~= nr_points || size(W,2) ~= nr_points
    error('Incorrect size of W');
end

%% Compute the difference of means curve (DOM)
curves.dom = zeros(1, nr_points);
for i=1:nr_groups-1
    for j=i+1:nr_groups
        ds = (M(i,:) - M(j,:));
        curves.dom = curves.dom + abs(ds);
    end
end

%% Compute the difference of means squared (SOSD)
curves.sosd = zeros(1, nr_points);
for i=1:nr_groups-1
    for j=i+1:nr_groups
        ds = (M(i,:) - M(j,:)).^2;
        curves.sosd = curves.sosd + ds;
    end
end

%% Compute the SNR curve
V = diag(W) / (nr_groups*(nr_traces-1));
V = V(:)';
M_var = var(M, 0, 1);
curves.snr = M_var ./ V;

%% Compute the SOST curve
curves.sost = zeros(1, nr_points);
for i=1:nr_groups-1
    for j=i+1:nr_groups
        curves.sost = curves.sost + ...
            ((M(i,:)-M(j,:)).^2) ./ V;
    end
end

%% Compute the F-test curve
if ~isempty(B)
    curves.ftest = zeros(1, nr_points);
    for k=1:nr_points
        curves.ftest(k) = (B(k,k)/(nr_groups-1)) / (W(k,k)/(nr_groups*(nr_traces-1)));
    end
end

end

