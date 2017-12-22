function [signal] = get_signal_strength_coef(X, coef, base)
%GET_SIGNAL_STRENGTH_COEF Returns signal strength estimates
%   [signal] = GET_SIGNAL_STRENGTH_COEF(X, coef)
%   returns the signal strength estimates for stochastic attacks, based on
%   the leakage matrix X and estimated base coefficients coef.
%
%   X is the leakge matrix of size nr_traces x nr_points.
%
%   coef is the matrix of base coeficients for each point, of size
%   nr_bases x nr_points.
%
%   base is a string specifying which base was used to compute the
%   coefficients, since this is important. Currently supported bases are:
%   - 'F9': constant power consumption plus individual 8 bits.
%   - 'F17': constant power consumption plus individual 16 bits.
%   - 'F17xor'
%   - 'F17tran'
%
%   This method returns a structure, with the following data:
%   - 'bnorm': the signal strength estimate based on the squared euclidian
%   norm of the base vecctors. See the Stochastic paper, CHES 2005.
%   - 'bnorm_std': uses the norm of base vectors and variance of each
%   sample point.
%
%   See also get_selection, get_signal_strength_ssp.


%% Check and initialize stuff
nr_points = size(X, 2);

%% Compute for selected base
if strcmp(base, 'F9') || strcmp(base, 'F17') || strcmp(base, 'F17xor') || strcmp(base, 'F17tran')
    %% Compute bnorm
    signal.bnorm = zeros(1, nr_points);
    for j=1:nr_points
        signal.bnorm(j) = coef(2:end,j)'*coef(2:end,j);
    end

    %% Compute bnorm_std
    vx = var(X, 0, 1);
    signal.bnorm_std = zeros(1, nr_points);
    for j=1:nr_points
        signal.bnorm_std(j) = signal.bnorm(j) ./ vx(j);
    end
else
    error('Unsupported base: %s', base);
end

end

