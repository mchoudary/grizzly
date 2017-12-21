function [data_out] = prepare_data_template(data_in, interest_points, ...
                                            varargin)
%PREPARE_DATA_TEMPLATE Prepares data for leakage analysis
%   [data_out] = PREPARE_DATA_TEMPLATE(data_in, interest_points, varargin)
%   transforms the input trace(s) into the required data for leakage
%   analysis by means of the given parameters.
%
%   This function simply selects the interesting points from the traces.
%   
%   interest_points is a vector of length nr_interest_points, containing
%   the positions in the traces containing the interesting points that must
%   be processed (e.g. for templates).
%
%   This method should be called after reading each trace from the file and
%   before passing the data to get_leakage_info.

%% Check and initialize stuff
if isempty(interest_points)
    error('interest_points is empty');
end

%% Process the data
data_out = data_in(:,interest_points);

