function [g] = get_ge_from_success_info(sinfo, nr_traces_vec)
%GET_GE_FROM_SUCCESS_INFO Returns guessing entropy data.
%   [g] = GET_GE_FROM_SUCCESS_INFO(sinfo, nr_traces_vec)
%   returns guessing entropy data from the given success info structure.
%
%   sinfo must be a structure as returned by get_success_info_like,
%   containing depth data for both ".avg" and ".joint" scores.
%
%   nr_traces_vec should be a vector of length nr_test_groups, with the
%   number of attack traces used in the experiments for which sinfo is
%   provided.
%
%   This method returns a structure ge, with 2 fields:
%   - g.avg: guessing entropy data from the ".avg" scores in sinfo.
%   - g.joint: guessing entropy data from the ".joint" scores in sinfo.
%   Both g.avg and g.joint are vectors of length nr_test_groups.
%
%   Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)

%% Check and initialize stuff
d.avg = sinfo.depth.avg;
d.joint = sinfo.depth.joint;
nr_test_groups = length(nr_traces_vec);
g.avg = zeros(nr_test_groups, 1);
g.joint = zeros(nr_test_groups, 1);

%% Compute the mean guessing entropy for all data
for k=1:nr_test_groups
    nr_iter = size(d.avg.group1, 2);
    avg = 0;
    joint = 0;
    for j = 1:nr_iter
        avg = avg + log2(mean(d.avg.(['group' num2str(k)])(:,j)));
        joint = joint + log2(mean(d.joint.(['group' num2str(k)])(:,j)));        
    end
    g.avg(k) = avg / nr_iter;
    g.joint(k) = joint / nr_iter;
end

end