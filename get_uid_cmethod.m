function [uid] = get_uid_cmethod(cmethod, cparams)
%GET_UID_CMETHOD Returns a Unique ID (UID) for a compression method/params
%   [uid] = GET_UID_CMETHOD(cmethod, cparams)
%   returns a unique ID (UID) for a particular compression method and
%   parameters.
%
%   The unique ID for each cmethod/cparams combination can be used for
%   example to use the same set plotting parameters.
%
%   Current UIDs:
%   1: LDA, m=4
%   2: PCA, m=4
%   3: sample select, 1ppc
%   4: sample select, 3ppc
%   5: sample select, 20ppm
%   6: sample select, allap
%   7: LDA, m=40
%   8: LDA, m=100 or m=5
%   9: PCA, m=40
%   10: PCA, m=100 or m=5
%   11: LDA, m=6
%   12: PCA, m=6
%   13: LDA, m=3
%
%   If the cmethod/cparams combination is unknown this method will return 0;

%% Return UID
uid = 0;
if strcmp(cmethod, 'LDA')
    if cparams.lda_dimensions == 4
        uid = 1;
    elseif cparams.lda_dimensions == 40
        uid = 7;
    elseif cparams.lda_dimensions == 100 || cparams.lda_dimensions == 5
        uid = 8;
    elseif cparams.lda_dimensions == 6
        uid = 11;
    elseif cparams.lda_dimensions == 3
        uid = 13;
    end
elseif strcmp(cmethod, 'PCA')
    if cparams.pca_dimensions == 4
        uid = 2;
    elseif cparams.pca_dimensions == 40
        uid = 9; 
    elseif cparams.pca_dimensions == 100 || cparams.pca_dimensions == 5
        uid = 10;    
    elseif cparams.pca_dimensions == 6
        uid = 12;
    end
elseif strcmp(cmethod, 'sample')
    if strcmp(cparams.sel, '1ppc')
        uid = 3;
    elseif strcmp(cparams.sel, '3ppc')
        uid = 4;
    elseif strcmp(cparams.sel, '20ppc')
        uid = 5;
    elseif strcmp(cparams.sel, 'allap')
        uid = 6;
    end
end

end