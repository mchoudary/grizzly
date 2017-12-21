function [mmap, metadata] = get_mmap(fname)
%GET_MMAP Returns a memmapfile object for the given data
%   [mmap, metadata] = GET_MMAP(fname)
%   returns a memmapfile object for the given data file.
%
%   fname should be the filename (including path) of the file containing
%   the data to be read. See 'write_mat_data' for an example of data
%   file structure.
%
%   Currently supported formats are:
%   - raw
%   - tpl
%   - mat2
%   - rawe2
%
%   This method retunrs two objects:
%   - mmap: the actual memmapfile object for the given data.
%   - metadata: the metadata structure associated with the data.
%
%   See also read_metadata, write_leakage_data_load, swapbytes.


%% Initialize and check parameters

%% Load metadata
metadata = read_metadata(fname);

%% Create memorymap object
if strcmp(metadata.format, 'raw')
    mmap = memmapfile(fname, 'Offset', metadata.offset, ...
                   'Repeat', metadata.nr_groups, ...
                   'Format',{metadata.precision, ...
                    [metadata.nr_traces metadata.nr_points],'X'});
elseif strcmp(metadata.format, 'tpl')
    mmap = memmapfile(fname, 'Offset', metadata.offset, ...
                      'Repeat', metadata.nr_groups, ...
                      'Format', { ...
                      metadata.precision, [1, metadata.nr_points], 'tmiu'; ...
                      metadata.precision, ...
                      [metadata.nr_points, metadata.nr_points],'tsigma'});    
elseif strcmp(metadata.format, 'mat2')
    mmap = memmapfile(fname, 'Offset', metadata.offset, ...
                      'Format', { ...
                      metadata.precision, ...
                      [metadata.nr_cols, metadata.nr_rows],'data'});
elseif strcmp(metadata.format, 'rawe2')
    mmap = memmapfile(...
        fname, ...
        'Offset', metadata.ridxoffset, ...
        'Format', { ...
          metadata.rifmt, [metadata.nr_bytes, 1],'rindex';
          metadata.xfmt, [metadata.nr_points, metadata.nr_trials],'X';
          metadata.bfmt, [metadata.nr_bytes, metadata.nr_trials], 'B';
          metadata.bfmt, [metadata.nr_bytes, metadata.nr_trials], 'R'});
else
    error('Unknown metadata format');
end

end
