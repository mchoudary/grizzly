function [metadata] = read_metadata(fname)
%READ_METADATA Reads metadata from a given data file.
%   [metadata] = READ_METADATA(fname)
%   reads the metadata from a given data file and returns this as a
%   structure.
%
%   fname should be the filename (including path) of the file containing
%   the data to be read. See 'write_leakage_data' for an example of data
%   file structure.
%
%   This method will return a metadata structure, containing information
%   about the data file. The only mandatory and common value is the
%   'format' field, which should be used to determine the kind of data
%   available in the file and the other fields that are returned. See below
%   the exact fields for each type of supported format.
%
%   You can use the returned metadata with get_mmap to obtain a memmapfile
%   object to access the data.
%
%   See also write_leakage_data_load, get_mmap, swapbytes.


%% Initialize and check parameters
fid = fopen(fname, 'rb');
if fid < 0
    error('Could not open given file');
end

%% Read format
fs = fread(fid, 1, 'uint8');
metadata.format = char(fread(fid, fs, '*uchar')');
for i=1:(7-fs)
    fread(fid, 1, 'uint8'); % ignore padded zeros
end

%% Read rest of data based on format
if strcmp(metadata.format, 'raw')
    machinefmt = 'l';
    metadata.machinefmt = machinefmt;
    metadata.nr_traces = fread(fid, 1, 'uint64', 0, machinefmt);
    metadata.nr_points = fread(fid, 1, 'uint64', 0, machinefmt);
    metadata.nr_groups = fread(fid, 1, 'uint64', 0, machinefmt);
    ps = fread(fid, 1, 'uint8', 0, machinefmt);
    metadata.precision = char(fread(fid, ps, '*uchar', 0, machinefmt)');
    for i=1:(7-ps)
        fread(fid, 1, 'uint8', 0, machinefmt); % ignore padded zeros
    end
    metadata.offset = 40;
elseif strcmp(metadata.format, 'tpl')
    machinefmt = 'l';
    metadata.machinefmt = machinefmt;
    metadata.np = fread(fid, 1, 'uint64', 0, machinefmt);
    metadata.nr_points = fread(fid, 1, 'uint64', 0, machinefmt);
    metadata.nr_groups = fread(fid, 1, 'uint64', 0, machinefmt);
    ps = fread(fid, 1, 'uint8', 0, machinefmt);
    metadata.precision = char(fread(fid, ps, '*uchar', 0, machinefmt)');
    for i=1:(7-ps)
        fread(fid, 1, 'uint8', 0, machinefmt); % ignore padded zeros
    end
    metadata.ni = fread(fid, metadata.np, metadata.precision);
    metadata.offset = 40 + (8*metadata.np);
elseif strcmp(metadata.format, 'mat2')
    machinefmt = 'l';
    metadata.machinefmt = machinefmt;
    metadata.nr_rows = fread(fid, 1, 'uint64', 0, machinefmt);
    metadata.nr_cols = fread(fid, 1, 'uint64', 0, machinefmt);
    ps = fread(fid, 1, 'uint8', 0, machinefmt);
    metadata.precision = char(fread(fid, ps, '*uchar', 0, machinefmt)');
    for i=1:(7-ps)
        fread(fid, 1, 'uint8', 0, machinefmt); % ignore padded zeros
    end
    metadata.offset = 32;
elseif strcmp(metadata.format, 'rawe2')
    machinefmt = 'l';
    metadata.machinefmt = machinefmt;
    metadata.nr_trials = fread(fid, 1, 'uint64', 0, machinefmt);
    metadata.nr_groups = fread(fid, 1, 'uint64', 0, machinefmt);
    metadata.nr_points = fread(fid, 1, 'uint64', 0, machinefmt);
    s_xfmt = fread(fid, 1, 'uint8', 0, machinefmt);
    metadata.xfmt = char(fread(fid, s_xfmt, '*uchar', 0, machinefmt)');
    for i=1:(7-s_xfmt)
        fread(fid, 1, 'uint8', 0, machinefmt); % ignore padded zeros
    end
    metadata.samplingrate = fread(fid, 1, 'double', 0, machinefmt);
    metadata.fclock = fread(fid, 1, 'double', 0, machinefmt);
    metadata.tscale = fread(fid, 1, 'double', 0, machinefmt);
    metadata.toffset = fread(fid, 1, 'double', 0, machinefmt);
    metadata.vscale = fread(fid, 1, 'double', 0, machinefmt);
    metadata.voffset = fread(fid, 1, 'double', 0, machinefmt);
    metadata.rvalue = fread(fid, 1, 'double', 0, machinefmt);
    metadata.dccoupling = fread(fid, 1, 'int64', 0, machinefmt);
    metadata.nr_bytes = fread(fid, 1, 'uint64', 0, machinefmt);
    s_bfmt = fread(fid, 1, 'uint8', 0, machinefmt);
    metadata.bfmt = char(fread(fid, s_bfmt, '*uchar', 0, machinefmt)');
    for i=1:(7-s_bfmt)
        fread(fid, 1, 'uint8', 0, machinefmt); % ignore padded zeros
    end
    metadata.address = fread(fid, 1, 'uint64', 0, machinefmt);
    s_rifmt = fread(fid, 1, 'uint8', 0, machinefmt);
    metadata.rifmt = char(fread(fid, s_rifmt, '*uchar', 0, machinefmt)');
    for i=1:(7-s_rifmt)
        fread(fid, 1, 'uint8', 0, machinefmt); % ignore padded zeros
    end
    metadata.ridxoffset = 136;
    ribs = get_bytes_class(metadata.rifmt);
    metadata.xoffset = metadata.ridxoffset + ...
        metadata.nr_bytes*ribs;
    xbs = get_bytes_class(metadata.xfmt);
    metadata.boffset = metadata.xoffset + ...
        metadata.nr_trials*metadata.nr_points*xbs;
    rbs = get_bytes_class(metadata.bfmt);
    metadata.roffset = metadata.boffset + ...
        metadata.nr_trials*metadata.nr_bytes*rbs;
else
    fprintf('Unknown format\n');
end

%% Close file
fclose(fid);

end
