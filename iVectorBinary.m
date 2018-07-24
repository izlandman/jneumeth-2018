% given the file sizes are increasing, begin reading and writing as binary
% files when handling I-Vectors

% pass in an array of I-Vectors (IVECTORS) and the depth of the I-Vectors
% when writing. Nothing needs to be passed when reading as the first
% element of the binary file will be the depth of the I-Vectors.
function i_vectors = iVectorBinary(varargin)

% if two arguments, write to file
% if one arguments, read from file
if( nargin == 2 )
    i_vectors = varargin{2};
    file_name = varargin{1};
    [depth,~] = size(i_vectors);
    i_vectors = i_vectors(:);
    i_vectors = [depth; i_vectors];
    fid = fopen(file_name,'w');
    fwrite(fid,i_vectors,'double');
    fclose(fid);
elseif( nargin == 1 )
    file_name = varargin{1};
    try
        fid = fopen(file_name,'r');
        i_vectors = fread(fid,'double');
        fclose(fid);
        depth = i_vectors(1);
        n = numel(i_vectors) - 1;
        i_vectors = reshape(i_vectors(2:end),depth,n/depth);
    catch
        fprintf('Unable to read file: %s\n', file_name);
        i_vectors = -1;
    end
else
    fprintf('iVectorBinary has failed. Check # of arguments.\n');
end