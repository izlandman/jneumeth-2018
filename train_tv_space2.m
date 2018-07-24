function T = train_tv_space2(dataList, ubmFilename, tv_dim, niter, ...
    nworkers, tvFilename)
% uses statistics in dataLits to train the i-vector extractor with tv_dim
% factors and niter EM iterations. The training process can be parallelized
% via parfor with nworkers. The output can be optionally saved in tvFilename.
%
% Technically, assuming a factor analysis (FA) model of the from:
%
%           M = m + T . x
%
% for mean supervectors M, the code computes the maximum likelihood
% estimate (MLE)of the factor loading matrix T (aka the total variability
% subspace). Here, M is the adapted mean supervector, m is the UBM mean
% supervector, and x~N(0,I) is a vector of total factors (aka i-vector).
%
% Inputs:
%   - dataList    : ASCII file containing stats file names (1 file per line)
%                   or a cell array of concatenated stats (i.e., [N; F])
%   - ubmFilename : UBM file name or a structure with UBM hyperparameters
%   - tv_dim      : dimensionality of the total variability subspace
%   - niter       : number of EM iterations for total subspace learning
%   - nworkers    : number of parallel workers
%   - tvFilename  : output total variability matrix file name (optional)
%
% Outputs:
%   - T 		  : total variability subspace matrix
%
% References:
%   [1] D. Matrouf, N. Scheffer, B. Fauve, J.-F. Bonastre, "A straightforward
%       and efficient implementation of the factor analysis model for speaker
%       verification," in Proc. INTERSPEECH, Antwerp, Belgium, Aug. 2007,
%       pp. 1242-1245.
%   [2] P. Kenny, "A small footprint i-vector extractor," in Proc. Odyssey,
%       The Speaker and Language Recognition Workshop, Singapore, Jun. 2012.
%   [3] N. Dehak, P. Kenny, R. Dehak, P. Dumouchel, and P. Ouellet, "Front-end
%       factor analysis for speaker verification," IEEE TASLP, vol. 19, pp.
%       788-798, May 2011.
%
%
% Omid Sadjadi <s.omid.sadjadi@gmail.com>
% Microsoft Research, Conversational Systems Research Center

if ischar(tv_dim), tv_dim = str2double(tv_dim); end
if ischar(niter), niter = str2double(niter); end
if ischar(nworkers), nworkers = str2double(nworkers); end

if ischar(ubmFilename)
    tmp  = load(ubmFilename);
    ubm  = tmp.gmm;
elseif isstruct(ubmFilename)
    ubm = ubmFilename;
else
    error('Oops! ubmFilename should be either a string or a structure!');
end
[ndim, nmix] = size(ubm.mu);
S = reshape(ubm.sigma, ndim * nmix, 1);

[N, F] = load_data(dataList, ndim, nmix);
if iscell(dataList), clear dataList; end

% fprintf('\n\nRandomly initializing T matrix ...\n\n');
% suggested in jfa cookbook
% T = randn(tv_dim, ndim * nmix) * sum(S) * 0.001;
% or how about we let it scale each column by S?
T = bsxfun( @times, randn(tv_dim, ndim * nmix), S');
T_best = T;
% fprintf('Re-estimating the total subspace with %d factors ...\n', tv_dim);
min_diff = [100,100,100];
min_update = 5;
iter = 1;
max_iter = 50;
while iter <= niter
    % fprintf('EM iter#: %d \t', iter);
    tim = tic;
    
    [LU, RU] = expectation_tv(T, N, F, S, tv_dim, nmix, ndim, nworkers);
    
    % LU and RU seem to have a problem with accuracy within Matlab. Often
    % times values will come back as very small -/+e-04/5/6 when the
    % marojirty of the matrix is not as such. This may make Matlab upset
    % and cause issues with resolving inversions.
    % RU = matrixAdjust(RU,100000);
    % LU = matrixAdjust(LU,1000);
    T_new = maximization_tv(LU, RU, ndim, nmix);
    % how much has T changed?
    PD = abs( (T_new-T)./T );
    column_diff = mean(PD)*100;
    row_diff = mean(PD,2)*100;
    full_mean = mean(row_diff);
    row_std = std(row_diff);
    col_std = std(column_diff);
    T = T_new;
    tim = toc(tim);
    % fprintf('[elaps = %.2f s]\n', tim);
    % fprintf('Percent Diff [mean | row sigma col sigma] %3.2f | %3.2f %3.2f\n', ...
    %    full_mean, row_std, col_std);
    if( min_diff(1) > full_mean )
        min_diff = [full_mean,row_std,col_std];
        T_best = T;
    end
    if((iter == niter) && (min_diff(1) > min_update) && (niter < max_iter))
        niter = niter + 1;
    end
    iter = iter + 1;
    % check if T is NaN or Inf
        IN = isnan(T);
        II = isinf(T);
        new_T = randn(tv_dim, ndim * nmix) * sum(S) * 0.001;
        if( sum(sum(IN)) > 0 )
            fprintf('T has NaN values!\n');
            % reset bad values to random state
            T(IN) = new_T(IN);
        end
        if( sum(sum(II)) > 0 )
            fprintf('T has Inf values!\n');
            % reset bad values to random state
            T(II) = new_T(II);
        end
end

fprintf('Best Percent Diff of mixture %d: %3.2f %3.2f %3.2f\n', ...
    nmix, min_diff);
T = T_best;

if ( nargin == 6 )
    fprintf('\nSaving T matrix to file %s\n', tvFilename);
    % create the path if it does not exist and save the file
    path = fileparts(tvFilename);
    if ( exist(path, 'dir')~=7 && ~isempty(path) ), mkdir(path); end
    save(tvFilename, 'T');
end

end

function [N, F] = load_data(datalist, ndim, nmix)
% load all data into memory
% !!! WARNING !!!
% the original version of this code called for the N & F matrices to be
% setup as 'single' precision. doing saves space in memory. however, this
% causes problems with the matrix accuracy when estimating the T matrix.
% changing them to 'double' which is the default in Matlab appears to
% correct matrix errors.
if ischar(datalist),
    fid = fopen(datalist, 'rt');
    filenames = textscan(fid, '%s');
    fclose(fid);
    filenames = filenames{1};
    nfiles = size(filenames, 1);
    N = zeros(nfiles, nmix, 'double');
    F = zeros(nfiles, ndim * nmix, 'double');
    for file = 1 : nfiles,
        tmp = load(filenames{file});
        N(file, :) = tmp.N;
        F(file, :) = tmp.F;
    end
else
    nfiles = length(datalist);
    N = zeros(nfiles, nmix, 'double');
    F = zeros(nfiles, ndim * nmix, 'double');
    for file = 1 : nfiles,
        N(file, :) = datalist{file}(1:nmix);
        F(file, :) = datalist{file}(nmix + 1 : end);
    end
end

end

function [LU, RU] = expectation_tv(T, N, F, S, tv_dim, nmix, ndim, nworkers)
% compute the posterior means and covariance matrices of the factors
% or latent variables
idx_sv = reshape(repmat(1 : nmix, ndim, 1), ndim * nmix, 1);
nfiles = size(F, 1);

LU = cell(nmix, 1);
LU(:) = {zeros(tv_dim)};

RU = zeros(tv_dim, nmix * ndim);
I = eye(tv_dim);
T_invS =  bsxfun(@rdivide, T, S');

parts = 250; % modify this based on your resources
nbatch = floor( nfiles/parts + 0.99999 );

% attempt to catch the NaN culprit
if( sum(sum(isnan(T))) > 0 || sum(sum(isinf(T))) > 0 )
    fprintf('Something is wrong with expecation_tv setup.\n');
end

for batch = 1 : nbatch
    start = 1 + ( batch - 1 ) * parts;
    fin = min(batch * parts, nfiles);
    len = fin - start + 1;
    index = start : fin;
    N1 = N(index, :);
    F1 = F(index, :);
    % check F1 for zero values
    F1(F1==0) = eps();
    Ex = zeros(tv_dim, len);
    Exx = zeros(tv_dim, tv_dim, len);
    parfor (ix = 1 : len, nworkers)
        %     for ix = 1 : len,
        L = I +  bsxfun(@times, T_invS, N1(ix, idx_sv)) * T';
        % pinv keeps returning errors that it is operating on NaN or
        % zero. To correct this check for NaNs and turn them into very
        % small, nearly zero, values using eps. check for infs and turn
        % them into
        
        % check if L is NaN or Inf
%         IN = isnan(L);
%         II = isinf(L);
%         
%         if( sum(sum(IN)) > 0 )
%             fprintf('L has NaN values!\n');
%         end
%         if( sum(sum(II)) > 0 )
%             fprintf('L has Inf values!\n');
%         end
        Cxx = pinv(L); % this is the posterior covariance Cov(x,x)
        
        B = T_invS * F1(ix, :)';
        Ex(:, ix) = Cxx * B; % this is the posterior mean E[x]
        Exx(:, :, ix) = Cxx + Ex(:, ix) * Ex(:, ix)';
    end
    % it appears that sometimes near zero values in Ex ( 4*10E-3) are
    % triggering values of RU to be orders of magnitude larger than they
    % should be which is leading to problems further on in the code.
    RU = RU + Ex * F1;
    parfor (mix = 1 : nmix, nworkers)
        %     for mix = 1 : nmix,
        tmp = bsxfun(@times, Exx, reshape(N1(:, mix),[1 1 len]));
        LU{mix} = LU{mix} + sum(tmp, 3);
        % fprintf('Check rank of LU: %d\n', rank(LU{mix}));
    end
end

end

function RU = maximization_tv(LU, RU, ndim, nmix)
% ML re-estimation of the total subspace matrix or the factor loading
% matrix
% first_flag = 0;
% second_flag = 0;
for mix = 1 : nmix
    idx = ( mix - 1 ) * ndim + 1 : mix * ndim;
    RU(:,idx) = LU{mix} \ RU(:, idx);
    % this seems to return poorly conditioned matrices, catch and 'correct'
    % these when this happens
    %     [RU(:,idx), inital_error, final_error] = matrixMassage(lastwarn, ...
    %         RU_temp, LU, RU, idx, mix);
    %     first_flag = first_flag + inital_error;
    %     second_flag = second_flag + final_error;
end

% fprintf('There were %d RCOND warnings, corrected to %d RCOND warnings.\n', ...
%     first_flag, second_flag);
end