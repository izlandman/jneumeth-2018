% instead of running this for each mixture, run it once and collect the
% mixtures as they are produced. this means that NMIX is now a vector built
% from all of the mixutres that will be SAVED. It must be in order, but
% does not need to be complete as the largest value will be used to create
% the incremental steps.

function ubms = allUBMs(data_list, nmix, final_niter, ds_factor, ...
    folder, verbose)
% fits a nmix-component Gaussian mixture model (GMM) to data in dataList
% using niter EM iterations per binary split. The process can be
% parallelized in nworkers batches using parfor.
%
% Inputs:
%   - dataList    : ASCII file containing feature file names (1 file per
%                   line) or a cell array containing features (nDim x
%                   nFrames). Feature files must be in uncompressed
%                   HTK format.
%   - nmix        : vector of Gaussian components (must be a power of 2)
%   - final_iter  : number of EM iterations in the final split
%   - ds_factor   : feature sub-sampling factor (every ds_factor frame)
%   - nworkers    : number of parallel workers
%   - filepath    : location to save UBMs
%
% Outputs:
%   - ubms		  : a cell of each UBM requested in NMIX
%					(ubm.mu: means, ubm.sigma: covariances, ubm.w: weights)
%
%
% Omid Sadjadi <s.omid.sadjadi@gmail.com>
% Microsoft Research, Conversational Systems Research Center

[ispow2, ~] = log2(nmix);
if ( sum(ispow2==0.5) ~= numel(nmix) )
    error('NMIX must contain powers of 2.');
end

if ischar(data_list) || iscellstr(data_list)
    data_list = load_data(data_list);
end
if ~iscell(data_list)
    error('Oops! dataList is incorrect format!');
end

nfiles = length(data_list);
if( verbose )
    fprintf('\n\nInitializing the GMM hyperparameters ...\n');
end
[gm, gv] = comp_gm_gv(data_list);
gmm = gmm_init(gm, gv);

% process NMIX to find max and saving mixtures
nmix_max = max(nmix);
[~,nmix_count] = log2(nmix_max);
ubms = cell(numel(nmix),1);
ubms_iter = 1;
% gradually increase the number of iterations per binary split, but make it
% automatic incase the request is larger than some hard coded limit!
niter = ones(nmix_count,1);
niter(2:end) = floor( (2:nmix_count) / 2)*2;
% those requested to be saved will have iterations overwritten?
niter(log2(nmix) + 1) = final_niter;
a0 = tic;
for m=1:nmix_count
    gmm_old = gmm;
    % not for the last two splits!
    if ( m >= nmix_count-2 )
        ds_factor = 1;
    end
    mix = 2^(m-1);
    if( verbose )
        fprintf('\nRe-estimatGMM hyperparameters for %d components ...\n',mix);
    end
    for iter = 1 : niter(m)
        if( verbose )
            fprintf('EM iter#: %d \t', iter);
        end
        N = 0;
        F = 0;
        S = 0;
        L = 0;
        nframes = 0;
        tim = tic;
        parfor ix = 1 : nfiles
            [n, f, s, l] = ...
                expectation(data_list{ix}(:, 1:ds_factor:end), gmm);
            N = N + n;
            F = F + f;
            S = S + s;
            L = L + sum(l);
            nframes = nframes + length(l);
        end
        tim = toc(tim);
        if( verbose )
            fprintf('[llk = %.2f] \t [elaps = %.2f s]\n', L/nframes, tim);
        end
        gmm = maximization(N, F, S);
        % what happens if NaN?
        if( sum(sum(isnan(gmm.mu))) > 0 || sum(sum(isnan(gmm.sigma))))
            gmm = gmm_old;
        end     
    end
    % brief problem with NaN being written into mu, and thus sigma. This
    % made the eer function return all NaNs which isn't helpful. added this
    % to help debug, but it seemed to 'correct' the error on its own or
    % with this brief 'pause' in the code?
    if( verbose )
        fprintf('NaNs in mean: %d \n', sum(sum(isnan(gmm.mu))));
        fprintf('Nans in sigma: %d \n', sum(sum(isnan(gmm.sigma))));
    end
    if( ismember(mix,nmix) )
        ubms{ubms_iter} = gmm;
        ubms_iter = ubms_iter + 1;
    end
    % if this mixture is flagged, save it! only if folder is valid!
    if( ismember(mix,nmix) && ischar(folder) )
        file_name = ['ubm_m' num2str(mix) '_i' num2str(final_niter) ...
            '_f' num2str(ds_factor) '.mat'];
        file_name = [folder filesep file_name];
        fprintf('Saving UBM as %s\n', file_name);
        % save it!
        save(file_name,'gmm');
    end
    % this performs a binary split to go to the next mixture size!
    gmm = gmm_mixup(gmm);
end
a1 = toc(a0);
if( verbose )
    fprintf('UBMs completed in %f seconds.\n', a1);
end
end

function data = load_data(datalist)
% load all data into memory
if ~iscellstr(datalist)
    fid = fopen(datalist, 'rt');
    filenames = textscan(fid, '%s');
    fclose(fid);
    filenames = filenames{1};
else
    filenames = datalist;
end
nfiles = size(filenames, 1);
data = cell(nfiles, 1);
parfor ix = 1 : nfiles
    data{ix} = htkread(filenames{ix});
end
end

function [gm, gv] = comp_gm_gv(data)
% computes the global mean and variance of data
nframes = cellfun(@(x) size(x, 2), data, 'UniformOutput', false);
nframes = sum(cell2mat(nframes));
gm = cellfun(@(x) sum(x, 2), data, 'UniformOutput', false);
gm = sum(cell2mat(gm'), 2)/nframes;
gv = cellfun(@(x) sum(bsxfun(@minus, x, gm).^2, 2), data, ...
    'UniformOutput', false);
gv = sum(cell2mat(gv'), 2)/( nframes - 1 );
end

function gmm = gmm_init(glob_mu, glob_sigma)
% initialize the GMM hyperparameters (Mu, Sigma, and W)
gmm.mu    = glob_mu;
gmm.sigma = glob_sigma;
gmm.w     = 1;
end

function [N, F, S, llk] = expectation(data, gmm)
% compute the sufficient statistics
[post, llk] = postprob(data, gmm.mu, gmm.sigma, gmm.w(:));
N = sum(post, 2)';
F = data * post';
S = (data .* data) * post';
end

function [post, llk] = postprob(data, mu, sigma, w)
% compute the posterior probability of mixtures for each frame
post = lgmmprob(data, mu, sigma, w);
llk  = logsumexp(post, 1);
post = exp(bsxfun(@minus, post, llk));
end

function logprob = lgmmprob(data, mu, sigma, w)
% compute the log probability of observations given the GMM
ndim = size(data, 1);
s_recip = 1./sigma;
s_mu = mu.*s_recip;
% C = sum(mu.*mu./sigma) + sum(log(sigma));
C = sum(mu.*s_mu) + sum(log(sigma));
% D = (1./sigma)' * (data .* data) - 2 * (mu./sigma)' * data  + ...
%     ndim * log(2 * pi);
D = s_recip'*(data.*data);
D = D - 2*s_mu'*data;
D = D + ndim*log(2*pi);
logprob = -0.5 * (bsxfun(@plus, C',  D));
logprob = bsxfun(@plus, logprob, log(w));
% logprob = log(w) - 0.5 * (C' + D);
end

function y = logsumexp(x, dim)
% compute log(sum(exp(x),dim)) while avoiding numerical underflow
xmax = max(x, [], dim);
y    = xmax + log(sum(exp(bsxfun(@minus, x, xmax)), dim));
ind  = find(~isfinite(xmax));
if ~isempty(ind)
    y(ind) = xmax(ind);
end
end

function gmm = maximization(N, F, S)
% ML re-estimation of GMM hyperparameters which are updated from accumulators
w  = N / sum(N);
mu = bsxfun(@rdivide, F, N);
sigma = bsxfun(@rdivide, S, N) - (mu .* mu);
sigma = apply_var_floors(w, sigma, 0.1);
gmm.w = w;
gmm.mu= mu;
gmm.sigma = sigma;
end

function sigma = apply_var_floors(w, sigma, floor_const)
% set a floor on covariances based on a weighted average of component
% variances
vFloor = sigma * w' * floor_const;
sigma  = bsxfun(@max, sigma, vFloor);
% sigma = bsxfun(@plus, sigma, 1e-6 * ones(size(sigma, 1), 1));
end

function gmm = gmm_mixup(gmm)
% perform a binary split of the GMM hyperparameters
mu = gmm.mu; sigma = gmm.sigma; w = gmm.w;
[ndim, nmix] = size(sigma);
[sig_max, arg_max] = max(sigma);
eps = sparse(0 * mu);
eps(sub2ind([ndim, nmix], arg_max, 1 : nmix)) = sqrt(sig_max);
% only perturb means associated with the max std along each dim
mu = [mu - eps, mu + eps];
% mu = [mu - 0.2 * eps, mu + 0.2 * eps]; % HTK style
sigma = [sigma, sigma];
w = [w, w] * 0.5;
gmm.w  = w;
gmm.mu = mu;
gmm.sigma = sigma;
end
