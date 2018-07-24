function new_llr = score_gmm_trials2(models, testFiles, ubmFilename)
% computes the log-likelihood ratio of observations given the UBM and
% speaker-specific (MAP adapted) models.
%
% Inputs:
%   - models      : a cell array containing the speaker specific GMMs
%   - testFiles   : a cell array containing feature matrices or file names
%   - trials      : a two-dimensional array with model-test verification
%                   indices (e.g., (1,10) means model 1 against test 10)
%   - ubmFilename : file name of the UBM or a structure containing
%					the UBM hyperparameters that is,
%					(ubm.mu: means, ubm.sigma: covariances, ubm.w: weights)
%
% Outputs:
%   - llr		  : log-likelihood ratios for all trials (one score per trial)
%

if ~iscell(models)
    error('Oops! models should be a cell array of structures!');
end

if ischar(ubmFilename)
    tmp = load(ubmFilename);
    ubm = tmp.gmm;
elseif isstruct(ubmFilename)
    ubm = ubmFilename;
else
    error('oh dear! ubmFilename should be either a string or a structure!');
end

if iscellstr(testFiles)
    tlen = length(testFiles);
    tests = cell(tlen, 1);
    for ix = 1 : tlen
        tests{ix} = htkread(testFiles{ix});
    end
elseif iscell(testFiles)
    tests = testFiles;
else
    error('Oops! testFiles should be a cell array!');
end

% number of tests
n_tests = numel(tests);

% ensure enough entries for the data
n_models = numel(models);
% this WAS a generic interger array, but that caused size/memory
% issues. instead we try a cell array to work around the size/memory error
% result = cell(n_tests,1);
% result(:) = {zeros(n_models,1)};
result2 = zeros(n_models,n_tests);

% build out variables from UBM data, never changes so only calculate once!
[ubm_SR, ubm_SMU, ubm_C, ubm_w, f_pi] = compute_ubm_components(ubm);

% this appears to be a slightly better way to calculate the llks. the model
% data is moved the most, which should be the smallest set of data when
% tests becomes large. the UBM is constant so its llk can be known for each
% test. this should work best when the inner for loop takes longer than the
% speed at which the slice of test data + model data can be sent to worker
% nodes. assuming all the nodes get loaded/primed before starting this
% should be a quick way with minimal data movement. testing locally showed
% this is was at least a 50% reduction in time using 4 workers on a laptop.
% an added benefit is that far fewer lines of code are also used now when
% compared to previous versions of the software.
parfor j=1:n_tests
    llk_ubm = compute_llk_ubm( tests{j},ubm_SR,ubm_SMU,ubm_C,ubm_w,f_pi );
    for i=1:n_models
        result2(i,j) = mean(compute_llk(tests{j},models{i},f_pi)-llk_ubm);
    end
end

new_llr = reshape(result2,n_models*n_tests,1);

end

% these can be static, as the UBM does not change once loaded!
function [s_recip,s_mu,C,w,f_pi] = compute_ubm_components(GMM)
s_recip = 1./GMM.sigma;
s_mu = GMM.mu.*s_recip;
C = sum(GMM.mu.*s_mu) + sum(log(GMM.sigma));
w = GMM.w(:);
f_pi = log(2*pi);
end

function llk = compute_llk_ubm(data,s_recip,s_mu,C,w,f_pi)
ndim = size(data,1);
D = s_recip'*(data.*data)-2*s_mu'*data+ndim*f_pi;
logprob = -0.5 * (bsxfun(@plus, C',  D));
logprob = bsxfun(@plus, logprob, log(w));
llk  = logsumexp(logprob, 1);
end

function llk = compute_llk(data, GMM, f_pi)
% compute the posterior probability of mixtures for each frame
post = lgmmprob(data, GMM.mu, GMM.sigma, GMM.w(:), f_pi);
llk  = logsumexp(post, 1);
end

function logprob = lgmmprob(data, mu, sigma, w, f_pi)
% compute the log probability of observations given the GMM
ndim = size(data, 1);
s_recip = 1./sigma;
s_mu = mu.*s_recip;
% C = sum(mu.*mu./sigma) + sum(log(sigma));
C = sum(mu.*s_mu) + sum(log(sigma));
% D = (1./sigma)' * (data .* data) - 2 * (mu./sigma)' ...
% * data  + ndim * log(2 * pi);
D = s_recip'*(data.*data)-2*s_mu'*data+ndim*f_pi;
logprob = -0.5 * (bsxfun(@plus, C',  D));
logprob = bsxfun(@plus, logprob, log(w));
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