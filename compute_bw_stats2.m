%%% EDITED!
% This returns values of 0 for probability which causes all sorts of
% problems in the code later during matrix math. Instead, nothing will be
% assigned a probability of 0 but instead be given a probability equal to
% half of the smallest seen by the other mixtures. This should improve the
% the robustness of the math while sacrificing a minimal hit on accuracy.

function [N, F] = compute_bw_stats2(feaFilename, ubmFilename, statFilename)
% extracts sufficient statistics for features in feaFilename and GMM 
% ubmFilename, and optionally save the stats in statsFilename. The 
% first order statistics are centered.
%
% Inputs:
%   - feaFilename  : input feature file name (string) or a feature matrix 
%					(one observation per column)
%   - ubmFilename  : file name of the UBM or a structure with UBM 
%					 hyperparameters.
%   - statFilename : output file name (optional)   
%
% Outputs:
%   - N			   : mixture occupation counts (responsibilities) 
%   - F            : centered first order stats
%
% References:
%   [1] N. Dehak, P. Kenny, R. Dehak, P. Dumouchel, and P. Ouellet, "Front-end 
%       factor analysis for speaker verification," IEEE TASLP, vol. 19, pp. 788-798,
%       May 2011. 
%   [2] P. Kenny, "A small footprint i-vector extractor," in Proc. Odyssey, 
%       The Speaker and Language Recognition Workshop, Jun. 2012.
%
%
% Omid Sadjadi <s.omid.sadjadi@gmail.com>
% Microsoft Research, Conversational Systems Research Center

if ischar(ubmFilename),
	tmp  = load(ubmFilename);
	ubm  = tmp.gmm;
elseif isstruct(ubmFilename),
	ubm = ubmFilename;
else
    error('Oops! ubmFilename should be either a string or a structure!');
end
[ndim, nmix] = size(ubm.mu);
m = reshape(ubm.mu, ndim * nmix, 1);
idx_sv = reshape(repmat(1 : nmix, ndim, 1), ndim * nmix, 1);

if ischar(feaFilename),
    data = htkread(feaFilename);
else
    data = feaFilename;
end

[N, F] = expectation(data, ubm);
F = reshape(F, ndim * nmix, 1);
F = F - N(idx_sv) .* m; % centered first order stats

if ( nargin == 3)
	% create the path if it does not exist and save the file
	path = fileparts(statFilename);
	if ( exist(path, 'dir')~=7 && ~isempty(path) ), mkdir(path); end
	parsave(statFilename, N, F);
end
end

function parsave(fname, N, F) %#ok
save(fname, 'N', 'F')
end

function [N, F] = expectation(data, gmm)
% compute the sufficient statistics
post = postprob(data, gmm.mu, gmm.sigma, gmm.w(:));
N = sum(post, 2);
F = data * post';
end

function [post, llk] = postprob(data, mu, sigma, w)
% compute the posterior probability of mixtures for each frame
post = lgmmprob(data, mu, sigma, w);
llk  = logsumexp(post, 1);
post = exp(bsxfun(@minus, post, llk));
% this post value can produce very small values when a given mixture is
% seen as a 'correct match. eventually this is problematic because it
% introduces zero or very small near zero vaues which are troublesome.

% tweek is %!
tweek = 0.00001;
% realize that some values are 1.000 and not 1! Those given as 1.000 appear
% to be naturally calculated and not as poorly conditioned in later steps
post_flag = post==1;
col_flag = sum(post_flag)>0;
post_flag2 = ~post_flag .* repmat(col_flag,numel(w),1);
% take % from perfect match
post(post_flag==1) = post(post_flag==1)-tweek;
% spread to remaining mixtures
post(post_flag2==1) = post(post_flag2==1)+(tweek/(numel(w)-1));
end

function logprob = lgmmprob(data, mu, sigma, w)
% compute the log probability of observations given the GMM
ndim = size(data, 1);
C = sum(mu.*mu./sigma) + sum(log(sigma));
D = (1./sigma)' * (data .* data) - 2 * (mu./sigma)' * data  + ndim * log(2 * pi);
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