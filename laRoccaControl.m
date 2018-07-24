% just make ONE function to call that parses ALL La Rocca Paper based
% experiments so everything can be streamlined. based upon the arguments
% passed into the function it will sort out what needs to be executed.

function [CCR, mCCR, EER, mEER, mFPR, mFNR, f_CCR, f_mCCR, CCR_bench] ...
    = laRoccaControl(varargin)

if( nargin == 5 )
    fprintf('Passed %d arguments. La Rocca Features & MAHAL.\n', nargin);
    file_list = varargin{1};
    fusion_flag = varargin{2};
    save_folder = varargin{3};
    covar_flag = varargin{4};
    workers = varargin{5};
    block_size = 100;
    window_overlap = 0;
elseif( nargin == 9 )
    fprintf('Passed %d arguments. La Rocca Features & GMM/IVECTOR.\n', ...
        nargin);
    file_list = varargin{1};
    fusion_flag = varargin{2};
    save_folder = varargin{3};
    covar_flag  = varargin{4};
    workers = varargin{5};
    %
    mixtures = varargin{6};
    iterations = varargin{7};
    ds_factor = varargin{8};
    eval_flag = varargin{9};
elseif( nargin == 12 )
    fprintf('Passed %d arguments. Cepstrum Features.\n', nargin);
    file_list = varargin{1};
    fusion_flag = varargin{2};
    save_folder = varargin{3};
    covar_flag = varargin{4};
    workers = varargin{5};
    %
    mixtures = varargin{6};
    iterations = varargin{7};
    ds_factor = varargin{8};
    eval_flag = varargin{9};
    %
    block_size = varargin{10};
    window_overlap = varargin{11};
    block_style = varargin{12};
else
    fprintf('Incorrect number of arguments [%d] passed!\n', nargin);
    return;
end

%% Mimic four argument functions

fprintf('Starting classificaiton!\n');
tic
% enhanced debug feature
% setSchedulerMessageHandler(@disp);
% clear existing cluster
delete(gcp('nocreate'))
% setup parallel environment
cluster = parcluster('local');
cluster.NumWorkers = workers;
% this needs to be ON for cluster operation and OFF for local debugging
cluster.JobStorageLocation = getenv('HOME');
parpool(cluster,workers, 'IdleTimeout', Inf);

% read file list and import data
[file_index, subject_count, session_count, full_index] ...
    = fileListConvert(file_list);

% build variables
channel_count = file_index{1,3};

% based upon channels. 22 for cepstrum, 56 & 1540 for PSD & COH
if channel_count == 22
    data_peek = htkread(full_index{1});
    [feature_count,samples] = size(data_peek);
    epochs = floor((samples-window_overlap) / (block_size-window_overlap));
else
    data_peek = iVectorBinary(full_index{1});
    epochs = channel_count;
    [feature_count,channel_count] = size(data_peek);
end
clear data_peek;

% how do we produce the data?
if( nargin == 5 )
    fprintf('Mahal Evaluation.\n');
    [CCR, EER, mahal_dist, FPR, FNR] = mahalResults(channel_count, ...
        epochs, full_index, subject_count, covar_flag, feature_count);
    % find average (mean) for everything
    mCCR = mean(CCR,2);
    mEER = mean(EER,2);
    mFPR = mean(mean(FPR,3),2);
    mFNR = mean(mean(FNR,3),2);
    a0 = toc;
    fprintf('Data, CCR, and mCRR completed in %f seconds.\n',a0);
elseif( nargin == 9 )
    fprintf('La Rocca GMM/I-VECTOR Evaluation.\n');
    % use La Rocca Feature sets
    [CCR, EER, scores, FPR, FNR] = gmmIvecResults(channel_count, ...
        epochs, full_index, subject_count, mixtures, eval_flag, ...
        iterations, ds_factor, covar_flag);
    % average over the CV_STEPS and channels for error rates
    mCCR = mean(CCR,3);
    mEER = mean(EER,3);
    mFPR = mean(mean(FPR,4),3);
    mFNR = mean(mean(FNR,4),3);
    a0 = toc;
    fprintf('Data, CCR, and mCRR completed in %f seconds.\n',a0);
elseif( nargin == 12 )
    fprintf('Cepstrum GMM/I-VECTOR Evaluation.\n');
    % use Cepstrum Feature Sets
    [CCR, EER, scores, FPR, FNR] = gmmIvecResults(channel_count, ...
        epochs, full_index, subject_count, mixtures, eval_flag, ...
        iterations, ds_factor, covar_flag, block_style, session_count, ...
        block_size, window_overlap);
    % average over the CV_STEPS and channels for error rates
    mCCR = mean(CCR,3);
    mEER = mean(EER,3);
    mFPR = mean(mean(FPR,4),3);
    mFNR = mean(mean(FNR,4),3);
    a0 = toc;
    fprintf('Data, CCR, and mCRR completed in %f seconds.\n',a0);
end

% now save the first set of results
% write results to file!
iVectorBinary([save_folder filesep 'mCCR.bin'],mCCR);
iVectorBinary([save_folder filesep 'CCR.bin'],CCR);
iVectorBinary([save_folder filesep 'mEER.bin'],mEER);
iVectorBinary([save_folder filesep 'EER.bin'],EER);
iVectorBinary([save_folder filesep 'mFPR.bin'],mFPR);
iVectorBinary([save_folder filesep 'mFNR.bin'],mFNR);

%% fusion match scoring
if( nargin == 5 )
    [f_CCR, f_mCCR, CCR_bench, a1] = fusionMatchScoring(fusion_flag, mCCR, ...
        epochs, subject_count, mahal_dist, save_folder);
    fprintf('Classification completed in %f seconds.\n', a0+a1);
else
    [f_CCR, f_mCCR, CCR_bench, a1] = fusionMatchScoring(fusion_flag, mCCR, ...
        epochs, subject_count, scores, save_folder, mixtures);
    fprintf('Classification completed in %f seconds.\n', a0+a1);
end

end