% just make ONE function to call that parses ALL La Rocca Paper based
% experiments so everything can be streamlined. based upon the arguments
% passed into the function it will sort out what needs to be executed.

function [CRR, mCRR, EER, mEER, mFPR, mFNR, f_CRR, f_mCRR, CRR_bench] ...
    = laRoccaControl(varargin)

if( nargin == 5 )
    fprintf('Passed %d arguments -> MAHAL.\n', nargin);
    file_list = varargin{1};
    fusion_flag = varargin{2};
    save_folder = varargin{3};
    workers = varargin{4};
elseif( nargin == 8 )
    fprintf('Passed %d arguments -> GMM/IVECTOR.\n', ...
        nargin);
    file_list = varargin{1};
    fusion_flag = varargin{2};
    save_folder = varargin{3};
    workers = varargin{4};
    %
    mixtures = varargin{5};
    iterations = varargin{6};
    ds_factor = varargin{7};
    eval_flag = varargin{8};
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
if isunix
    cluster.JobStorageLocation = getenv('HOME');
end
parpool(cluster,workers, 'IdleTimeout', Inf);

% read file list and import data
[file_index, subject_count, session_count, full_index] ...
    = fileListConvert(file_list);

% build variables
channel_count = file_index{1,3};

% based upon channels. 22 for cepstrum, 56 & 1540 for PSD & COH
data_peek = iVectorBinary(full_index{1});
epochs = channel_count;
[feature_count,channel_count] = size(data_peek);
clear data_peek;

% how do we produce the data?
if( nargin == 5 )
    fprintf('Mahal Evaluation.\n');
    [CRR, EER, mahal_dist, FPR, FNR] = mahalResults(channel_count, ...
        epochs,full_index,subject_count,feature_count);
    % find average (mean) for everything
    mCRR = mean(CRR,2);
    mEER = mean(EER,2);
    mFPR = mean(mean(FPR,3),2);
    mFNR = mean(mean(FNR,3),2);
    a0 = toc;
    fprintf('Data, CRR, and mCRR completed in %f seconds.\n',a0);
elseif( nargin == 8 )
    fprintf('GMM/I-VECTOR Evaluation.\n');
    % use La Rocca Feature sets
    [CRR, EER, scores, FPR, FNR] = gmmIvecResults(channel_count, ...
        epochs, full_index, subject_count, mixtures, eval_flag, ...
        iterations, ds_factor);
    % average over the CV_STEPS and channels for error rates
    mCRR = mean(CRR,3);
    mEER = mean(EER,3);
    mFPR = mean(mean(FPR,4),3);
    mFNR = mean(mean(FNR,4),3);
    a0 = toc;
    fprintf('Data, CRR, and mCRR completed in %f seconds.\n',a0);
else
    fprintf('Incorrect number of input arguments!\n');
end

% now save the first set of results
% write results to file!
directoryCheck(save_folder);
iVectorBinary([save_folder filesep 'mCRR.bin'],mCRR);
iVectorBinary([save_folder filesep 'CRR.bin'],CRR);
iVectorBinary([save_folder filesep 'mEER.bin'],mEER);
iVectorBinary([save_folder filesep 'EER.bin'],EER);
iVectorBinary([save_folder filesep 'mFPR.bin'],mFPR);
iVectorBinary([save_folder filesep 'mFNR.bin'],mFNR);

%% fusion match scoring
if( nargin == 4 )
    [f_CRR, f_mCRR, CRR_bench, a1] = fusionMatchScoring(fusion_flag, mCRR, ...
        epochs, subject_count, mahal_dist, save_folder);
    fprintf('Classification completed in %f seconds.\n', a0+a1);
else
    [f_CRR, f_mCRR, CRR_bench, a1] = fusionMatchScoring(fusion_flag, mCRR, ...
        epochs, subject_count, scores, save_folder, mixtures);
    fprintf('Classification completed in %f seconds.\n', a0+a1);
end

end