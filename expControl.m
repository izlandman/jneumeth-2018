function expControl(varargin)

list_folder = varargin{1};
exp_folder = varargin{2};
trial = varargin{3};
exp_flag = varargin{4};
covar_flag = varargin{5};
if( nargin > 5 )
    % store ubm/iv data
    ubmIvStats = zeros(2,2);
    mixtures = varargin{6};
    % ubm iters
    ubmIvStats(1,1) = varargin{7};
    % ds_factor
    ubmIvStats(1,2) = varargin{8};
    % iv iters
    ubmIvStats(2,1) = varargin{9};
    % iv depth
    ubmIvStats(2,2) = varargin{10};
end

% find files of file_lists
file_types = {'CEP','PSD','COH'};
n_types = numel(file_types);
% build listing prior to loop, in case some are missing!
final_list = cell(n_types,1);
final_types = zeros(n_types,1);
for i=1:n_types
    file_list_folder = [list_folder filesep file_types{i} ];
    final_list{i} = findMatchingFiles(file_list_folder,trial);
    final_types(i) = ~isempty(final_list{i});
end
file_list = [final_list{final_types>0}];
n_types = sum(final_types);
% operate on all three feature sets
for i=1:n_types
    fprintf('Using file list: %s\n', file_list{i});
    
    % read file list and import data
    [file_index, subject_count, session_count, full_index] ...
        = fileListConvert(file_list{i});
    
    % find the number of epochs
    epochs = file_index{1,3};
    
    % peek the data to learn more about it
    data_peek = iVectorBinary(full_index{1});
    [feature_count,channel_count] = size(data_peek);
    clear data_peek;
    
    
    %% EXPERIMENTS
    tic
    exp_count = numel(exp_flag);
    for j=1:exp_count
        if( exp_flag(j) == 0 )
            % MAHALANOBIS experiment
            extn = 'MAHAL';
            save_folder = [exp_folder filesep '_RESULT' filesep ...
                file_types{i} filesep extn];
            directoryCheck(save_folder);
            fprintf('Mahalanobis Evaluation.\n');
            [CRR,EER,FPR,FNR] = evalMahal(channel_count,epochs,full_index,...
                subject_count,covar_flag,feature_count);
            % find average (mean) for everything
            mCCR = mean(CRR,2);
            mEER = mean(EER,2);
            mFPR = mean(mean(FPR,3),2);
            mFNR = mean(mean(FNR,3),2);
            a0 = toc;
            fprintf([extn ...
                ' : Data, CCR, and mCRR completed in %f seconds.\n'],a0);
        elseif( exp_flag(j) > 0 )
            if( exp_flag(j) == 1 )
                % GMMUBM experiment
                extn = 'GMMUBM';
            else
                % I-Vector experiment
                extn = 'IVECTOR';
            end
            % produce save folder
            save_folder = [exp_folder filesep '_RESULT' filesep ...
                file_types{i} filesep extn];
            directoryCheck(save_folder);
            fprintf(['La Rocca ', extn ,' Evaluation.\n']);
            % use La Rocca Feature sets
            [CRR, EER, FPR, FNR] = evalGMMiVEC(full_index,subject_count,...
                channel_count,mixtures,epochs,covar_flag,ubmIvStats,...
                exp_flag,save_folder);
            % average over the CV_STEPS and channels for error rates
            mCCR = mean(CRR,3);
            mEER = mean(EER,3);
            mFPR = mean(mean(FPR,4),3);
            mFNR = mean(mean(FNR,4),3);
            a0 = toc;
            fprintf([extn ...
                ' : Data, CCR, and mCRR completed in %f seconds.\n'],a0);
        end
        
        %% save results
        
        % now save the first set of results
        
        % write results to file!
        iVectorBinary([save_folder filesep 'mCCR.bin'],mCCR);
        iVectorBinary([save_folder filesep 'CCR.bin'],CRR);
        iVectorBinary([save_folder filesep 'mEER.bin'],mEER);
        iVectorBinary([save_folder filesep 'EER.bin'],EER);
        iVectorBinary([save_folder filesep 'mFPR.bin'],mFPR);
        iVectorBinary([save_folder filesep 'mFNR.bin'],mFNR);
    end
end

end