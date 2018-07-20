% carry out GMM-UBM and I-Vector evaluation of rocca style features
function [CRR,EER,FPR,FNR] = evalGMMiVEC(full_index,subject_count,...
    channel_count,mixtures,epochs,covar_flag,ubmIvStats,eval_flag,...
    save_folder)

n_mixtures = numel(mixtures);
x_subjects = subject_count * subject_count;
CRR = zeros(channel_count, n_mixtures, epochs);
EER = zeros(channel_count, n_mixtures, epochs);
FPR = zeros(x_subjects, channel_count, n_mixtures, epochs);
FNR = FPR;

% break apart the iterations
ubm_iters = ubmIvStats(1,1);
ds_factor = ubmIvStats(1,2);
iv_iters =  ubmIvStats(2,1);
iv_depth = ubmIvStats(2,2);

% flag training for all mixtures
ubm_index = 0;
iv_index = 0;

parfor cv_step = 1:epochs
    az = tic;
    [training_data, testing_data] = laRoccaDataBuilderUBM( ...
        full_index, subject_count, epochs, cv_step);
    % build object!
    UBM = UniversalBackgroundModel(mixtures,ubm_iters,ds_factor,...
        subject_count,channel_count);
    % build UBMs
    train(UBM,training_data(:),0);
    % build TVMs if evaluating I-Vectors
    if( eval_flag == 2 )
        % produce IVECTOR based results through TVM class
        TVM = TotalVariabilityMatrix(UBM.UBMs,iv_depth,iv_iters,...
            channel_count);
        % train the models
        train(TVM,training_data,iv_index);
    end
    %%% remember ubm/gmm/ivector wants cell array of subject x channels
    %%% with each entry being the features x samples
    for c=1:channel_count
        % the ENROLLMENT data can come from the full set of channels OR the
        % specific channel used in TESTING
        if( covar_flag == 0 )
            % this is the defaul because it should be harder given less
            % samples
            enrollment_data = training_data(:,c);
        else
            enrollment_data = training_data;
        end
        if( eval_flag == 1 )
            % handle evaluation internal to UBM class!
            evaluate(UBM,enrollment_data,testing_data(:,c),c,ubm_index);
        elseif( eval_flag == 2 )
            % handle evaluation internal to TVM class!
            evaluate(TVM,enrollment_data,testing_data(:,c),c,iv_index);
        else
            fprintf('Bad eval_flag!\n');
        end
    end
    % save UBM regardless of GMMUBM or I-Vector evaluation
    save_name = ['UBM_e' num2str(cv_step) '.mat'];
    saveobj(UBM,save_folder,save_name);
    if( eval_flag == 1 )
        CRR(:,:,cv_step) = UBM.crr;
        EER(:,:,cv_step) = UBM.eer;
        FPR(:,:,:,cv_step) = UBM.fpr;
        FNR(:,:,:,cv_step) = UBM.fnr;
    else
        save_name = ['TVM_e' num2str(cv_step) '.mat'];
        saveobj(TVM,save_folder,save_name);
        CRR(:,:,cv_step) = TVM.crr;
        EER(:,:,cv_step) = TVM.eer;
        FPR(:,:,:,cv_step) = TVM.fpr;
        FNR(:,:,:,cv_step) = TVM.fnr;
    end
    fprintf('>>> Iteration %d of %d completed in %f seconds.\n', ...
        cv_step, epochs, toc(az));
end

end