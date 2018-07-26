function [SUB_SES, SUB, SES, scores] = ...
    iVectorEvalRocca2(ubms, training_data, testing_data, subjects, sessions)

mixture_count = numel(ubms);

% ubms = allUBMs(training_data(:), mixtures, iterations_ubm, ds_factor,0,0);

% produce I-Vectors for each given UBM mixture
[train_major, train_minor] = size(training_data);
[test_major, test_minor] = size(testing_data);
SUB_SES = zeros(2,mixture_count);
SUB = zeros(subjects,mixture_count);
SES = zeros(sessions,sessions,2,mixture_count);
% dim 1 is cosine, dim 2 is pLDA
scores = zeros(train_major,test_major*test_minor,mixture_count);

% tell it the max workers, just in case?
for m=1:mixture_count
    % extract specific model, NO SAVING!
    ubm = ubms{m};
    stats = cell(train_major,train_minor);
    subject_id = zeros(train_major,train_minor);
    
    parfor sub=1:train_major
        for ses=1:train_minor
            [N,F] = compute_bw_stats2(training_data{sub,ses}, ubm);
            % models are for SUBJECTS only
            stats{sub,ses} = [N;F];
            subject_id(sub,ses) = sub;
        end
    end
    
    % Learn the total variability subspace from all the speaker data.
    tv_dim = 100;
    % given that LDA keeps failing, don't let the tv_dim be LARGE
    tv_dim = min(tv_dim,train_major-1);
    
    % this probably needs to be set higher?
    iterations = 20;
    % full_data = cat(2,stats{:});
    % T cannot be trainined if the the # mixtures * # features > rows of
    % full_data
    
    % stats contains an F array with many ZEROS which is breaking
    % train_tv_space so it cannot produce functional T matrices
    T = train_tv_space2(stats(:), ubm, tv_dim, iterations, 0);
    
    % fprintf('\n ***extract_ivector: enrollment_data *** \n');
    % This builds a cell for each target which contains the I-Vectors for each
    % of its blocks. Since the blocks are internal to the enrollment data, they
    % need to be broken free of the cell structure.
    dev_IVs = zeros(tv_dim, train_major, train_minor);
    
    parfor sub=1:train_major
        for ses=1:train_minor
            dev_IVs(:,sub,ses) = extract_ivector2(stats{sub,ses},ubm,T);
        end
    end
    
    % Now do LDA on the iVectors to find the dimensions that matter.
    lda_dim = min(100, train_major-1);
    devIVbySpeaker = reshape(dev_IVs, tv_dim, train_major*train_minor);
    [V,D] = lda(devIVbySpeaker, subject_id(:));
    if( sum(isnan(D)) > 0 || sum(isinf(D)) > 0 )
        fprintf('\n*** LDA may have failed: eigenvectors are non-real ***\n');
    end
    % individual models
    final_develop_IVs = V(:, 1:lda_dim)' * devIVbySpeaker;
    
    % save model_IVs as those are being compared to the test data
    % fprintf('\n *** compute_bw_stats: testing_data *** \n');
    % Now compute the ivectors for the test set
    test_IVs = zeros(tv_dim, test_major, test_minor);
    parfor sub=1:test_major
        for ses=1:test_minor
            [N, F] = compute_bw_stats2(testing_data{sub,ses}, ubm);
            test_IVs(:,sub,ses) = extract_ivector2([N; F], ubm, T);
        end
    end
    test_IV_Speaker = reshape(test_IVs, tv_dim, test_major*test_minor);
    final_test_IVs = V(:, 1:lda_dim)' * test_IV_Speaker;
    
    % produce scored matrix and evaluate!
    scores(:,:,m) = cosineDistance(final_develop_IVs,final_test_IVs);
    [SUB_SES(:,m), SUB(:,m), SES(:,:,:,m)] = ...
        scoreReport(scores(:,:,m),subjects, sessions,0);
end

end