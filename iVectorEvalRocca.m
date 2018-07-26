function [CRR,EER,scores,FPR,FNR] = ...
    iVectorEvalRocca(training_data,testing_data,ubms,channel)

mixture_count = numel(ubms);

% build training_data from only the channel or from all the subject's
% channels.
training_data = training_data(:,channel);

% produce I-Vectors for each given UBM mixture
[speakers, channels] = size(testing_data);
[features,~] = size(testing_data{1});
EER = zeros(1,mixture_count);
CRR = zeros(1,mixture_count);
scores = zeros(speakers,speakers,mixture_count);
FPR = zeros(speakers*speakers,mixture_count);
FNR = FPR;

% tell it the max workers, just in case?
for m=1:mixture_count
    % extract specific model, NO SAVING!
    ubm = ubms{m};
    stats = cell(speakers,1);
    
    % this throws an error about 'too many arguments'
    parfor t=1:speakers
        stats{t} = computeBwStatsFromBlocks(training_data{t,:}, ubm);
    end
    
    % Learn the total variability subspace from all the speaker data.
    % fprintf('\n *** train_tv_space: enrollment_data *** \n');
    tic
    tv_dim = 100;
    % given that LDA keeps failing, don't let the tv_dim be LARGE
    tv_dim = min(tv_dim,speakers-1);
    
    % this probably needs to be set higher?
    iterations = 20;
    full_data = cat(2,stats{:});
    % T cannot be trainined if the the # mixtures * # features > rows of
    % full_data
    T = train_tv_space4(full_data, ubm, tv_dim, iterations, 0);
    
    % fprintf('\n ***extract_ivector: enrollment_data *** \n');
    % This builds a cell for each target which contains the I-Vectors for each
    % of its blocks. Since the blocks are internal to the enrollment data, they
    % need to be broken free of the cell structure.
    enrollment_stats = cat(2,stats{:});
    [~,enrollment_blocks] = size(enrollment_stats);
    develop_IVs = cell(enrollment_blocks,1);
    
    parfor t=1:enrollment_blocks
        develop_IVs{t} = extract_ivector3(enrollment_stats(:,t), ubm, T);
    end
    
    % save develop_IVs
    % Now do LDA on the iVectors to find the dimensions that matter.
    lda_dim = min(tv_dim, enrollment_blocks-1);
    % find blocks in each target
    block_counts = cellfun(@size,stats(:),'UniformOutput',0);
    block_counts = sum(cell2mat(block_counts));
    speaker_id = 1:block_counts(2);
    develop_IV_Speaker = cat(2,develop_IVs{:});
    [V,D] = lda(develop_IV_Speaker, speaker_id(:)');
    % If D contains non-real values, lda probably failed to adequately resolve
    % itself from the given data
    if( sum(isnan(D)) > 0 || sum(isinf(D)) > 0 )
        fprintf('\n*** LDA may have failed: eigenvectors are non-real ***\n');
    end
    % apparently V can come back with imaginary values?
    V = real(V);
    final_develop_IVs = V(:, 1:lda_dim)' * develop_IV_Speaker;    
    
    % save model_IVs as those are being compared to the test data
    % fprintf('\n *** compute_bw_stats: testing_data *** \n');
    % Now compute the ivectors for the test set
    [test_blocks,~] = size(testing_data);
    test_IVs = cell(test_blocks, 1);
    parfor t=1:test_blocks
        temp_stats = computeBwStatsFromBlocks(testing_data{t}, ubm);
        test_IVs{t} = extract_ivector3(temp_stats, ubm, T);
    end
    test_IV_Speaker = cat(2,test_IVs{:});
    final_test_IVs = V(:, 1:lda_dim)' * test_IV_Speaker;
    
    % produce scored matrix and evaluate!
    scores = cosineDistance(final_develop_IVs, final_test_IVs);
    [~,ind] = sort(scores,'descend');
    CRR(m) = sum( ind(1,:) == [1:speakers] ) / speakers;
    answers = eye(speakers);
    [EER(m), FPR(:,m), FNR(:,m)] = compute_eer_2(scores(:),answers(:),0);  
    scores(:,:,m) = scores;
end

end