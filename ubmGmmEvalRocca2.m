function [SUB_SES, SUB, SES, scores] = ...
    ubmGmmEvalRocca2(ubms, training_data, testing_data, subjects, sessions)

mixture_count = numel(ubms);
% build UBM from all data. other variations would be to build it from all
% subjects but the given channel or build one for each subject using all
% their channels?
% ubms = allUBMs(training_data(:), mixtures, iterations_ubm, ds_factor,0,0);

[subject_tra, session_tra] = size(training_data);
[subject_tes, session_tes] = size(testing_data);
SUB_SES = zeros(2,mixture_count);
SUB = zeros(subjects,mixture_count);
SES = zeros(sessions,sessions,2,mixture_count);
scores = zeros(subject_tra,subject_tes,mixture_count);
for m=1:mixture_count
    % ubms contains a mixture of size 1! ignore it!
    ubm = ubms{m};
    % build gmms for each speaker in relation to the ubm
    map_tau = 10;
    config = 'mwv';
    gmm_subjects = cell(subject_tra,1);
    
    % should exist since it is built above by gmm_em?
    for ses=1:subject_tra
        % mapAdapt wants a cell of subjects x channels!
        gmm_subjects{ses,1} = mapAdapt2(training_data(ses,:), ...
            ubm, map_tau,config);
    end
      
    % this returns a cell array, n_subjects and then n_trials
    gmm_scores = zeros(subject_tra,subject_tes);
    for sub=1:subject_tes
        gmm_scores(:,sub) = score_gmm_trials2(gmm_subjects, ...
            testing_data(sub,:), ubm);
    end
    scores(:,:,m) = gmm_scores;
    % average models, find subject CCR
    [SUB_SES(:,m), SUB(:,m), SES(:,:,:,m)] = ...
        scoreReport(gmm_scores,subjects,sessions,0);    
end

end