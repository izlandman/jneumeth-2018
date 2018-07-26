function [CRR,EER,scores,FPR,FNR] = ...
    ubmGmmEvalRocca(training_data,testing_data,ubms,channel)

mixture_count = numel(ubms);
% build UBM from all data. other variations would be to build it from all
% subjects but the given channel or build one for each subject using all
% their channels?

% build training_data from only the channel or from all the subject's
% channels.
training_data = training_data(:,channel);


[speakers, channels] = size(testing_data);
[features,epochs] = size(training_data{1});
EER = zeros(1,mixture_count);
CRR = zeros(1,mixture_count);
scores = zeros(speakers, speakers, mixture_count);
FNR = zeros(speakers*speakers,mixture_count);
FPR = FNR;

for m=1:mixture_count
    %a0 = tic;
    % ubms contains a mixture of size 1! ignore it!
    ubm = ubms{m};
    % build gmms for each speaker in relation to the ubm
    map_tau = 19;
    config = 'mwv';
    
    % should exist since it is built above by gmm_em?
    gmm_speakers = cell(speakers,1);
    for i=1:speakers
        % mapAdapt wants a cell of subjects x channels!
        gmm_speakers{i} = mapAdapt2(training_data(i,:),ubm,map_tau,config);
    end
    
    % these become too large when dealing with 30 second epochs!
    % trials = zeros(speakers*channels*speakers, 2);
    answers = zeros(speakers*channels*speakers, 1);
    
    for ix = 1 : speakers
        % need to make 'results' on a subject by subject basis to reduce
        % memory overhead
        % these shouldn't break the memory limits
        % trials = zeros(channels*speakers, 2);
        % answers = zeros(channels*speakers, 1);
        % old stuff, that now must be changed, or at least remember this is
        % happened for only one speaker at a time!
        b = (ix-1)*speakers*channels + 1;
        % trials(b:e, :)  = [ix * ones(speakers*channels, 1),...
        %    (1:speakers*channels)'];
        answers((ix-1)*channels+b : (ix-1)*channels+b+channels-1) = 1;
    end
    
    % this returns a cell array, n_subjects and then n_trials
    %a1= tic;
    gmm_scores = score_gmm_trials2(gmm_speakers, ...
        reshape(testing_data', speakers*channels,1), ubm);
    %fprintf('score_gmm_trials2: %f s\n', toc(a1));
    
    scores(:,:,m) = reshape(gmm_scores,speakers,speakers);
    [~,ind] = sort(scores(:,:,m),'descend');
    CRR(m) = sum( ind(1,:) == [1:speakers] ) / speakers;
    [EER(m), FPR(:,m), FNR(:,m)] = compute_eer_2(gmm_scores,answers,0);
    %fprintf('mixture %d: %f s\n', m, toc(a0));
end

end