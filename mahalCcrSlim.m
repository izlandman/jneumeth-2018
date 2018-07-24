% produce CCR from training_data and testing_data, assumes both data
% structures are arrays of (features, channels/elements, epochs). a slim
% version that does not focus on channels, but instead builds models on the
% subject level from the training data and evaluates on the 'subjects' in
% the testing data. this data most likely comes in as cell arrays

% training_data and testing_data are cells
%   size SUBJECT x SESSION
%   cell elements -> FEATURES x EPOCHS
% covar_flag is an integer

function [SUB_SES, SUB, SES, mahal_dist] = ...
    mahalCcrSlim(training_data, testing_data, subjects, sessions)

% setup variables
[n_subject_tra, n_session_tra] = size(training_data);
[n_feat_tra,~] = size(training_data{1,1});
[n_subject_tes, n_session_tes] = size(testing_data);
mahal_dist = zeros(n_subject_tra, n_subject_tes, n_session_tes);

% build mahal model
class_mu  = zeros(n_subject_tra, n_feat_tra);
class_sigma = zeros(n_feat_tra, n_feat_tra, n_subject_tra);
for sub=1:n_subject_tra
    class_mu(sub,:) = mean(cell2mat(training_data(sub,:)),2);
    class_sigma(:,:,sub) = cov(cell2mat(training_data(sub,:))');
end

parfor s=1:n_subject_tra
    % produce statistical object
    mahal_obj = gmdistribution(class_mu(s,:), class_sigma(:,:,s));
    % evaluate each test vector to object
    for t_ses=1:n_session_tes
        for t_sub=1:n_subject_tes
            mahal_dist(s,t_sub,t_ses) = mahal(mahal_obj, ...
                mean(cell2mat(testing_data(t_sub,t_ses)),2)');
        end
    end
end

% CCR for averaged subject, individual subject, and individual session
[SUB_SES, SUB, SES] = scoreReport(mahal_dist,subjects, sessions, 1);
end