% produce all necessary CCR, EER, and ROC data
function [SUB_SES_res, SUB_res, SES_res] = scoreReport(scores,subjects,...
    sessions,sort_flag)

[n_subject_tra,n_subject_tes] = size(scores);
if( sort_flag == 0 )
    % GMM/IVECTOR
    sort_flag = 'descend';
else
    % MAHAL
    sort_flag = 'ascend';
end

% find subject-session ccr
[~, subject_data_sort] = sort(scores,sort_flag);
sub_ses_ccr = sum( subject_data_sort(1,:) == 1:n_subject_tra)/...
    n_subject_tra;
% accuracy for a subject independent of session
subject_ccr = zeros(subjects,1);
for sub=1:subjects
    sub_ind = sub:subjects:n_subject_tra;
    subject_ccr(sub) = sum(ismember(sub_ind,...
        subject_data_sort(1,sub_ind)))/sessions;
end
% subject-session eer
answers_1 = eye(n_subject_tra);
if( strcmp(sort_flag,'descend') )
    [sub_ses_eer, sub_ses_FPR, sub_ses_FNR] = ...
        compute_eer_2(scores(:),answers_1(:),0);
else
    [sub_ses_eer, sub_ses_FPR, sub_ses_FNR] = ...
        compute_eer_descend(scores(:),answers_1(:),0);
end
% reshape into session segments for session based ccr and eer
scores = reshape(scores,subjects,sessions,subjects,sessions);
session_ccr = zeros(sessions,sessions);
session_eer = zeros(sessions,sessions);
answers_2 = eye(subjects);
for ses_o=1:sessions
    for ses_i=1:sessions
        data = squeeze(scores(:,ses_o,:,ses_i));
        [~,data_ind] = sort( data, sort_flag );
        session_ccr(ses_o,ses_i) = sum(data_ind(1,:)==1:subjects)/...
            subjects;
        if( strcmp(sort_flag,'descend') )
            session_eer(ses_o,ses_i) = compute_eer_2(data(:),...
                answers_2(:),0);
        else
            session_eer(ses_o,ses_i) = compute_eer_descend(data(:),...
                answers_2(:),0);
        end
    end
end

SUB_SES_res = [sub_ses_ccr;sub_ses_eer];
SUB_res = subject_ccr;
SES_res = cat(3,session_ccr,session_eer);

end