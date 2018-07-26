% handle all mahalanobis result production
function [CRR, EER, scores, FPR, FNR] = mahalResults(channel_count, ...
    epochs, full_index, subject_count, feature_count)

CRR = zeros(channel_count, epochs);
EER = CRR;
x_subjects = subject_count * subject_count;
FPR = zeros(x_subjects, channel_count, epochs);
FNR = FPR;
scores = zeros(subject_count, subject_count, channel_count, epochs);

az = tic;
parfor cv_step = 1:epochs
    % fprintf('>>> Iteration %d of %d.\n', cv_step, epochs);
    % create training and testing indexes
    [tra_data, tes_data] = ...
        laRoccaDataBuilder(full_index,feature_count,channel_count, ...
        subject_count,epochs,cv_step);
    % calculate CCR of testing data from training data
    for c=1:channel_count
        [CRR(c,cv_step),EER(c,cv_step),scores(:,:,c,cv_step), ...
            FPR(:,c, cv_step),FNR(:, c, cv_step)] = ...
            laRoccaCCR(tra_data,tes_data,c);
    end
    fprintf('>>> Iteration %d of %d completed in %f seconds.\n', ...
        cv_step, epochs, toc(az));
end

end