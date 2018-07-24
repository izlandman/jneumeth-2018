% handle all mahalanobis result production
function [CCR, EER, mahal_dist, FPR, FNR] = mahalResults(channel_count, ...
    epochs, full_index, subject_count, covar_flag, feature_count)

CCR = zeros(channel_count, epochs);
EER = CCR;
x_subjects = subject_count * subject_count;
FPR = zeros(x_subjects, channel_count, epochs);
FNR = FPR;
mahal_dist = zeros(subject_count, subject_count, channel_count, epochs);

if( channel_count == 22 )
    
    az = tic;
    for cv_step = 1:epochs
        % fprintf('>>> Iteration %d of %d.\n', cv_step, epochs);
        % create training and testing indexes
        [tra_data, tes_data] = ...
            laRoccaDataBuilderCepstrum(full_index, channel_count, ...
            subject_count, epochs, cv_step);
        parfor c=1:channel_count
            % calculate CCR of testing data from training data
            [CCR(c,cv_step), EER(c,cv_step), mahal_dist(:,:,c,cv_step), ...
                FPR(:,c,cv_step), FNR(:,c,cv_step)] = ...
                laRoccaCCR(tra_data, tes_data, c, covar_flag);
        end
        fprintf('>>> Iteration %d of %d completed in %f seconds.\n', ...
            cv_step, epochs, toc(az));
    end
    
elseif( channel_count == 56 || channel_count == 1540 )
    
    az = tic;
    for cv_step = 1:epochs
        % fprintf('>>> Iteration %d of %d.\n', cv_step, epochs);
        % create training and testing indexes
        [tra_data, tes_data] = ...
            laRoccaDataBuilder(full_index,feature_count,channel_count, ...
            subject_count,epochs,cv_step);
        % calculate CCR of testing data from training data
        parfor c=1:channel_count
            [CCR(c,cv_step), EER(c,cv_step), mahal_dist(:,:,c,cv_step), ...
                FPR(:,c, cv_step), FNR(:, c, cv_step)] = ...
                laRoccaCCR(tra_data, tes_data, c, covar_flag);
        end
        fprintf('>>> Iteration %d of %d completed in %f seconds.\n', ...
            cv_step, epochs, toc(az));
    end
    
end

end