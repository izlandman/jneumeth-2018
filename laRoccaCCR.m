% produce CCR from training_data and testing_data, assumes both data
% structures are arrays of (features, channels/elements, epochs)

function [CCR,EER,mahal_dist,FPR,FNR] = ...
    laRoccaCCR(training_data,testing_data,channel,covar_flag)

% setup variables
[features, epochs, subject_count, channel_count] = size(training_data);
mahal_dist = zeros(subject_count, subject_count);

% pooled covariance is either built from the OTHER channels of the subject
% OR all channels from all subjects (minus the subject's channel)


if( covar_flag == 0 )
    % pooled from 'class' which we take to mean subject
    for s=1:subject_count
        % mu is built from the channel-subject specific training data
        class_mu = mean( training_data(:,:,s,channel),2 );
        % sigma is built from all channels of the subject
        class_sigma = cov( (reshape( squeeze(training_data(:,:,s,:)), ...
            features, epochs*channel_count)-class_mu)' );
        % produce statistical object
        mahal_obj = gmdistribution(class_mu', class_sigma);
        % evaluate each test vector to object
        for t=1:subject_count
            mahal_dist(t,s) = mahal(mahal_obj, testing_data(:,t,channel)');
        end
    end
elseif( covar_flag == 1 )
    % pooled from ALL subjects channel data
    for s=1:subject_count
        % data from the subject, but sourced from alternative channels
        class_mu = mean( reshape(training_data(:,:,s,:), features, ...
            epochs*channel_count),2 );
        % zero center with respect to class mu
        class_zero_centered = permute(training_data(:,:,:,channel), ...
            [1 3 2]) - class_mu;
        % data from all subjects of the same channel
        class_sigma = cov( (reshape( class_zero_centered, ...
            features, epochs*subject_count) - class_mu)' );
        mahal_obj = gmdistribution(class_mu', class_sigma);
        for t=1:subject_count
            mahal_dist(s,t) = mahal(mahal_obj, testing_data(:,t,channel)');
        end
    end
elseif( covar_flag == - 1 )
    % this works the best, why?
    
    % class_mu is an average over all subjects for the given channel
    class_mu = squeeze(mean(training_data(:,:,:,channel),2));
    % double zero centered as cov zero centers as well, but it works since
    % this is zero centering the data against each the mean for each epoch
    class_zero_centered = permute(training_data(:,:,:,channel),[1 3 2]) ...
        - class_mu;
    % find the pooled covariance AND mean center again
    class_sigma = cov(reshape(class_zero_centered, features, ...
        subject_count*epochs)');
    for s=1:subject_count
        mahal_obj = gmdistribution(class_mu(:,s)', class_sigma);
        for t=1:subject_count
            mahal_dist(s,t) = mahal(mahal_obj, testing_data(:,t,channel)');
        end
    end
end

% calculate CCR
[m_val,m_index] = sort(mahal_dist);
CCR = (1/subject_count) * sum( [1:subject_count] == m_index(1,:));
answers = eye(subject_count);
[EER, FPR, FNR] = compute_eer_descend(mahal_dist(:),answers(:),0);
end