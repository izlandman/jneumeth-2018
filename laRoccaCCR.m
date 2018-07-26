% produce CCR from training_data and testing_data, assumes both data
% structures are arrays of (features, channels/elements, epochs)

function [CRR,EER,mahal_dist,FPR,FNR] = ...
    laRoccaCCR(training_data,testing_data,channel)

% setup variables
[features, epochs, subject_count, channel_count] = size(training_data);
mahal_dist = zeros(subject_count, subject_count);

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

% calculate CCR
[m_val,m_index] = sort(mahal_dist);
CRR = (1/subject_count) * sum( [1:subject_count] == m_index(1,:));
answers = eye(subject_count);
[EER, FPR, FNR] = compute_eer_descend(mahal_dist(:),answers(:),0);
end