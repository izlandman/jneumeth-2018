% gather the data from the  La Rocca feature files
% output data is an array!
%   size is FEATURE x EPOCHS x SUBJECTS x CHANNELS
function [training_data, testing_data] = laRoccaDataBuilderCepstrum ...
    (full_index, channel_count, subject_count, epochs, cv_step)

data_peek = htkread(full_index{1});
[features,~] = size(data_peek);

feature_ind = 1:features;
block_size = 100;
window_overlap = 0;

training_data = zeros(features, epochs-1, subject_count, channel_count);
testing_data = zeros(features, subject_count, channel_count);

% create training and testing indexes
training = 1:epochs;
testing = cv_step;
training = training( ~ismember(training,testing) );

% load the training data and testing data
parfor s=1:subject_count
    % needed for debugging error in epoch production
    % fprintf('Build epochs for subject: %d\n', s);
    for c=1:channel_count
    f_index = c + (s-1)*channel_count;
    [training_data(:,:,s,c),testing_data(:,s,c)] = eegTestTrainSplit( ...
        full_index{f_index}, feature_ind, block_size, window_overlap, ...
        training, testing);
    end
end

end