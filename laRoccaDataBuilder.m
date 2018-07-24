% gather the data from the  La Rocca feature files
function [training_data, testing_data] = ...
    laRoccaDataBuilder(full_index, feature_count, channel_count, ...
    subject_count, epochs, cv_step)

% output data
training_data = zeros(feature_count,epochs-1,subject_count,channel_count);
testing_data = zeros(feature_count,subject_count,channel_count);

% create training and testing indexes
training = 1:epochs;
testing = cv_step;
training = training( ~ismember(training,testing) );

% load the training data and testing data
for s=1:subject_count
    test_index = (s-1)*epochs + testing;
    train_index = (s-1)*epochs + training;
    testing_data(:,s,:) = iVectorBinary(full_index{test_index});
    for e=1:epochs-1
        training_data(:,e,s,:) = ...
            iVectorBinary(full_index{train_index(e)});
    end
end

end