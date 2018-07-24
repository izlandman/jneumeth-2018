% prep data allow for altering how UBM is made!
function [testing_data, enrollment_data] = ...
    dataUbmPrep(subject_count, session_count, full_index, features, ...
    block_size, window_overlap, epochs, cv_step)

% size everything out!
channel_count = numel(full_index)/session_count;
enrollment_data = cell(subject_count,channel_count);
testing_data = cell(subject_count,channel_count);

% organize epoch structure
enrollment_blocks = 1:epochs;
testing_block = cv_step;
enrollment_blocks = enrollment_blocks( ...
    ~ismember(enrollment_blocks,testing_block) );

% process enrollmenta and testing data!
for s=1:subject_count
    for c=1:channel_count
        index = (s-1)*channel_count + (c-1) + 1;
        [enrollment_data{s,c}, testing_data{s,c}] = eegTestTrainSplit( ...
            full_index{index}, features, block_size, window_overlap, ...
            enrollment_blocks, testing_block);
    end
end

% training_data is built from a different set of files, ignore epochs as
% the ubm won't care about them!
end