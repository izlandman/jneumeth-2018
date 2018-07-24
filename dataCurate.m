% split data from listAgg.
% given input_data and epoch parameters in epoch_listing, produce a
% training and testing dataset. i suppose this should be run each time
% epoch_listing changes

% input_data is a cell!
%   size SUBJECT x SESSION
%   cell elements -> FEATURES x CHANNELS x EPOCHS

% epoch_listing is a cell!
%   {1} : level of permutation | epochs, channels, sessions, subjects
%   {2} : -1 (random number) OR [y...z] set of test index
function [training_data, testing_data] = ...
    dataCurate(input_data, element_listing)

[n_subject, n_session] = size(input_data);
% find minimum number of epochs per session
full_subject = n_subject * n_session;
min_epoch = zeros(n_session,1);
rand_epochs = zeros(full_subject, element_listing);
for ses=1:n_session
    min_epoch(ses) = min(min(cell2mat(cellfun(@size, input_data(:,ses), ...
        'UniformOutput', false)),[],2));
    epoch_init = (ses-1)*n_subject + 1;
    epoch_fin = epoch_init + n_subject - 1;
    rand_epochs (epoch_init:epoch_fin,:) = repmat(randperm( ...
        min_epoch(ses),element_listing), n_subject, 1 );
end

% setup output variables, combine subject X epoch size!
training_data = cell(full_subject, 1);
testing_data = cell(full_subject, element_listing);

for ses=1:n_session
    for sub=1:n_subject
%         fprintf('Building data for sub: %d ses: %d\n', sub, ses);
%         if( sub == 45 )
%             fprintf('Error here!');
%         end
        full_index = (sub-1) + (ses-1)*n_subject + 1;
        test_index = rand_epochs(full_index,:);
        [features, channels, epochs] = size(input_data{sub,ses});
        % check if all testing_epochs are within data space
        all_index = 1:epochs;
        train_index = all_index(~ismember(all_index, ...
            test_index));
        training_data{full_index} = reshape(...
            input_data{sub,ses}(:,:,train_index), features, channels * ...
            numel(train_index));
        for e=1:element_listing
        testing_data{full_index,e} = reshape(...
            input_data{sub,ses}(:,:,test_index(e)), features, channels * ...
            1);
        end
    end
end

end
