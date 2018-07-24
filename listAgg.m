% produce a function that takes in N file lists and aggregates them in a
% way that can be used for producing UBMs

function [output_data, index_str, index_list] = listAgg(varargin)

% sort out how many lists are being given and what to do with them, as HTK
% files must be split while COH/PSD come in their epoch based format!

n_inputs = nargin;
if( n_inputs == 0 )
    fprintf('No arguments passed!\n');
    return;
elseif( ischar(varargin{1}) )
    % first arguement is a character, making this PSD or COH files
    fprintf('Pass in %d arguments for PSD/COH.\n',nargin);
    % read from the files!
    for i=1:n_inputs
        [file_index, subject_count, session_count, full_index] ...
            = fileListConvert(varargin{i});
        if( i == 1 )
            index_full = full_index;
            index_epoch = cell2mat(file_index(:,3));
            epoch_sort = file_index(:,4);
        else
            index_full = cat(1,index_full, full_index);
            index_epoch = cat(1,index_epoch, cell2mat(file_index(:,3)));
            epoch_sort = cat(1,epoch_sort,file_index(:,4));
        end
    end
    % produce final epoch listing
    [e_label, e_sort] = sort(epoch_sort);
    index_epoch = index_epoch(e_sort);
    epoch_tags = zeros(numel(index_epoch),2);
    epoch_tags(:,1) = 1+[0 cumsum(index_epoch(1:end-1))'];
    epoch_tags(:,2) = cumsum(index_epoch);
    % reorder according to files names! (they better be correct!)
    [index_str, order_ind] = sort(index_full);
    % parse names to make index listing!
    n_files = numel(index_str);
    % subject by session by epoch!
    index_list = zeros(n_files,3);
    for i=1:n_files
        [~,f_name,~] = fileparts(index_str{i});
        split_name = strsplit(f_name,'_');
        channel = strsplit(split_name{3},'h');
        index_list(i,1:3) = [ str2double(split_name{1}), ...
            str2double(split_name{2}) ,str2double(channel{2})];
    end
    
    n_subject = numel(unique(index_list(:,1)));
    n_session = numel(unique(index_list(:,2)));
    data_peek = iVectorBinary(index_str{1});
    [n_feature, n_channel] = size(data_peek);
    % now read the data into a very large cell array?
    output_data = cell(n_subject,n_session);
    parfor sub=1:n_subject
        for ses=1:n_session
            l_index = (sub-1)*n_session + ses;
            n_epoch = epoch_tags(l_index,2)-epoch_tags(l_index,1)+1;
            data_peek = zeros(n_feature,n_channel,n_epoch);
            e_index = epoch_tags(l_index,1):epoch_tags(l_index,2);
            for e=1:n_epoch
                data_peek(:,:,e) = iVectorBinary(index_str{e_index(e)});
            end
            output_data{sub,ses} = data_peek;
        end
    end
    
elseif( isnumeric(varargin{1}) && n_inputs > 2 )
    % first argument is a number, making this HTK files
    fprintf('Pass in %d arguments for HTK.\n',nargin);
    block_size = varargin{1};
    % should be %!
    window_overlap = varargin{2};
    start_offset = ceil(block_size/2);
    step_offset = ceil(block_size * window_overlap);
    for i=3:n_inputs
        [file_index, subject_count, session_count, full_index] ...
            = fileListConvert(varargin{i});
        if( i == 3 )
            index_full = full_index;
            index_channel = cell2mat(file_index(:,3));
            channel_sort = file_index(:,4);
        else
            index_full = cat(1,index_full, full_index);
            index_channel = cat(1,index_channel, cell2mat(file_index(:,3)));
            channel_sort = cat(1,channel_sort,file_index(:,4));
        end
    end
    % produce final epoch listing
    [c_label, c_sort] = sort(channel_sort);
    index_channel = index_channel(c_sort);
    channel_tags = zeros(numel(index_channel),2);
    channel_tags(:,1) = 1+[0 cumsum(index_channel(1:end-1))'];
    channel_tags(:,2) = cumsum(index_channel);
    % reorder according to files names! (they better be correct!)
    [index_str, order_ind] = sort(index_full);
    n_files = numel(index_str);
    % subject by session by epoch!
    index_list = zeros(n_files,3);
    for i=1:n_files
        [~,f_name,~] = fileparts(index_str{i});
        split_name = strsplit(f_name,'_');
        channel = strsplit(split_name{3},'h');
        index_list(i,1:3) = [ str2double(split_name{1}), ...
            str2double(split_name{2}) ,1+str2double(channel{2})];
    end
    n_subject = numel(unique(index_list(:,1)));
    n_session = numel(unique(index_list(:,2)));
    data_peek = htkread(index_str{1});
    [n_feature, ~] = size(data_peek);
    % now read the data into a very large cell array?
    output_data = cell(n_subject,n_session);
    parfor sub=1:n_subject
        for ses=1:n_session
            l_index = (sub-1)*n_session + ses;
            n_channel = channel_tags(l_index,2)-channel_tags(l_index,1)+1;
            c_index = channel_tags(l_index,1):channel_tags(l_index,2);
            data_peek = htkread(index_str{c_index(1)});
            [~, n_sample] = size(data_peek);
            % build epochs on the fly!
            blocks = start_offset:(block_size-step_offset):n_sample;
            n_blocks = numel(blocks);
            data_peek = zeros(n_feature,n_channel,n_blocks);
            for c=1:n_channel
                all_data = htkread(index_str{c_index(c)});
                data_peek(:,c,:) = all_data(:,blocks);
            end
            % ensure the set has the same number of epochs!
            output_data{sub,ses} = data_peek;
        end
    end
    
end
end