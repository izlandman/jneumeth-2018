% take file name and use file list to generate cell listing of edf files
% for processing in the programs i already wrote

% subject_loc  -  indicate where to find the subject ID in the string of
% the filename

function [ordered_result, subject_count, session_total, result_copy] ...
    = fileListConvert(file_name)

try
    fprintf('Attempting to open file: %s\n', file_name);
    open_file = fopen(file_name);
        % use textscan to grab all filenames
    result = textscan(open_file,'%s');
    result = result{:};
    count = numel(result);
catch
    fprintf('FOPEN failed!\n');
    if( iscell(file_name) )
        count = numel(file_name);
        fprintf('cell count: %d\n', count);
        if( count > 1 )
            result = file_name;
        else
            result = file_name{:};
            count = numel(result);
        end
    else
        fprintf('The passed variable is not a cell!\n');
        result = -1;
        count = -1;
    end
end

if( count == -1 )
    fprintf('No data returned from the file list. Abort.\n');
else
    
    subject_tags = zeros(count,1);
    session_tags = zeros(count,1);
    display(['file count: ' num2str(count)]);
    
    for i=1:count
        % means new file type, physionet type file structure
        [f_path, f_name, ~] = fileparts(result{i});
        % assume the last folder is the subject ID, filsep
        folder = strsplit(f_path,'/');
        if( iscell(folder)~= 1)
            folder = strsplit(f_path,'\');
        end
        if( numel(folder) == 1 )
            folder = strsplit(f_path,'\');
        end
        % check to see what type of file path that needs to be parsed!
        if( ~isempty( strfind(result{i},'/data/isip/data/') ) )
            % means it is of type old file type
            parts = strsplit(f_name,'_');
            % update subject_tag
            subject_tags(i) = str2double(folder(end));
            % check for leading 'a', which ruins everything
            if( strcmp(parts(1),'a' ) )
                parts = strsplit(f_path,'/');
                % find the section with 'sXX_' use fancy strfind
                parts_index = find(~cellfun(@isempty,...
                    strfind(parts,'s')),1,'last');
                subject_tags(i) = str2double(folder(parts_index-1));
                % find correct session ID, check betwen a_x_
                session_split = strsplit(f_name,'_');
                if( numel(session_split) == 2 )
                    session_tag_a = 00;
                else
                    session_tag_a = str2double(session_split{end-1});
                end
                % then find session tag from previous folder!
                session_tag_s = strsplit(parts{parts_index},'_');
                session_tag_s = strsplit(session_tag_s{1},'s');
                session_tag_s = str2double(session_tag_s{end});
                % now cat them and turn them into the true session_id
                session_tags(i) = session_tag_s + session_tag_a/1000;
            elseif( isnan(str2double( parts(end-1)) ) )
                session_tags(i) = 0;
            else
                session_tags(i) = str2double( parts(end-1) );
            end
        else
            % means new file type, physionet type file structure
            parts = strsplit(f_name,'_');
            % filter out any letters from subject string
            sub_id = folder(end);
            subject_tags(i) = str2double(sub_id);
            if( isnan( subject_tags(i) ) )
                temp_tag = strsplit( sub_id{1}, 'S' );
                if( strcmp(temp_tag{1},sub_id{1}) )
                    temp_tag = strsplit( sub_id{1}, 's' );
                    temp_tag = strsplit( temp_tag{2}, '_');
                    temp_tag = strjoin(temp_tag);
                    temp_tag = temp_tag(find(~isspace(temp_tag)));
                    subject_tags(i) = str2double( temp_tag );
                else
                    subject_tags(i) = str2double( temp_tag{2} );
                end
            end
            if( sum( isnan(str2double(parts{2})) ) > 0 )
                session_tags(i) = subject_tags(i);
            else
                session_tags(i) = str2double(parts{2});
            end
        end
    end
    
    % write new cell array that stores the grouped file names by subject
    % and session to make it easier to build models. first figure out how
    % many unique entries in the data are present
    ids = unique(subject_tags);
    subject_count = numel(ids);
    session_uni = cell(subject_count,1);
    session_count = zeros(subject_count,1);
    channels = zeros(count,1);
    
    % there seems to be an error with the sorting of the channels causing
    % the 000 channel to slide between channel groups. to fix this sort
    % each subject's sessions and restructure the input file list.
    result_copy = result;
    session_copy = session_tags;
    for i=1:subject_count
        sub_index = (subject_tags == i);
        [~,index] = sort(session_tags(sub_index),'ascend');
        result_copy(min(index):max(index)) = result(index);
        session_copy(min(index):max(index)) = session_tags(index);
    end
    session_tags = session_copy;
    
    % find session count for each subject
    offset = 0;
    for i=1:subject_count
        ind =  subject_tags == ids(i) ;
        session_uni{i} = unique(session_tags(ind));
        session_count(i) = numel(session_uni{i});
        for r=1:session_count(i)
            offset = 1 + offset;
            channels(offset) = sum(session_uni{i}(r) == session_tags(ind));
        end
        session_count(i) = numel(session_uni{i});
    end
    session_total = sum(session_count);
    % prune empty entries
    channels = channels(1:session_total);
    
    channel_based_index = [0 cumsum(channels(1:end-1))'] + 1;
    result_temp = cell(length(channel_based_index),4);
    
    if( sum( channels/min(channels) == channels/max(channels) )== session_total )
        channel_count = max(channels);
        fprintf('Channels are uniform for all subjects: %d\n',channel_count);
        % order files, THIS IS NO GOOD. WHY DOES IT SORT THEM AND LOSE
        % TRACK OF THE TRUE FILENAMES! BAD BAD BAD! FIX!
        for i=1:length(channel_based_index)
            result_temp{i,1} = subject_tags(channel_based_index(i));
            result_temp{i,2} = session_tags(channel_based_index(i));
            result_temp{i,3} = channels(i);
            result_temp{i,4} = result{channel_based_index(i)};
        end
    else
        min_channel = min(channels);
        max_channel = max(channels);
        fprintf('Channels are not uniform! Min: %d Max: %d\n',...
            min_channel, max_channel);
        for i=1:length(channel_based_index)
            result_temp{i,1} = subject_tags(channel_based_index(i));
            result_temp{i,2} = session_tags(channel_based_index(i));
            result_temp{i,3} = channels(i);
            result_temp{i,4} = result{channel_based_index(i)};
        end
    end
    ordered_result = result_temp;
end

end
