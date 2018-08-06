% write the desired features to files. read each file and write out channel
% data to new files for the given features

% feature_type: 0 for PSD | 1 for COH | 2 for CEP

function laRoccaFeatureBuilder(file_list,feature_type,final_folder,...
    dur_epoch,n_epochs,workers)

% enhanced debug feature
% setSchedulerMessageHandler(@disp);
% clear existing cluster
delete(gcp('nocreate'))
% setup parallel environment
cluster = parcluster('local');
cluster.NumWorkers = workers;
% this needs to be ON for cluster operation and OFF for local debugging
if isunix
    cluster.JobStorageLocation = getenv('HOME');
end
parpool(cluster,workers);

[~,name,~] = fileparts(file_list);
fid = fopen(file_list);
txt = textscan(fid,'%s');
fclose(fid);
file_list = txt{1};
n_files = numel(file_list);

% break out each feature generation
if ( feature_type == 2 )
    fprintf('Building CEP features.\n');
    % build out from the already existing cepstrum features from which
    % there are only 22 channels
    name = name(end-2:end);
    extn = 'CEP';
    data_peek = htkread(file_list{1});
    [dp_feats, dp_samples] = size(data_peek);
    n_channels = 22;
    % 10 htks per second of data, offset being 'middle' frame
    dur_window = 10;
    offset_window = 5;
    n_subjects = n_files/n_channels;
    % resolve floor of samples per epoch. htk data is already stored as
    % features per tenth second. will have to average samples over range to
    % produce final features. use FULL window of each second, no overlaps!
    epoch_samples = floor(dp_samples/n_epochs/dur_window);
    e_valid = reshape(1:n_epochs*epoch_samples,epoch_samples,[]);
    e_map = (e_valid-1)*dur_window+offset_window;
    file_list = reshape(file_list,n_channels,n_subjects);
    output_data = zeros(dp_feats,n_channels,n_epochs,n_subjects);
    parfor n=1:n_subjects
        data_temp = zeros(dp_feats,n_channels,n_epochs);
        fprintf('Building CEP features: %d of %d.\n',n,n_subjects);
        for c=1:n_channels
            data_peek = htkread(file_list{c,n});
            % check for sample underfill!
            underfill = max(e_map(:))/size(data_peek,2);
            if( underfill > 1 )
                data_peek = repmat(data_peek,1,ceil(underfill));
                data_peek = data_peek(:,1:max(e_map(:)));
            end
            data_shift = reshape(...
                data_peek(:,(e_valid-1)*dur_window+offset_window),...
                dp_feats,epoch_samples,[]);
            data_temp(:,c,:) = squeeze(mean(data_shift,2));
        end
        output_data(:,:,:,n) = data_temp;
    end
elseif( sum(feature_type == [0 1 3 4] == 1 ) )
    
    fprintf('Building COH features.\n');
    % Rocca Filter Prep
    hann_win = hann(100);
    overlap = 50;
    nfft = 100;
    index = 2:41;
    feats = numel(index);
    n_subjects = n_files;
    
    % peek at data, with channel list
    if( feature_type == 0 || feature_type == 1 )
        % LaRocca's channel list
        valid_channels = 1:64;
        remove_channels = [25,29,39,40,43,44,61,64];
        channel_list = valid_channels(~ismember(valid_channels,...
            remove_channels));
        mean_channels = [43, 44];
    elseif( feature_type == 3 || feature_type == 4)
        % TUH EEG TCP Montage channel list
        channel_list = [9,11,13,22,24,30,32,36,38,41,42,43,44,47,49,53,...
            55,61,63];
        mean_channels = 0;
    else
        fprintf('Incorrect feature_type value.\n');
    end
    data_peek = laRoccaFeature(file_list{1},channel_list,mean_channels,...
        dur_epoch,n_epochs);
    [dp_feat, dp_channels, dp_epochs] = size(data_peek);
    name = name(end-2:end);
    
    % unique for COH and PSD
    if( feature_type == 1 || feature_type == 4 )
        extn = 'COH';
        dp_channels = ( dp_channels * (dp_channels-1) )/2;
        output_data = zeros(feats,dp_channels,n_epochs,n_subjects);
        parfor n=1:n_subjects
            channel_data = laRoccaFeature(file_list{n},channel_list,...
                mean_channels,dur_epoch,n_epochs);
            [samples, channels, epochs] = size(channel_data);
            % apparently some data doesn't have enough of  a recording?
            diff_epoch = n_epochs  - epochs;
            if( diff_epoch > 0 )
                cat_data = channel_data(:,:,end);
                for z=1:diff_epoch
                    channel_data = cat(3,channel_data,cat_data);
                end
            end
            % COH features
            % tic
            fprintf('Building COH features: %d of %d.\n',n,n_subjects);
            data_temp = zeros(feats,dp_channels,n_epochs);
            for e=1:n_epochs
                % they scale them with atan in their work
                data_temp(:,:,e) = createCOH(channel_data(:,:,e),...
                    hann_win,overlap,nfft,index);
            end
            % this represents all possible combinations of coherence,
            % however for case 4 we only want those relationships that
            % mimic the TCP layout of TUH EEG
            output_data(:,:,:,n) = data_temp;
        end
        % thus it must be pruned
        if( feature_type == 4 )
            TCP_montage = [53 84 130 160 70 118 141 170 128 9 1 19 43 ...
                137 54 6 14 164 69 40 48 168];
            output_data = output_data(:,TCP_montage,:,:);
        end
    elseif( feature_type == 0 || feature_type == 3 )
        extn = 'PSD';
        output_data = zeros(feats,dp_channels,n_epochs,n_subjects);
        for n=1:n_subjects
            channel_data = laRoccaFeature(file_list{n},channel_list,...
                mean_channels,dur_epoch,n_epochs);
            [samples, channels, epochs] = size(channel_data);
            % apparently some data doesn't have enough of  a recording?
            diff_epoch = n_epochs  - epochs;
            if( diff_epoch > 0 )
                cat_data = channel_data(:,:,end);
                for z=1:diff_epoch
                    channel_data = cat(3,channel_data,cat_data);
                end
            end
            % PSD features, average the data over each epoch
            % tic
            fprintf('Building PSD features: %d of %d.\n',n,n_subjects);
            data_temp = zeros(feats, channels, n_epochs);
            for e=1:n_epochs
                for c=1:channels
                    % perhaps this should be pwelch not spectrogram
                    data_full = pwelch(channel_data(:,c,e), ...
                        hann_win, overlap, nfft);
                    data_temp(:,c,e) = data_full(index);
                end
            end
            output_data(:,:,:,n) = log10(data_temp);
        end
    else
    end
else
    fprintf('Wrong feature_type selection.\n');
end

% save file listings
listing = [final_folder filesep extn];
directoryCheck(listing);
listing = [final_folder filesep extn filesep 'features_' name '.list'];
fid = fopen(listing,'w');

for n=1:n_subjects
    % write channels to file! use the binary tool! Each 'channel' is given
    % a file which is a bit silly, but maintains the consistency of their
    % experiment and how the I-Vector/GMM kits presently function.
    % tic
    [features, channels, epochs] = size(output_data(:,:,:,n));
    if( n_files > n_subjects )
        [path, name, ext] = fileparts(file_list{1,n});
        % read htks, different filenames
        sub_fldr = path(end-3:end);
        sub = name(1:3);
        ses = name(5:6);
    else
        [path, name, ext] = fileparts(file_list{n});
        sub_fldr = name(1:4);
        sub = name(2:4);
        ses = name(6:7);
    end
    sub_save = [final_folder filesep extn filesep sub_fldr];
    % open output file
    directoryCheck(sub_save);
    for e=1:epochs
        data_out = squeeze(output_data(:,:,e,n));
        file_out = [sub_save filesep sub '_' ses '_' 'epoch' ...
            sprintf('%04d',e) '.bin'];
        iVectorBinary(file_out,data_out);
        % write file listing to new file_list
        fprintf(fid,'%s\n',file_out);
    end
    % a1 = toc;
    % fprintf('Wrote subject %s session %s in %f seconds.\n',sub,ses,a1);
end
fclose(fid);

% make plots?
% stats_file = [final_folder filesep extn];
% laRoccaDataProfile(output_data,feature_type,stats_file);
end