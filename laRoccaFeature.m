% turn the raw edf data into the features used in the 2014 La Rocca paper.
% this requires them to be built into 10 second epochs from which PSD
% features and spectral coherence connectivity features are drawn. they onl
% use 56 of the 64 electrodes!

% the returned 4D matarix can be used to build an averaged PSD for the
% epoch OR compute the COH for the each channel in the epoch. this requires
% a function wrapper on this to address the type of feature to be built and
% how that feature is saved.

% pass in a file_name and a list of channels to process.
% !!! Be sure true_channels includes channels 43 and 44 to follow 

function data_100 = laRoccaFeature(file,true_channels,mean_channels,...
    epoch_dur, n_epoch)

% read the edf file!
[header, data] = edfread(file);

% channel parameters
n_channels = numel(true_channels);

% they convery them from 160 samples per second to 100 samples per second.
% up sample, then downsample!

[~,samples] = size(data(1,:));
t_epoch_size = epoch_dur * header.samples(1);
t_epoch_n = samples / t_epoch_size;
% is there enough data? otherwise copy and paste!
if( t_epoch_n < n_epoch )
    new_data = repmat(data,1,ceil(n_epoch/t_epoch_n));
    data = new_data(:,1:t_epoch_size*n_epoch);
end
new_sample = 100;
epoch_size = epoch_dur * new_sample;
data_100 = zeros(epoch_size,n_channels,n_epoch);

% build out the PSD based on the hanning window samples. the returned
% matrix is 4D, containing the windows for each channel/frequency
% combination for each epoch.

% filter before resampling?
% fc = 40;
% fs = 100;
% [b,a] = butter(6,fc/(fs/2));
% harder edge!
[b,a] = ellip(6,5,40,49/100*2);
data_f = filter(b,a,data);
data_down = resample(data_f',5,8);

% allow for coherence normalization to 'mean' channels
if( mean_channels ~= 0 )
    data_mean = mean(data_down(:,mean_channels),2);
else
    data_mean = 0;
end
data_down = data_down - data_mean;

for c=1:n_channels
    % La Rocca said use resample!
    % data_up = interp(data(c,:),5);
    % data_down = decimate(data_up,8);
    chan_length = numel(data_down(:,true_channels(c)));
    for e=1:n_epoch
        index_s = (e-1)*epoch_size + 1;
        index_f = e*epoch_size;
        index_v = index_s:index_f;
        if( index_f > chan_length )
            index_f = chan_length;
            index_s = index_f-epoch_size+1;
            index_v = index_s:index_f;
        end
        data_100(:,c,e) = data_down(index_v,true_channels(c));
    end
end

end