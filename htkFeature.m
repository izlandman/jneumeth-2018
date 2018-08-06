function htk_data = htkFeature(file)

% read the edf file!
[header, data] = edfread(file);

%% Parameter Organization

% setup montage channels
montage_channels = [22 30; 30 41; 41 47; 47 61; 24 38; 38 42; 42 55; ...
    55 63; 43 41; 41 9; 9 11; 11 13; 13 42; 42 44; 22 32; 32 9; 9 49; ...
    49 61; 24 36; 36 13; 13 53; 53 63];
n_channel = numel(unique(montage_channels(:)));
chan_list = unique(montage_channels(:));
% interp channels as A1 (43) and A2 (44) do not exist in physioNet
A1 = [41 47 30];
A2 = [42 55 38];
data(43,:) = mean( data(A1,:) );
data(44,:) = mean( data(A2,:) );

% remove empty columns of data
data = data(:,sum(data)~=0);

% mel options, ensure linear frequency scaling, hence the f
mel_options = 'fm';
min_freq = 0.0031;
max_freq = 0.1563;

% mel parameters
sam_per_sec = 10;
sample_rate = max(header.samples);
mel_filters = 8;
mel_feat = ((mel_filters+1)*3)-1;
window = 32;
frame = 16;
nfft = 512;
% bonus feature params
d_win = 9;
dd_win = 3;

%% Build HTK data
duration = round(size(data,2)/max(header.samples)*sam_per_sec);
E_d = zeros(1,duration);
delta = zeros(mel_filters+1,duration);
delta_delta = zeros(mel_filters+1,duration);
output_data = zeros( mel_feat, duration, n_channel );
for c=1:n_channel
    % pad leading and trailing half-windows with data
    [cepstrum_data,~,~,E_f] = myCEP([data(chan_list(c),1:frame/2) ...
        data(chan_list(c),:) data(chan_list(c),end-frame+1:end)], ...
        sample_rate,window,frame,nfft,mel_filters,mel_options,min_freq,...
        max_freq);
    % replace 1st with peak energy, 0th isn't returned!
    output_data(1,:,c) = E_f;
    output_data(2:mel_filters,:,c) = cepstrum_data(2:end,:);
    % produce max diff between E_f
    for i=1:duration
        % windowReturn provides filler zeros to pad ends, remove these!
        E_d(i) = max(windowReturn(E_f,d_win,i,0)) - ...
            min(windowReturn(E_f,d_win,i,0));
    end
    output_data(mel_filters+1,:,c) = E_d;
    % build deltas
    for i=1:duration
        delta(:,i) = cepDelta(output_data(1:mel_filters+1,:,c),d_win,i,1);
    end
    % build delta-deltas
    for i=1:duration
        delta_delta(:,i) = cepDelta(delta,dd_win,i,1);
    end
    % build cepstrum channel data, exclude 0th terms and replace with
    % delta, and delta-delta resulting in 9*3-3-2=26 features
    output_data(10:18,:,c) = delta;
    output_data(19:end,:,c) = delta_delta([1:mel_filters],:);
end

%% Package and write out as .HTK
% build montage
htk_data = zeros(mel_feat,duration,n_channel);
for c=1:n_channel
    htk_data(:,:,c) = output_data(:,:,chan_list==montage_channels(c,1))-...
        output_data(:,:,chan_list==montage_channels(c,2));
end

end

function result = cepDelta(data,window,index,filler)

result = windowReturn(data,window,index,filler);

numer = sum(result,2);
denom = 2*sum(result(:,(window-1)/2+1).^2,2);

result = numer ./ denom;

end

function result = windowReturn(data,window,index,filler)
% return data centered around frame index
fill = 0;
middle = index;
start = middle - floor((window-1)/2) - 1;
finish = middle + floor((window-1)/2);
if( start < 1 )
    start = 1;
    fill = -1;
end
if( finish > size(data,2) )
    finish = size(data,2);
    fill = 1;
end
% zero pad if necessary?
if( filler )
    missing = window - size(data(:,start:finish),2);
    if( fill == -1 )
        result = [zeros(size(data,1),missing) data(:,start:finish)];
    elseif( fill == 1)
        result = [data(:,start:finish) zeros(size(data,1),missing)];
    else
        result = data(:,start:finish);
    end
else
    result = data(:,start:finish);
end

end