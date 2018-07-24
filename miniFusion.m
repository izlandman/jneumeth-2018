% split up the fusion-match scoring procedure
function [f_CCR, f_mCCR, CCR_bench]=miniFusion(dist,mCCR_val, mCCR_ind,...
    epochs,fusion_index,subject_count,eval_flag)

f_CCR = zeros(epochs,1);
% channels will be added to chan_fusion as they improve algorithm
% performance, start with the 'best' individual element
chan_fusion = mCCR_ind(1);
channel_count = numel(fusion_index);
CCR_bench = zeros(channel_count,4);
CCR_index = 1;
CCR_bench(1,:) = mCCR_val(1);
CCR_bench(1,[2 4]) = chan_fusion;
if( eval_flag == 1 )
    % check every channel
    for chan_ind=2: channel_count
        new_channel = mCCR_ind(chan_ind);
        valid_channels = cat(2,chan_fusion,new_channel);
        parfor cv_step=1:epochs
            fusion_dist = zeros(subject_count,subject_count);
            for s=1:subject_count
                for t=1:subject_count
                    fusion_dist(s,t) = sum(1./dist(s, t, ...
                        valid_channels, cv_step), 3);
                end
            end
            [~,f_index] = sort(fusion_dist,'descend');
            f_CCR(cv_step) = (1/subject_count) * ...
                sum( [1:subject_count] == f_index(1,:) );
        end
        [f_mCCR, CCR_bench, chan_fusion, CCR_index] = ...
            ccrTrack(epochs, f_CCR, new_channel, CCR_index, ...
            CCR_bench, valid_channels, chan_ind, chan_fusion);
    end
else
    % check every channel
    for chan_ind=2: channel_count
        new_channel = mCCR_ind(chan_ind);
        valid_channels = cat(2,chan_fusion,new_channel);
        parfor cv_step=1:epochs
            fusion_dist = zeros(subject_count,subject_count);
            for s=1:subject_count
                for t=1:subject_count
                    fusion_dist(s,t) = sum(dist(s, t, ...
                        valid_channels, cv_step), 3);
                end
            end
            [~,f_index] = sort(fusion_dist,'descend');
            f_CCR(cv_step) = (1/subject_count) * ...
                sum( [1:subject_count] == f_index(1,:) );
        end
        [f_mCCR, CCR_bench, chan_fusion, CCR_index] = ...
            ccrTrack(epochs, f_CCR, new_channel, CCR_index, ...
            CCR_bench, valid_channels, chan_ind, chan_fusion);
    end
end
end

function [f_mCCR, CCR_bench, chan_fusion, CCR_index] = ...
    ccrTrack(epochs, CCR, channel, ccr_index, ccr_bench, ...
    v_channels, chan_ind, c_fusion)

f_mCCR = (1/epochs)*sum(CCR);
ccr_bench(chan_ind,3) = f_mCCR;
ccr_bench(chan_ind,4) = channel;

if( f_mCCR > ccr_bench(ccr_index,1) )
    ccr_index = ccr_index + 1;
    c_fusion = v_channels;
    ccr_bench(ccr_index,1) = f_mCCR;
    ccr_bench(ccr_index,2) = channel;
end

CCR_bench = ccr_bench;
chan_fusion = c_fusion;
CCR_index = ccr_index;

end