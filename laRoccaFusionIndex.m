% curate an appropriate index of the channels/elements
function [mCCR_val, mCCR_ind, fusion_index] = ...
    laRoccaFusionIndex(fusion_flag, mCCR)

n_mCCR = numel(mCCR);

% they  run three sections of features, enable a way to prune index
% listing of mCCR
if( fusion_flag == 1 )
    fusion_index = 1:n_mCCR;
elseif( fusion_flag == 2 )
    % Frontal Channels
    fusion_index = [22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, ...
        34, 35, 36];
elseif( fusion_flag == 3 )
    % Central Channels
    fusion_index = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, ...
        15, 16, 17, 18, 19, 20, 21, 37, 38, 39, 40];
elseif( fusion_flag == 4 )
    % Parieto-occipital Channels
    fusion_index = [41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, ...
        52, 53, 54, 55, 56];
end
% deal with the 1540 COH channels!
if( n_mCCR == 1540 )
    coh_index = zeros(56,56,2);
    coh_index(:,:,1) = repmat([1:56],56,1);
    coh_index(:,:,2) = repmat([1:56]',1,56);
    coh_struct = zeros(56*56,2);
    for i=1:56
        for j=i+1:56
            coh_struct( (i-1)*56+j,: ) = coh_index(i,j,:);
        end
    end
    coh_struct = coh_struct(coh_struct(:,1)~=0,:);
    [ind_i,~] = ismember(coh_struct(:,1), fusion_index);
    [ind_j,~] = ismember(coh_struct(:,2), fusion_index);
    fusion_index = 1:1540;
    fusion_index = fusion_index(logical(ind_i.*ind_j));
end

[mCCR_val, mCCR_ind] = sort(mCCR(fusion_index),'descend');

end