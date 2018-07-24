% singular fusion matching function

function [f_CCR, f_mCCR, CCR_bench, a1] = fusionMatchScoring(varargin)

% define variables from those passed in
fusion_flag = varargin{1};
mCCR = varargin{2};
epochs = varargin{3};
subject_count = varargin{4};
dist = varargin{5};
save_folder = varargin{6};
eval_flag = 1;
if( nargin == 7 )
    n_mixtures = numel(varargin{7});
    eval_flag = 0;
end

%% enable feature fusion algorithm
if( fusion_flag ~= 0 )
    tic
    if( eval_flag == 1 )
        % distances from Mahal can be scaled like this because smaller
        % values indicate a better match
        fprintf('Fusion classification beings!\n');
        [mCCR_val, mCCR_ind, fusion_index] = laRoccaFusionIndex(...
            fusion_flag, mCCR);
        [f_CCR, f_mCCR, CCR_bench] = miniFusion(dist,mCCR_val,mCCR_ind,...
            epochs,fusion_index,subject_count,eval_flag);
    else
        % handle those producing gmm and i-vector scores
        fprintf('Fusion (mixture) classification beings!\n');
        [channel_count,~,~] = size(mCCR);
        f_CCR = zeros(epochs, n_mixtures);
        f_mCCR = zeros(n_mixtures, 1);
        % depth of 4 to produce mean CCR, best CCR, ideal CCR, CCR index
        CCR_bench = zeros(channel_count, 4, n_mixtures);
        parfor m=1:n_mixtures
            tic
            ops_data = squeeze(dist(:,:,m,:,:));
            [mCCR_val, mCCR_ind, fusion_index] = laRoccaFusionIndex(...
                fusion_flag, mCCR(:,m));
            [f_CCR(:,m), f_mCCR(m), CCR_bench(:,:,m)] = miniFusion(...
                ops_data,mCCR_val,mCCR_ind,epochs,fusion_index,...
                subject_count,eval_flag);
            a0 = toc;
            fprintf('Resolved fusion of mixture size %d in %f seconds\n',...
                m, a0);
        end
    end
    % write results to file
    iVectorBinary([save_folder filesep 'f_mCCR.bin'], f_mCCR);
    iVectorBinary([save_folder filesep 'f_CCR.bin'], f_CCR);
    iVectorBinary([save_folder filesep 'CCR_bench.bin'], CCR_bench);
    a1 = toc;
    fprintf('Fusion classificaiton complete in %f seconds.\n',a1);
else
    tic
    f_mCCR = -1;
    f_CCR = -1;
    CCR_bench = -1;
    a1 = toc;
end
end