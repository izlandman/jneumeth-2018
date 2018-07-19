% build a class for the TVM, defined as type 'handle' means only one can
% exist at a time

classdef TotalVariabilityMatrix < handle
    properties
        % store depth of matrix
        Depth
        % store iteration limit
        Iterations
        % Should the output be verbose?
        Verbosity
        % Provide an output file name
        FileName
        % store the matrix condition number
        MatCond
        % store the number of mixtures
        Mixtures
        % store the number of features
        Features
        % this is a structure with three types. mu Sigma w
        UBM
        % UBM could be a cell array contain multiple models~
        n_UBM
        % this is the actual matrix
        TVM
        % store error between the mixtures
        Error
        % TVM sets
        Sets
        % TVM sets threshol
        SetsThreshold
        % store column based mse results for each mixture
        fullMSE
        % store the 'learning rate' of new/old TVM
        LearningRate
        % store similarity based upon normalized row weights
        RowSim
        % LDA transform?
        V
        % and
        D
        % store the LDA norm'd TVM
        TVM_norm
        % store the sim of TVM_norm
        RowSim_norm
        % store correct recognition rate
        crr
        % store equal error rate
        eer
        % channels
        n_Channels
        % fnr values
        fnr
        % fpr values
        fpr
    end
    methods
        % instantiate the object by providing a universal background model
        function obj = TotalVariabilityMatrix(ubm,depth,iterations,...
                channels)
            % store depth (size) of ivs
            obj.Depth = depth;
            % store iterations of training process
            obj.Iterations = iterations;
            % store number of channels to be processed
            obj.n_Channels = channels;
            % allow more than a single model to be passed into the class!
            obj.n_UBM = numel(ubm);
            obj.UBM = ubm;
            % determine each mixture size
            obj.Mixtures = zeros(obj.n_UBM,1);
            for i=1:obj.n_UBM
                obj.Mixtures(i) = numel(ubm{i}.w);
            end
            % assume all mixtures have equal number of features
            obj.Features = size(ubm{1}.mu,1);
            % default these settings to off
            obj.Verbosity = 0;
            obj.FileName = 0;
            % TVM must be a cell array  to handle various sized matrices
            % based upon the size of the mixtures
            obj.TVM = cell(obj.n_UBM,obj.n_Channels);
            % set threshold
            obj.SetsThreshold = 0.89;
            % condition matrix
            obj.MatCond = zeros(obj.n_UBM,obj.n_Channels);
            % full mse results
            obj.fullMSE = cell(obj.n_UBM,obj.n_Channels);
            % setup learning rate cell, full matrix and sub-matrices
            obj.LearningRate = cell(obj.n_UBM,2,obj.n_Channels);
            % setup Row Sim results
            obj.RowSim = cell(obj.n_UBM,obj.n_Channels);
            % setup V LDA transformation matrix
            obj.V = cell(obj.n_UBM,obj.n_Channels);
            % setup D weights vector
            obj.D = cell(obj.n_UBM,obj.n_Channels);
            % store the transformed TVM
            obj.TVM_norm = cell(obj.n_UBM,obj.n_Channels);
            % store the transformed TVM similarity
            obj.RowSim_norm = cell(obj.n_UBM,obj.n_Channels);
            % store CRR results
            obj.crr = zeros(obj.n_Channels,obj.n_UBM);
            % store EER results
            obj.eer = zeros(obj.n_Channels,obj.n_UBM);
        end
        % produce the TVM
        function train(tvm,data,index)
            % enable all training if index == 0
            if( index == 0 )
                index_arr = 1:tvm.n_UBM;
            else
                index_arr = index;
            end
            n_data = size(data,1);
            for m=1:numel(index_arr)
                index = index_arr(m);
                % produce statistics
                stats = cell(n_data,1);
                for t=1:n_data
                    stats{t} = computeBwStatsFromBlocks(data{t},tvm.UBM{index});
                end
                % the version of train_tv_space4 wants a matrix, not a cell
                % array!
                stats = cat(2,stats{:});
                % make sure we don't build a depth that exceeds the
                % number of input elements
                if( tvm.Depth >= n_data )
                    tvm.Depth = n_data - tvm.Mixtures(index) - 1;
                end
                TVM = buildTVM(stats,tvm.UBM{index},tvm.Depth,...
                    tvm.Iterations,tvm.Verbosity);
                tvm.MatCond(index) = cond(TVM(:,:,end));
                tvm.TVM{index} = TVM;
                % produce learning rates
                rate = zeros(size(TVM,3)-1,1);
                rate_split = zeros(size(TVM,3)-1,tvm.Mixtures(index));
                for LR = 1:size(TVM,3)-1
                    rate(LR) = squeeze(mean(mean(abs(...
                        TVM(:,:,LR+1)./TVM(:,:,LR)))));
                    % split it into individual mixture bases
                    sliced = reshape(TVM,tvm.Depth,tvm.Features,...
                        tvm.Mixtures(index),[]);
                    for o=1:tvm.Mixtures(index)
                        rate_split(LR,o) = squeeze(mean(mean(abs(...
                            sliced(:,:,o,LR+1)./sliced(:,:,o,LR)))));
                    end
                end
                tvm.LearningRate{index,1} = rate;
                tvm.LearningRate{index,2} = rate_split;
            end
        end
        % evaluate enrollment and testing data
        function evaluate(tvm,enrollment_data,testing_data,channel,mixture)
            if( mixture == 0 )
                sub_total = size(enrollment_data,1)*size(testing_data,1);
                tvm.fpr = zeros(sub_total,tvm.n_Channels,tvm.n_UBM);
                tvm.fnr = zeros(sub_total,tvm.n_Channels,tvm.n_UBM);
                [tvm.crr(channel,:),tvm.eer(channel,:), ...
                    tvm.fnr(:,channel,:),tvm.fpr(:,channel,:)] = ...
                    iVectorEval(enrollment_data,testing_data,tvm.UBM);
            else
                sub_total = size(enrollment_data,1)*size(testing_data,1);
                tvm.fpr(:,:,mixture) = zeros(sub_total,tvm.n_Channels);
                tvm.fnr(:,:,mixture) = zeros(sub_total,tvm.n_Channels);
                [tvm.crr(channel,mixture),tvm.eer(channel,mixture), ...
                    tvm.fnr(:,channel,mixture),...
                    tvm.fpr(:,channel,mixture)] = ...
                    iVectorEval(enrollment_data,testing_data,...
                    tvm.UBM(mixture));
            end
        end
        % order the top two mixture matches
        function result = pair(tvm, mixture, depth)
            [vals,inds] = sort(tvm.RowSim{mixture,1}(:,:,end),'descend');
            if(depth >= tvm.Mixtures(mixture))
                depth = tvm.Mixtures(mixture)-1;
            end
            vals = [repmat(inds(1,:),1,depth-1); reshape(vals(2:depth,:)',1,[])];
            inds = [repmat(inds(1,:),1,depth-1); reshape(inds(2:depth,:)',1,[])];
            [~,best] = sort(vals(2,:));
            result = [inds(2,best); vals(:,best)];
            % how do the UBMs compare to these TVM pairings?
        end
        % compare relative positions of UBM, cosine distance, to TVM, MSE
        function result = shift(tvm, mixture)
            [v_tvm, i_tvm] = sort(tvm.RowSim{mixture,1}(:,:,end),'descend');
            [v_ubm, i_ubm] = sort(cosineDistance(tvm.UBM{mixture,1}.mu),...
                'descend');
            result = i_tvm==i_ubm;
            figure('numbertitle','off','name',...
                ['Mixture ' num2str(tvm.Mixtures(mixture)) ' Shift']);
            subplot(211);
            bar(sum(result(2:end,:),1)/(tvm.Mixtures(mixture)-1));
            grid on; ylim([0 1]);
            title('Unchanged order - entire list');
            subplot(212);
            s_mix = ceil(tvm.Mixtures(mixture)/2);
            bar(sum(result(2:s_mix,:),1)/(s_mix-1));
            grid on; ylim([0 1]);
            xlabel('UBM Mixture');
            title('Unchanged ordered - top half');
        end
        % use normalized rows to inform a decision about the similarity
        % between the sub-matrices
        function similarity(tvm,flag)
            result = cell(tvm.n_UBM,1);
            if( flag == 0 )
                parfor m=1:tvm.n_UBM
                    result{m} = TotalVariabilityMatrix.rowSim(...
                        tvm.TVM{m}(:,:,end),m);
                end
                tvm.RowSim = result;
            else
                parfor m=1:tvm.n_UBM
                    result{m} = (TotalVariabilityMatrix.rowSim(...
                        tvm.TVM_norm{m},m) + tvm.Depth)./(tvm.Depth*2);
                end
                tvm.RowSim_norm = result;
            end
            
        end
        % produce mappings of the similarity
        function [ubm_map, tvm_map, gradient] = map(tvm,mixture)
            ubm_dist = cosineDistance(tvm.UBM{mixture,1}.mu);
            gradient = (tvm.RowSim{mixture,1}(:,:,end)+20)/40 - ...
                (ubm_dist+1)/2;
            [vals_ubm, inds_ubm] = sort(ubm_dist,'descend');
            [vals_tvm, inds_tvm] = sort(tvm.RowSim{mixture,1}(:,:,end),...
                'descend');
            % cat together for output
            ubm_map = cat(3,inds_ubm,vals_ubm);
            tvm_map = cat(3,inds_tvm,vals_tvm);
        end
        % produce I-Vectors
        function [ldaIV,IV] = ivectorGen(tvm,data,index)
            % enable all training if index == 0
            if( index == 0 )
                index_arr = 1:tvm.n_UBM;
                IV = cell(numel(index_arr),1);
                ldaIV = IV;
            else
                index_arr = index;
                IV = cell(1,1);
                ldaIV = IV;
            end
            n_data = size(data,1);
            n_minor = size(data,2);
            for m = 1:numel(index_arr)
                ubm = tvm.UBM{m,1};
                index = index_arr(m);
                % produce statistics
                stats = cell(n_data,n_minor);
                for r=1:n_minor
                    for t=1:n_data
                        stats{t,r} = computeBwStatsFromBlocks(data{t,r},...
                            tvm.UBM{index});
                    end
                end
                % the version of train_tv_space4 wants a matrix, not a cell
                % array!
                major_iter = size(stats,1);
                minor_iter = size(stats,2);
                IVs = zeros(tvm.Depth,major_iter,minor_iter);
                for sub=1:major_iter
                    for ses=1:minor_iter
                        IVs(:,sub,ses) = extract_ivector2(stats{sub,ses},...
                            ubm,tvm.TVM{m,1}(:,:,end));
                    end
                end
                % LDA!
                lda_dim = min(tvm.Depth, major_iter-1);
                sID = (1:1:major_iter);
                subject_id = repmat(sID,1,minor_iter);
                IV{m} = reshape(IVs, tvm.Depth, major_iter*minor_iter);
                [tvm.V{m},tvm.D{m}] = lda(IV{m}, subject_id(:));
                ldaIV{m} = tvm.V{m}(:, 1:lda_dim)' * IV{m};
            end
        end
        % produce remapped ubm_map and sim_map if given key
        function [ubm,sim] = remap(tvm,key,mixture)
            if( isempty(tvm.RowSim) )
                similarity(tvm);
            end
            [ubm,sim] = map(tvm,mixture);
            ubm = ubm(:,:,1);
            sim = sim(:,:,1);
            [vals,inds] = sort(key,'descend');
            for r=1:tvm.Mixtures(mixture)
                ubm(ubm==r) = inds(1,r);
                sim(sim==r) = inds(1,r);
            end
        end
        % modify TVM through LDA matrix and compare against unmodified TVM
        function TVM_norm = ldaTVM(tvm)
            TVM_norm = cell(tvm.n_UBM,1);
            parfor m=1:tvm.n_UBM
                lda_dim = size(tvm.V{m},1);
                iters = size(tvm.TVM{m},2);
                temp = zeros(lda_dim,iters);
                for s=1:iters
                    temp(:,s) = tvm.V{m} * tvm.TVM{m}(:,s,end);
                end
                TVM_norm{m,1} = temp;
            end
            tvm.TVM_norm = TVM_norm;
        end
        % save the TVM!
        function s = saveobj(tvm,save_folder,save_name)
            full_save = [save_folder filesep save_name '.mat'];
            % build structure
            s.Depth = tvm.Depth;
            s.Iterations = tvm.Iterations;
            s.MatCond = tvm.MatCond;
            s.TVM = tvm.TVM;
            s.Error = tvm.Error;
            s.Sets = tvm.Sets;
            s.SetsThreshold = tvm.SetsThreshold;
            s.fullMSE = tvm.fullMSE;
            s.LearningRate = tvm.LearningRate;
            s.RowSim = tvm.RowSim;
            s.V = tvm.V;
            s.D = tvm.D;
            s.TVM_norm = tvm.TVM_norm;
            s.RowSim_norm = tvm.RowSim_norm;
            s.crr = tvm.crr;
            s.eer = tvm.crr;
            s.fnr = tvm.fnr;
            s.fpr = tvm.fpr;
            save(full_save,'s');
        end
    end
    methods (Static)
        % load a TVM from an object
        function obj = loadobj(s)
            if isstruct(s)
                newObj = TotalVariabilityMatrix(s.ubm,s.Depth,...
                    s.Iterations,s.n_Channels);
                newObj.MatCond = s.MatCond;
                newObj.TVM = s.TVM;
                newObj.Error = s.Error;
                newObj.Sets = s.Sets;
                newObj.SetsThreshold = s.SetsThreshold;
                newObj.fullMSE = s.fullMSE;
                newObj.LearningRate = s.LearningRate;
                newObj.RowSim = s.RowSim;
                newObj.V = s.V;
                newObj.D = s.D;
                newObj.TVM_norm = s.TVM_norm;
                newObj.RowSim_norm = s.RowSim_norm;
                newObj.crr = s.crr;
                newObj.eer = s.eer;
                newObj.fnr = s.fnr;
                newObj.fpr = s.fpr;
                obj = newObj;
            else
                obj = s;
            end
        end
        % perform row similarity calculation
        function result = rowSim(t_mats,m_count)
            mixtures = 2^m_count;
            t_dep = size(t_mats,1);
            iters = size(t_mats,2);
            t_fea = iters/mixtures;
            t_norm = normc(reshape(t_mats',t_fea,[]));
            t_norm = permute(reshape(t_norm,t_fea,t_dep,[]),[1 3 2]);
            inter_dist = zeros(mixtures,mixtures,t_dep);
            for r=1:t_dep
                inter_dist(:,:,r) = cosineDistance(t_norm(:,:,r));
            end
            result = sum(inter_dist,3);
        end
        % helper?
        function out = errorCalc(input,type)
            [n_dim,n_feat,n_mix] = size(input);
            if( type == 0 )
                out = zeros(n_mix,n_mix);
                % use this to produce final weights from TVM slice
                test_iv = ones(1,n_dim);
                for o=1:n_mix
                    slice = test_iv * input(:,:,o);
                    for i=1:n_mix
                        out(o,i) = [test_iv*input(:,:,i)] / slice;
                    end
                end
            end
        end
    end
end