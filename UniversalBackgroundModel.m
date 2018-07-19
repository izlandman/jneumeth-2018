% produce a class for Universal Background Models

classdef UniversalBackgroundModel < handle
    properties
        % basic properties of the data
        
        % number of features
        n_Features
        % number of mixtures
        n_Mixtures
        % which mixtures?
        Mixtures
        % modeling iterations
        Iterations
        % store ds_factor
        dsFactor
        % store number of subjects
        n_Subjects
        % store number of channels
        n_Channels
        % the actual results
        UBMs
        % classification rate
        crr
        % equal error rate
        eer
        % false negative rate
        fnr
        % false positive rate
        fpr
        % similarity between UBM.mus
        ubm_dist
    end
    
    methods
        % instantiate the object by providing training data and a list of
        % mixtures to model
        function obj = UniversalBackgroundModel(mixtures,iterations,...
                dsFactor,subjects,channels)
            obj.Iterations = iterations;
            obj.n_Mixtures = numel(mixtures);
            obj.Mixtures = mixtures;
            obj.dsFactor = dsFactor;
            obj.n_Subjects = subjects;
            obj.n_Channels = channels;
            % outputs~
            obj.crr = zeros(obj.n_Channels,obj.n_Mixtures);
            obj.eer = obj.crr;
            obj.fnr = zeros(obj.n_Subjects*obj.n_Subjects,...
                obj.n_Channels,obj.n_Mixtures);
            obj.fpr = obj.fnr;
        end
        
        % train the ubms!
        function train(ubm,training_data,verbose)
            ubm.UBMs = allUBMs(training_data,ubm.Mixtures,...
                ubm.Iterations,ubm.dsFactor,0,verbose);
            ubm.n_Features = size(ubm.UBMs{1}.mu);
        end
        % once built, evaluate the UBMs against enrollment and testing data
        function evaluate(ubm,enrollment_data,testing_data,channel,mixture)
            if( mixture == 0 )
                [ubm.crr(channel,:),ubm.eer(channel,:),...
                    ubm.fnr(:,channel,:),...
                    ubm.fpr(:,channel,:)] = ...
                    ubmGmmTest(enrollment_data,testing_data,ubm.UBMs);
            else
                [ubm.crr(channel,mixture),ubm.eer(channel,mixture),...
                    ubm.fnr(:,channel,mixture),...
                    ubm.fpr(:,channel,mixture)] = ...
                    ubmGmmTest(enrollment_data,testing_data,...
                    ubm.UBMs(mixture));
            end
        end
        % produce a mapping of a given mixture against a set of data
        function [label, labels, index, llks] = map(ubm, map_data, ...
                mixture, map_style, block_style, block_size)
            
            % build data to be mapped if map_data is string
            [map_data,subjects,queries] = ubm.readData(ubm,map_data,...
                map_style,block_style,block_size);
            
            % using chosen UBM, match blocks to models via numeric index
            model = ubm.UBMs{mixture,1};
            n_model = 2^mixture;
            
            % evaluate, somewhat /extra/
            [s_recip, s_mu, C, w, f_pi] = compute_ubm_components(model);
            llks = zeros(queries,n_model);
            
            for m=1:n_model
                gmm.w = model.w(m);
                gmm.mu = model.mu(:,m);
                gmm.sigma = model.sigma(:,m);
                for q=1:queries
                    llks(q,m) = mean(compute_llk(map_data{q},...
                        gmm,f_pi));
                end
            end
            
            % find best fit for given block
            [~, label] = sort(llks,2,'descend');
            labels = zeros(n_model,1);
            for l=1:n_model
                labels(l) = sum( label(:,1)==l );
            end
            % make it a percentage, match with .w values
            labels = labels./queries;
            index = reshape(label(:,1),subjects,[]);
        end
        
        % perform data verification, enrollment and test are same data,
        % using htk based data
        function verification(ubm,data,mixture,target_style,block_style,...
                block_size)
            % gather data
            if( ~iscell(data) )
                [data,subjects,queries] = ubm.readData(ubm,data,...
                    target_style,block_style,block_size);
            end
            % when mixutre is 0 build for all mixtures
            if( mixture == 0 )
                [ubm.crr,ubm.eer,ubm.gmmScores] = ...
                    ubmGmmEvalRocca(data,data,ubm.UBMs,0,0);
            else
                [ubm.crr(mixture),ubm.eer(mixture),...
                    ubm.gmmScores(:,:,mixture)] = ...
                    ubmGmmEvalRocca(data,data,ubm.UBMs(mixture),0,0);
            end            
        end
        
        % produce similarity measure between means of mixtures
        function similarity(ubm)
            for m=1:ubm.n_Mixtures
                ubm.ubm_dist{m} = cosineDistance(ubm.UBMs{m}.mu);
            end
        end
        
        % save the UBM as a mat file, large but what can one do!
        function s = saveobj(obj,save_folder,save_name)
            % build structure
            s.n_Features = obj.n_Features;
            s.Mixtures = obj.Mixtures;
            s.Iterations = obj.Iterations;
            s.dsFactor = obj.dsFactor;
            s.n_Subjects = obj.n_Subjects;
            s.n_Channels = obj.n_Channels;
            s.UBMs = obj.UBMs;
            s.crr = obj.crr;
            s.eer = obj.eer;
            s.fnr = obj.fnr;
            s.fpr = obj.fpr;
            % save structure
            full_save = [save_folder filesep save_name];
            save(full_save,'s');
        end
    end
    
    methods(Static)
        % process input data into blocks
        function [map_data,subjects,queries] = readData(ubm,map_data,...
                map_style,block_style,block_size)
            if( ischar(map_data) )
                feature_sel = 1:ubm.n_Features;
                map_points = blockPrep(map_data, map_style, ...
                    block_style, block_size, feature_sel);
                [queries,~] = size(map_points);
                % sort out number of subjects?
                subjects = strsplit(map_points{end,2},'_');
                subjects = str2double(subjects{2});
                map_data = map_points(:,1);
                clear map_points
            else
                subjects = size(map_data,1);
                queries = subjects;
            end
        end
        
        % rebuild the object from a saved structure
        function obj = loadobj(s)
            if isstruct(s)
               newObj = UniversalBackgroundModel(s.Mixtures,s.Iterations,...
                   s.dsFactor,s.n_Subjects,s.n_Channels);
               newObj.crr = s.crr;
               newObj.eer = s.eer;
               newObj.fnr = s.fnr;
               newObj.fpr = s.fpr;
               newObj.UBMs = s.UBMs;
               obj = newObj;
            else
                obj = s;
            end
        end
    end
end

% borrow a bunch of functions to produce LLKs
function [s_recip,s_mu,C,w,f_pi] = compute_ubm_components(GMM)
s_recip = 1./GMM.sigma;
s_mu = GMM.mu.*s_recip;
C = sum(GMM.mu.*s_mu) + sum(log(GMM.sigma));
w = GMM.w(:);
f_pi = log(2*pi);
end
function llk = compute_llk(data, GMM, f_pi)
% compute the posterior probability of mixtures for each frame
post = lgmmprob(data, GMM.mu, GMM.sigma, GMM.w(:), f_pi);
llk  = logsumexp(post, 1);
end
function logprob = lgmmprob(data, mu, sigma, w, f_pi)
% compute the log probability of observations given the GMM
ndim = size(data, 1);
s_recip = 1./sigma;
s_mu = mu.*s_recip;
% C = sum(mu.*mu./sigma) + sum(log(sigma));
C = sum(mu.*s_mu) + sum(log(sigma));
% D = (1./sigma)' * (data .* data) - 2 * (mu./sigma)' ...
% * data  + ndim * log(2 * pi);
D = s_recip'*(data.*data)-2*s_mu'*data+ndim*f_pi;
logprob = -0.5 * (bsxfun(@plus, C',  D));
logprob = bsxfun(@plus, logprob, log(w));
end
function y = logsumexp(x, dim)
% compute log(sum(exp(x),dim)) while avoiding numerical underflow
xmax = max(x, [], dim);
y    = xmax + log(sum(exp(bsxfun(@minus, x, xmax)), dim));
ind  = find(~isfinite(xmax));
if ~isempty(ind)
    y(ind) = xmax(ind);
end
end