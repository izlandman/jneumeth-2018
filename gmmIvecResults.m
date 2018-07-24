% handle GMM & I-Vector result production when running La Rocca
% experiments. probably a bit more difficult to clean up 2 feature sets and
% 2 classifiers, but it should make the previous function easier to read!

% should build data to work on and then pass that data to another function
% that handles evaluation

function [CCR, EER, scores, FPR, FNR] = gmmIvecResults(varargin)

if( nargin == 9 )
    channel_count = varargin{1};
    epochs = varargin{2};
    full_index = varargin{3};
    subject_count = varargin{4};
    mixtures = varargin{5};
    eval_flag = varargin{6};
    iterations = varargin{7};
    ds_factor = varargin{8};
    covar_flag = varargin{9};
else
    features = 1:26;
    channel_count = varargin{1};
    epochs = varargin{2};
    full_index = varargin{3};
    subject_count = varargin{4};
    mixtures = varargin{5};
    eval_flag = varargin{6};
    iterations = varargin{7};
    ds_factor = varargin{8};
    covar_flag = varargin{9};
    block_style = varargin{10};
    session_count = varargin{11};
    block_size = varargin{12};
    window_overlap = varargin{13};
end

n_mixtures = numel(mixtures);
x_subjects = subject_count * subject_count;
CCR = zeros(channel_count, n_mixtures, epochs);
EER = zeros(channel_count, n_mixtures, epochs);
FPR = zeros(x_subjects, n_mixtures, channel_count, epochs);
FNR = FPR;
scores = zeros(subject_count, subject_count, n_mixtures, channel_count, ...
    epochs);


if ( channel_count == 56 || channel_count == 1540 )
    % La Rocca Features then GMMs or I-Vectors
    fprintf('La Rocca Features: %d\n', channel_count);
    az = tic;
    
    for cv_step = 1:epochs
        % fprintf('>>> Iteration %d of %d.\n', cv_step, epochs);
        % build La Rocca features
        [training_data, testing_data] = laRoccaDataBuilderUBM( ...
            full_index, subject_count, epochs, cv_step);
        ubms = allUBMs(training_data(:), mixtures, iterations, ...
            ds_factor,0,0);
        CCR_temp = zeros(channel_count,n_mixtures);
        EER_temp = zeros(channel_count,n_mixtures);
        scores_temp = zeros(subject_count,subject_count,n_mixtures,...
            channel_count);
        FPR_temp = zeros(x_subjects,n_mixtures,channel_count);
        FNR_temp = FPR_temp;
        %%% remember ubm/gmm/ivector wants cell array of subject x channels
        %%% with each entry being the features x samples
        if( eval_flag == 0 )
            % produce GMM based results
            parfor c=1:channel_count
                [CCR_temp(c,:), EER_temp(c,:), ...
                    scores_temp(:,:,:,c), FPR_temp(:,:,c), ...
                    FNR_temp(:,:,c)] = ubmGmmEvalRocca(...
                    training_data, testing_data(:,c), ubms, covar_flag, c);
            end
        elseif( eval_flag == 1 )
            % produce IVECTOR based results
            parfor c=1:channel_count
                [CCR_temp(c,:), EER_temp(c,:), ...
                    scores_temp(:,:,:,c), FPR_temp(:,:,c), ...
                    FNR_temp(:,:,c)] = iVectorEvalRocca(...
                    training_data, testing_data(:,c), ubms, covar_flag, c);
            end
        else
            fprintf('Bad eval_flag!\n');
        end
        CCR(:,:,cv_step) = CCR_temp;
        EER(:,:,cv_step) = EER_temp;
        scores(:,:,:,:,cv_step) = scores_temp;
        FPR(:,:,:,cv_step) = FPR_temp;
        FNR(:,:,:,cv_step) = FNR_temp;
        fprintf('>>> Iteration %d of %d completed in %f seconds.\n', ...
            cv_step, epochs, toc(az));
    end
    
    
elseif( channel_count == 22 )
    % Cepstrum Features with GMM or I-Vectors
    
    fprintf('Cepstrum Features - ');
    az = tic;
    if( block_style < 1 )
        % randomize!
    elseif( block_style == 1 )
        for cv_step=1:epochs
            % build the features. ubms generated as needed based upon
            % channel
            [testing_data, training_data] = dataUbmPrep(subject_count, ...
                session_count, full_index, features, block_size, ...
                window_overlap, epochs, cv_step);
            CCR_temp = zeros(channel_count,n_mixtures);
            EER_temp = zeros(channel_count,n_mixtures);
            scores_temp = zeros(subject_count,subject_count,n_mixtures,...
                channel_count);
            FPR_temp = zeros(x_subjects,n_mixtures,channel_count);
            FNR_temp = FPR_temp;
            % UBMS!
            ubms = allUBMs(training_data(:), mixtures, iterations, ...
                ds_factor,0,0);
            if( eval_flag == 0 )
                % IVECTORS
                if( cv_step == 1 )
                    fprintf('GMMs\n');
                end
                parfor c=1:channel_count
                    [CCR_temp(c,:), EER_temp(c,:), ...
                        scores_temp(:,:,:,c), FPR_temp(:,:,c), ...
                        FNR_temp(:,:,c)] = ubmGmmEvalRocca(...
                        training_data,testing_data(:,c),ubms,covar_flag,c);
                end
            elseif( eval_flag == 1 )
                % GMMS
                if( cv_step == 1 )
                    fprintf('IVECTORs\n');
                end
                parfor c=1:channel_count
                    [CCR_temp(c,:), EER_temp(c,:), ...
                        scores_temp(:,:,:,c), FPR_temp(:,:,c), ...
                        FNR_temp(:,:,c)] = iVectorEvalRocca(...
                        training_data,testing_data(:,c),ubms,covar_flag,c);
                end
            else
                fprintf('Bad eval_flag for cepstrum features!\n');
            end
            CCR(:,:,cv_step) = CCR_temp;
            EER(:,:,cv_step) = EER_temp;
            scores(:,:,:,:,cv_step) = scores_temp;
            FPR(:,:,:,cv_step) = FPR_temp;
            FNR(:,:,:,cv_step) = FNR_temp;
            fprintf('>>> Iteration %d of %d completed in %f seconds\n', ...
                cv_step, epochs, toc(az));
        end
    else
        % error
        fprintf('Error. Incorrect block_style.\n');
    end
else
    fprintf('Incorrect number of channels!\n');
    return;
end

end