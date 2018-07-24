% overlord function for NEW experiment running with channel data

function experimentTrio(param_file, workers)

% enhanced debug feature
% setSchedulerMessageHandler(@disp);
% clear existing cluster
delete(gcp('nocreate'))
% setup parallel environment
cluster = parcluster('local');
cluster.NumWorkers = workers;
% this needs to be ON for cluster operation and OFF for local debugging
% cluster.JobStorageLocation = getenv('HOME');
parpool(cluster,workers, 'IdleTimeout', Inf);

% process parameters
[list_args, subjects, sessions, iterations, element_listing, mixtures, ...
    ubm_iters, ds_factor, save_location] = readParameters(param_file);

% build data from file lists
[data_source,~,~] = listAgg(list_args{:});

% save parameters
sub_ses_m = cell(iterations,1);
sub_ses_g = sub_ses_m;
sub_ses_i = sub_ses_m;
sub_m = sub_ses_m;
sub_g = sub_ses_m;
sub_i = sub_ses_m;
ses_m = sub_ses_m;
ses_g = sub_ses_m;
ses_i = sub_ses_m;
SCORES_m = sub_ses_m;
SCORES_g = sub_ses_m;
SCORES_i = sub_ses_m;

for i=1:iterations
    fprintf('Iteration %d of %d\n',i,iterations);
    
    % produce train/test split (train is enroll  for I-Vector)
    [train, test] = dataCurate(data_source,element_listing);
    
    % build ubms
    ubms = allUBMs(train(:), mixtures, ubm_iters, ds_factor,0,0);
    
    % perform experiments
    fprintf('Mahal: %d\n',i);
    [sub_ses_m{i}, sub_m{i}, ses_m{i}, SCORES_m{i}] = ...
        mahalCcrSlim(train,test,subjects,sessions);
    fprintf('GMM: %d\n',i);
    [sub_ses_g{i}, sub_g{i}, ses_g{i}, SCORES_g{i}] = ...
        ubmGmmEvalRocca2(ubms,train,test,subjects,sessions);
    fprintf('I-Vector: %d\n',i);
    [sub_ses_i{i}, sub_i{i}, ses_i{i}, SCORES_i{i}] = ...
        iVectorEvalRocca2(ubms,train,test,subjects,sessions);
end
% break cells into arrays
SUB_SES_m = cell2mat(sub_ses_m);
SUB_SES_g = cell2mat(sub_ses_g);
SUB_SES_i = cell2mat(sub_ses_i);
SUB_m = cell2mat(sub_m);
SUB_g = cell2mat(sub_g);
SUB_i = cell2mat(sub_i);
SES_m = cell2mat(ses_m);
SES_g = cell2mat(ses_g);
SES_i = cell2mat(ses_i);
SCORES_m = cell2mat(SCORES_m);
SCORES_g = cell2mat(SCORES_g);
SCORES_i = cell2mat(SCORES_i);
% save experiments
mahal_folder = [save_location filesep 'MAHAL' filesep];
gmm_folder = [save_location filesep 'GMM' filesep];
ivector_folder = [save_location filesep 'IVECTOR' filesep];
% make sure the directories exist
directoryCheck(mahal_folder);
directoryCheck(gmm_folder);
directoryCheck(ivector_folder);
% save the files
saveFiles(mahal_folder,SUB_SES_m,SUB_m,SES_m,SCORES_m);
saveFiles(gmm_folder,SUB_SES_g, SUB_g,SES_g,SCORES_g);
saveFiles(ivector_folder,SUB_SES_i,SUB_i,SES_i,SCORES_i);

end


function saveFiles(folder,sub_ses,sub,ses,scores)
iVectorBinary([folder 'sub_ses.bin'],sub_ses);
iVectorBinary([folder 'sub.bin'],sub);
iVectorBinary([folder 'ses.bin'],ses);
iVectorBinary([folder 'scores.bin'],scores);
end