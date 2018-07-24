% run over three feature sets using three classifiers for ONE EXPERIMENT
% lists folder contains the listings of the three feature sets
% results folder is tiered with:
%   PSD - MAHAL, GMM, IVECTOR
%   COH - MAHAL, GMM, IVECTOR
%   CEP - MAHAL, GMM, IVECTOR

% when running on a folder, if data already exists skip that experiment!

function oneLaRocca(exp_folder,list_folder,trial,block_size,covar_flag,...
    workers)

% mixtures = [2 4 8 16 32 64 128 256 512 1024 2048 4096];
% data appears to small for larger mixture sizes, esp on 1 minute
% recordings. instead only go up to 128?
mixtures = [2 4 8 16 32 64 128 256 512];
ds_factor = 1;
iterations = 10;
% block_size = 100;
window_overlap = 0;
block_style = 1;

% find files of file_lists
file_types = {'CEP','PSD','COH'};
n_types = numel(file_types);
for i=1:n_types
    file_list_folder = [list_folder filesep file_types{i} ];
    file_list = findMatchingFiles(file_list_folder,trial);
    fprintf('Using file list: %s\n', file_list{1});
    % now call all the experiments!
    if( strcmp(file_types{i},'PSD')==1 || strcmp(file_types{i},'COH')==1)
        % MAHAL!
        save_folder = [exp_folder filesep '_RESULT' filesep ...
            file_types{i} filesep '_MAHAL'];
        directoryCheck(save_folder);
        if( expSkip(save_folder,'.bin') == 0 )
            laRoccaControl(file_list{1}, 1, save_folder, covar_flag, ...
                workers);
        end
        % GMM!
        save_folder = [exp_folder filesep '_RESULT' filesep ...
            file_types{i} filesep '_GMM'];
        directoryCheck(save_folder);
        if( expSkip(save_folder,'.bin') == 0 )
            laRoccaControl(file_list{1}, 1, save_folder, covar_flag, ...
                workers, mixtures, iterations, ds_factor, 0);
        end
        % IVECTOR!
        save_folder = [exp_folder filesep '_RESULT' filesep ...
            file_types{i} filesep '_IVECTOR'];
        directoryCheck(save_folder);
        if( expSkip(save_folder,'.bin') == 0 )
            laRoccaControl(file_list{1}, 1, save_folder, covar_flag, ...
                workers, mixtures, iterations, ds_factor, 1);
        end
    else
        % MAHAL!
        save_folder = [exp_folder filesep '_RESULT' filesep ...
            file_types{i} filesep '_MAHAL'];
        directoryCheck(save_folder);
        if( expSkip(save_folder,'.bin') == 0 )
            laRoccaControl(file_list{1}, 1, save_folder, covar_flag, ...
                workers);
        end
        % GMM!
        save_folder = [exp_folder filesep '_RESULT' filesep ...
            file_types{i} filesep '_GMM'];
        directoryCheck(save_folder);
        if( expSkip(save_folder,'.bin') == 0 )
            laRoccaControl(file_list{1}, 1, save_folder,covar_flag, ...
                workers,mixtures,iterations,ds_factor,0,block_size,...
                window_overlap,block_style);
        end
        % IVECTOR!
        save_folder = [exp_folder filesep '_RESULT' filesep ...
            file_types{i} filesep '_IVECTOR'];
        directoryCheck(save_folder);
        if( expSkip(save_folder,'.bin') == 0 )
            laRoccaControl(file_list{1},1,save_folder,covar_flag,...
                workers,mixtures,iterations,ds_factor,1,block_size,...
                window_overlap,block_style);
        end
    end
end

end

function result = expSkip(folder,file_type)
result = numel(findMatchingFiles(folder,file_type));
end