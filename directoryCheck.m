function directoryCheck(directory_path)
display(['Checking directory: ' directory_path]);
if( exist(directory_path,'dir') == 0)
    % directy does not exist, makeone
    disp('Directory does not exist. One will be created.');
    mkdir(directory_path);
end
end