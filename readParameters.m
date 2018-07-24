function [list_args, subjects, sessions, iterations, element_listing, ...
    mixtures, ubm_iters, ds_factor, save_location] = ...
    readParameters(param_file)

% open file
try
    fid = fopen(param_file);
catch
    error('fopen(paramfile); failed!\n');
end

param_count = 9;
iter = 1;
while iter <= param_count
    new_line = fgetl(fid);
    if( ~strcmp(new_line(1),'#') )
        switch iter
            case 1
                % number of subjects
                subjects = str2double(new_line);
                iter = iter + 1;
            case 2
                % process filelist
                sessions = str2double(new_line);
                files = cell(sessions,1);
                for f=1:sessions
                    files{f} = fgetl(fid);
                end
                iter = iter + 1;
            case 3
                % handle blocks
                nums = strsplit(new_line,' ');
                b_size = str2double(nums{1});
                b_overlap = str2double(nums{2});
                iter = iter + 1;
            case 4
                % handle iterations
                iterations = str2double(new_line);
                iter = iter + 1;
            case 5
                % handle elements withheld
                element_listing = str2double(new_line);
                iter = iter + 1;
            case 6
                % handle mixtures
                mixtures = str2double(strsplit(new_line,' '));
                iter = iter + 1;
            case 7
                % handle ubm_iters
                ubm_iters = str2double(new_line);
                iter = iter + 1;
            case 8
                % handle ds_factor
                ds_factor = str2double(new_line);
                iter = iter + 1;
            case 9
                % handle save_location
                save_location = new_line;
                iter = iter + 1;
            otherwise
                fprintf('Something broke in switch statement.\n');
        end
    end
end
if( b_size == 0 )
    list_args = {files{:}};
else
    list_args = {b_size,b_overlap,files{:}};
end
fclose(fid);
end