function makedir(dir_path)
    % make dir if it does not exist
    if ~exist(dir_path,'dir')
        mkdir(dir_path);
    end    
end