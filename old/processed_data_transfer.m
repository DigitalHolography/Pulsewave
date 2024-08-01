function [] = processed_data_transfer()
% Function to transfer processed data from the original folder with raw
% measurements to the folder with only processed data (.avi or .png for
% instance)

rawFolderPath = "Z:\LDH Data\Raw";
processedFolderPath = "Z:\LDH Data\Processed";

% rawFolderPath = string
% processedFolderPath = string

%% Turn off the Warning message : "Warning: Directory already exists." because it will happen a lot
warning('off','MATLAB:MKDIR:DirectoryExists');

%% List of Subfolders with raw measurements
d = dir(rawFolderPath);
% remove all files (isdir property is 0)
subfoldersName = d([d(:).isdir]);
% remove '.' and '..' 
subfoldersName = subfoldersName(~ismember({subfoldersName(:).name},{'.','..'}));
subfoldersName = {subfoldersName.name};


%% Search for each folder and copy only the results
for ii=1:length(subfoldersName)
    mkdir(processedFolderPath,string(subfoldersName(ii))); %create the folder of results in the processed folder if not existing
    tmp_processed_folder = fullfile(rawFolderPath,string(subfoldersName(ii))); 
    % search for each subfolder in the processed data
    tmp_dir = dir(tmp_processed_folder);
    tmp_subfolders = tmp_dir([tmp_dir(:).isdir]);
    tmp_subfolders = tmp_subfolders(~ismember({tmp_subfolders(:).name},{'.','..'}));
    tmp_subfolders = {tmp_subfolders.name};

    for kk=1:length(tmp_subfolders)
        mkdir(fullfile(processedFolderPath,string(subfoldersName(ii))),string(tmp_subfolders(kk))); %for each processing folder

        if contains(string(tmp_subfolders(kk)),'preview')
            copyfile(fullfile(tmp_processed_folder,string(tmp_subfolders(kk))),fullfile(processedFolderPath,string(subfoldersName(ii)),string(tmp_subfolders(kk))));
        else
            % access to the folders 'avi', 'mp4', 'png, ...
            sub_d = dir(fullfile(tmp_processed_folder,tmp_subfolders(kk)));
            subsubfoldersName = sub_d([sub_d(:).isdir]);
            subsubfoldersName = subsubfoldersName(~ismember({subsubfoldersName(:).name},{'.','..'}));
            subsubfoldersName = {subsubfoldersName.name};

            for jj=1:length(subsubfoldersName)
                if string(subsubfoldersName(jj))=="raw" % we don't need to duplicate the raw videos created in holowaves
                    1;
                elseif contains(string(subsubfoldersName(jj)),'.') % we will not transfer the files .cine or .mat or .mp4
                    1;
                else
                    mkdir(fullfile(processedFolderPath,string(subfoldersName(ii)),string(tmp_subfolders(kk))),string(subsubfoldersName(jj))); %for each type of file
                    copyfile(fullfile(tmp_processed_folder,string(tmp_subfolders(kk)),string(subsubfoldersName(jj))),fullfile(processedFolderPath,string(subfoldersName(ii)),string(tmp_subfolders(kk)),string(subsubfoldersName(jj))));

                end
            end
        end
    end
end

%% Turn on the Warning message for other programs
warning('on','MATLAB:MKDIR:DirectoryExists');

end