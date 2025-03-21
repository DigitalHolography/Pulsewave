function ClearParams(paths)
% This function clears all the input parameters files

for ind = 1:length(paths)
    EF_path_json = fullfile(paths{ind}, 'eyeflow', 'json');
    list_dir = dir(fullfile(EF_path_json, 'InputEyeFlowParams*.json'));

    matchingFiles = {list_dir.name};

    for i = 1:length(matchingFiles)
        delete(fullfile(EF_path_json, matchingFiles{i}))
    end

    % if isfile(fullfile(EF_path,'mask','forceMaskArtery.png'))
    %     movefile(fullfile(EF_path,'mask','forceMaskArtery.png'),fullfile(EF_path,'mask','oldForceMaskArtery.png'));
    % end
end
