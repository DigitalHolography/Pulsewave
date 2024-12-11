function ClearParams(PW_paths)
% This function clears all the input parameters files 

for ind = 1:length(PW_paths)
    pw_path_json = fullfile(PW_paths{ind},'pulsewave','json');
    list_dir = dir(fullfile(pw_path_json, 'InputPulseWaveParams*.json'));

    matchingFiles = {list_dir.name};

    for i=1:length(matchingFiles)
        delete(fullfile(pw_path_json,matchingFiles{i}))
    end
    % if isfile(fullfile(pw_path,'mask','forceMaskArtery.png'))
    %     movefile(fullfile(pw_path,'mask','forceMaskArtery.png'),fullfile(pw_path,'mask','oldForceMaskArtery.png'));
    % end
end