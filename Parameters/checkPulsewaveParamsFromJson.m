function [PW_paramsNames] = checkPulsewaveParamsFromJson(path)

    % This function checks if in the folder you can find a Pulsewave
    % parameter file. If not it creates a Default one. If you find obsolete files it also
    % fills them with the new parameters

    % Additionally this function returns the list of the names of all valid
    % PulseWaveParameters files found (they must be in the form
    % 'InputPulsewaveParams*.json')

    jsonInput = fileread(fullfile("Parameters","DefaultPulsewaveParams.json"));
    init_data = jsondecode(jsonInput);

    %n_fields = count_fields_json(init_data);

    % [~,filename,~] = fileparts(path);
    dir_path_json = fullfile(path, 'pulsewave', 'json');
    jsonFiles = dir(fullfile(dir_path_json, 'InputPulsewaveParams*.json'));
    filename_json = jsonFiles(1).name;

    %filename_json = strcat(filename,filename_json);
    jsonFilePath = fullfile(dir_path_json, filename_json);
    json_exists = exist(jsonFilePath);

    filename_txt = 'InputPulsewaveParams.txt';
    dir_path_txt = fullfile(path, 'txt');

    txtFilePath = fullfile(dir_path_txt, filename_txt);
    txt_exists = exist(txtFilePath);

    if json_exists
        disp("Parameter file already exists, updating in process")

        jsonData = fileread(jsonFilePath);
        parsedData = jsondecode(jsonData);

        [correct_data, ~] = compare_json_data(init_data, parsedData);

        delete(jsonFilePath)

        jsonData = jsonencode(correct_data, PrettyPrint = true);

        jsonFilePath = fullfile(dir_path_json, filename_json);

        fileID = fopen(jsonFilePath, 'w');
        fprintf(fileID, jsonData);
        fclose(fileID);

    elseif txt_exists
        disp("Parameter file does not exist, found an old txt file. Writing a new json parameter file in process")

        fileContent = fileread(txtFilePath);

        json_data = txt2json_param(fileContent);
        [correct_data, ~] = compare_json_data(init_data, json_data);

        jsonData = jsonencode(correct_data, PrettyPrint = true);

        if ~isfolder(dir_path_json)
            mkdir(dir_path_json);
            disp(['Directory ', dir_path_json, ' has been created.']);
        end

        jsonFilePath = fullfile(dir_path_json, filename_json);

        fileID = fopen(jsonFilePath, 'w');
        fprintf(fileID, jsonData);
        fclose(fileID);

    else
        disp("Parameter file does not exist, writing in process")

        jsonData = jsonencode(init_data, PrettyPrint = true);

        if ~isfolder(dir_path_json)
            mkdir(dir_path_json);
            disp(['Directory ', dir_path_json, ' has been created.']);
        end

        jsonFilePath = fullfile(dir_path_json, filename_json);

        fileID = fopen(jsonFilePath, 'w');
        fprintf(fileID, jsonData);
        fclose(fileID);

    end


    % At this point at least one fparameter file exists. We collect the
    % names of all and return it.
    jsonFiles = dir(fullfile(dir_path_json, 'InputPulsewaveParams*.json'));

    % Initialize a cell array to store the file names
    PW_paramsNames = cell(1, numel(jsonFiles));
    
    % Store each file name in the cell array
    for i = 1:numel(jsonFiles)
        PW_paramsNames{i} = jsonFiles(i).name;
    end

end
