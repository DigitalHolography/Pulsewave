function [paramsNames] = checkEyeFlowParamsFromJson(path)

% This function checks if in the folder you can find a EyeFlow
% parameter file. If not it creates a Default one. If you find obsolete files it also
% fills them with the new parameters

% Additionally this function returns the list of the names of all valid
% EyeFlowParameters files found (they must be in the form
% 'InputEyeFlowParams*.json')

dir_path_json = fullfile(path, 'eyeflow', 'json');

% We first check if an old txt parameter exists and write a json parameter file instead
filename_txt = 'InputEyeFlowParams.txt';
filename_json = 'InputEyeFlowParams.json';
dir_path_txt = fullfile(path, 'txt');
txtFilePath = fullfile(dir_path_txt, filename_txt);
txt_exists = exist(txtFilePath);

jsonInput = fileread(fullfile("Parameters", "DefaultEyeFlowParams.json"));
init_data = jsondecode(jsonInput);

if txt_exists
    disp("Found an old txt file. Writing a new json parameter file instead")
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
end

% We now check all the existing json files named like 'InputEyeFlowParams*.json'
jsonFiles = dir(fullfile(dir_path_json, 'InputEyeFlowParams*.json'));

if ~isempty(jsonFiles)
    disp("Found parameter files : ")

    for i = 1:numel(jsonFiles)
        jsonFilePath = fullfile(dir_path_json, jsonFiles(i).name);
        disp(jsonFilePath)
        jsonData = fileread(jsonFilePath);
        parsedData = jsondecode(jsonData);

        [correct_data, ~] = compare_json_data(init_data, parsedData);

        delete(jsonFilePath)

        jsonData = jsonencode(correct_data, PrettyPrint = true);

        jsonFilePath = fullfile(dir_path_json, jsonFiles(i).name);

        fileID = fopen(jsonFilePath, 'w');
        fprintf(fileID, jsonData);
        fclose(fileID);

    end

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

% At this point at least one parameter file exists. We collect the
% names of all and return it.
jsonFiles = dir(fullfile(dir_path_json, 'InputEyeFlowParams*.json'));

% Initialize a cell array to store the file names
paramsNames = cell(1, numel(jsonFiles));

% Store each file name in the cell array
for i = 1:numel(jsonFiles)
    paramsNames{i} = jsonFiles(i).name;
end

end
