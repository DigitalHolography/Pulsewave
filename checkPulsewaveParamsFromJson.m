function [] =  checkPulsewaveParamsFromJson(path)

jsonInput = fileread("InputPulsewaveParams.json");
init_data = jsondecode(jsonInput);

%n_fields = count_fields_json(init_data);

[~,filename,~] = fileparts(path);
filename_json = 'DefaultsPulsewaveParams.json';
dir_path_json = fullfile(path,'json');


%filename_json = strcat(filename,filename_json);
jsonFilePath = fullfile(dir_path_json,filename_json);
json_exists = exist(jsonFilePath);

if json_exists 
    disp("Parameter file already exists")

    jsonData = fileread(jsonFilePath);
    parsedData = jsondecode(jsonData);

    [correct_data, ~] = compare_json_data(init_data, parsedData);

    disp("New parameter file created, writing in process")

    delete(jsonFilePath)

    jsonData = jsonencode(correct_data, PrettyPrint=true);

    jsonFilePath = fullfile(dir_path_json,filename_json);

    fileID = fopen(jsonFilePath, 'w');
    fprintf(fileID, jsonData);
    fclose(fileID);




else
    disp("Parameter file does not exist, writing in process")

    jsonData = jsonencode(init_data, PrettyPrint=true);

    jsonFilePath = fullfile(dir_path_json,filename_json);

    fileID = fopen(jsonFilePath, 'w');
    fprintf(fileID, jsonData);
    fclose(fileID);

end
end


