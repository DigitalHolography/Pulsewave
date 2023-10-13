function [] =  checkPulsewaveParamsFromJson(path)

jsonInput = fileread("DefaultsPulsewaveParams.json");
init_data = jsondecode(jsonInput);

%n_fields = count_fields_json(init_data);

[~,filename,~] = fileparts(path);
filename_json = 'InputPulsewaveParams.json';
dir_path_json = fullfile(path,'json');


%filename_json = strcat(filename,filename_json);
jsonFilePath = fullfile(dir_path_json,filename_json);
json_exists = exist(jsonFilePath);

if json_exists 
    disp("Parameter file already exists, updating in process")

    jsonData = fileread(jsonFilePath);
    parsedData = jsondecode(jsonData);

    [correct_data, ~] = compare_json_data(init_data, parsedData);

    delete(jsonFilePath)

    jsonData = jsonencode(correct_data, PrettyPrint=true);

    jsonFilePath = fullfile(dir_path_json,filename_json);

    fileID = fopen(jsonFilePath, 'w');
    fprintf(fileID, jsonData);
    fclose(fileID);




else
    disp("Parameter file does not exist, writing in process")

    jsonData = jsonencode(init_data, PrettyPrint=true);
    
    if ~isfolder(dir_path_json)
        mkdir(dir_path_json);
        disp(['Directory ', dir_path_json, ' has been created.']);
    end
    
    jsonFilePath = fullfile(dir_path_json,filename_json);

    fileID = fopen(jsonFilePath, 'w');
    fprintf(fileID, jsonData);
    fclose(fileID);

end
end


