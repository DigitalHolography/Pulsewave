classdef Parameters_json < handle
% Class for storing parameters from an input json file

properties
    path
    name
    json
    exportVideos
    veins_analysis
    timePeriodMin
    px_size
end

methods

    function obj = Parameters_json(dir_path, filename)
        obj.path = dir_path;
        obj.name = filename;

        obj = obj.GetParameters();
    end

    function obj = GetParameters(obj)
        % Constructor method
        %[~, filename, ~] = fileparts(obj.path);
        filename_json = obj.name;
        dir_path_json = fullfile(obj.path, 'eyeflow', 'json');
        jsonPath = fullfile(dir_path_json, filename_json);

        if isfile(jsonPath)
            jsonData = fileread(jsonPath);
            parsedData = jsondecode(jsonData);

            obj.json = parsedData;
            % Recherche de chaque paramètre dans le fichier json

            obj.veins_analysis = parsedData.VeinsAnalysis;
            obj.exportVideos = parsedData.ExportVideos;
            obj.timePeriodMin = parsedData.MinimumGifPeriod;
            obj.px_size = parsedData.BloodVolumeRateAnalysis.PixelSize / (2 ^ parsedData.Preprocess.InterpolationFactor);

        else
            error('The json file could not be found.');
        end

    end

    function WriteParametersToJson(obj, outputPath)

        % Convert the structure into a JSON string
        jsonString = jsonencode(obj.json, "PrettyPrint", true);

        % Write the JSON string to the specified output path
        fileID = fopen(outputPath, 'w');

        if fileID == -1
            error('Could not open the file for writing: %s', outputPath);
        end

        fwrite(fileID, jsonString, 'char');
        fclose(fileID);

        disp('New JSON file created successfully.');
    end

end

end
