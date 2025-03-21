function [data_ref, data_test] = compare_json_data(data_ref, data_test)

% jsonData = fileread('InputEyeFlowParams.json');
% data_ref = jsondecode(jsonData);
% jsonData = fileread('Copy_of_InputEyeFlowParams.json');
% data_test = jsondecode(jsonData);

if isstruct(data_ref)

    keys = fieldnames(data_ref);

    for i = 1:numel(keys)
        key = keys{i};

        try
            [data_ref.(key), data_test.(key)] = compare_json_data(data_ref.(key), data_test.(key));
        catch
            fprintf('Parameter added : %s (value set by default) \n', key)
        end

    end

else

    if ~isstruct(data_test)
        data_ref = data_test;
    end

end
