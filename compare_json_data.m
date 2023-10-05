    function [data_ref, old_data] = compare_json_data(data_ref, old_data)

% jsonData = fileread('InputPulsewaveParams.json');
% data_ref = jsondecode(jsonData);
% jsonData = fileread('Copy_of_InputPulsewaveParams.json');
% old_data = jsondecode(jsonData);

if isstruct(data_ref)

    keys = fieldnames(data_ref);

    for i=1:numel(keys)
        key = keys{i};
        try
            [data_ref.(key), old_data.(key)] = compare_json_data(data_ref.(key), old_data.(key));
        catch
            fprintf('Parameter added : %s (value set by default) \n', key)
        end
    end

else

    if ~isstruct(old_data)
        data_ref = old_data;
    end

end