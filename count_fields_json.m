function [n_fields] = count_fields_json(parsedData)

    n_fields = 0;

    try
        subfields = fieldnames(parsedData);

        for i = 1:numel(subfields)
            n_fields = n_fields + count_fields_json(parsedData.(subfields{i}));
        end

    catch
        n_fields = 1;
    end

end
