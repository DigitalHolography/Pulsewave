function data = clear_empty_fields(data)

if isstruct(data)

    keys = fieldnames(data);
    data = rmfield(data, keys(structfun(@isempty, data)));

    keys = fieldnames(data);

    for i=1:numel(keys)
        
        data.(keys{i}) = clear_empty_fields(data.(keys{i}));

    end    

end