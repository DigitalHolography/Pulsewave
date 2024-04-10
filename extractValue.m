function value = extractValue(fileContent, tag)

startIndex = strfind(fileContent, tag);

if isempty(startIndex)
    value = [];
else

    % Finding the value
    startIndex = startIndex + length(tag) + 1;
    endIndex = strfind(fileContent(startIndex:end), newline);

    if isempty(endIndex)
        % If the tag is the last element in the file
        endIndex = length(fileContent);
    else
        endIndex = startIndex + endIndex(1) - 2; % Exclude endline character
    end

    valueStr = fileContent(startIndex:endIndex);
    value = str2double(valueStr);

    if isnan(value)
        value=[];
    end

end
end