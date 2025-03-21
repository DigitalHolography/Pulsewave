function MEdisp(ME, path)
fprintf(2, "==========================================\nERROR\n==========================================\n")
fprintf(2, 'Error while loading : %s\n', path)
fprintf(2, "%s\n", ME.identifier)
fprintf(2, "%s\n", ME.message)

for stackIdx = 1:size(ME.stack, 1)
    fprintf(2, "%s : %s, line : %d\n", ME.stack(stackIdx).file, ME.stack(stackIdx).name, ME.stack(stackIdx).line);
end

if ME.identifier == "MATLAB:audiovideo:VideoReader:FileNotFound"

    fprintf(2, "No Raw File was found, please check 'save raw files' in HoloDoppler\n")

end

fprintf(2, "==========================================\n")
end
