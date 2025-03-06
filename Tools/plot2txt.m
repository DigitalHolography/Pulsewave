function [] = plot2txt(tabx, taby, filename)
TB = getGlobalToolBox;
tabx = reshape(tabx, [length(tabx) 1]);
taby = reshape(taby, [length(tabx) 1]);
tmp = [tabx, taby];
fileID = fopen(fullfile(TB.path_txt, strcat(TB.main_foldername, '_', filename, '.txt')), 'w');
fprintf(fileID, '%f %f \r\n', tmp');
fclose(fileID);
end
