function [] = plot2txt(tabx, taby, filename, ToolBox)

tmp = [tabx, taby];

fileID = fopen(fullfile(ToolBox.PW_path_txt, strcat(ToolBox.main_foldername,'_', filename, '.txt')),'w') ;
fprintf(fileID,'%f %f \r\n',tmp');
fclose(fileID);

end