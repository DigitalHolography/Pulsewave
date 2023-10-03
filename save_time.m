function [] = save_time (filename, function_name, time)

%get function name with coder.mfunctionname

fileID = fopen(filename,'a+') ;
fprintf(fileID,'%s : %f s \r\n',function_name,time);
fclose(fileID);

end