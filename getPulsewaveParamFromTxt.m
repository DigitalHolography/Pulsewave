function [flag,parameter] =  getPulsewaveParamFromTxt(path,str)

myPath = split(path,'\');
myWritePath = {myPath{1},myPath{2},'txt'};
myWritePath  = join(myWritePath,'\');
filename_txt = '_PulsewaveParams.txt';
filename_txt = strcat(myPath{2},filename_txt);
txt_exists = exist(fullfile(myWritePath{1},filename_txt));

if txt_exists
    fileID = fopen(fullfile(myWritePath{1}, filename_txt),'r');
    flag=0;
    while ~feof(fileID)
        line = fgetl(fileID);
        if strcmp(str, line)
            line = fgetl(fileID);
            parameter = str2double(line);
            flag = 1;
        end
    end
    if flag==0
        disp(join(['Parameter' str 'not found' ]))
    end
    fclose(fileID) ;

else
    disp("File not found")
    flag = 0;
    parameter = 1;
end
    