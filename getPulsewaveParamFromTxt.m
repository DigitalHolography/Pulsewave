
function [flag,parameter] =  getPulsewaveParamFromTxt(path,str)
myPath = split(path,'\');
myPath{size(myPath,1)} = '..';
myWritePath = myPath;
myPath{size(myPath,1)+1} = '..';

myWritePath{size(myWritePath,1)+1} = 'txt';
myWritePath = join(myWritePath,'\');

myPath{size(myPath,1)+1} = myPath{size(myPath,1)-3};
myPath = join(myPath,'\');


filename_txt = '_PulsewaveParams.txt';
[~,filename,~] = fileparts(myPath);
filename_txt = strcat(filename,filename_txt);
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
            disp(join(['Parameter' str 'found' ]))
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
    