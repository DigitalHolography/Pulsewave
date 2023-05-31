function [flag,k,radius_ratio] =  getPulsewaveParamsFromTxt(path)

myPath = split(path,'\');
myWritePath = {myPath{1},myPath{2},'txt'};
myWritePath  = join(myWritePath,'\');
filename_txt = '_PulsewaveParams.txt';
filename_txt = strcat(myPath{2},filename_txt);
txt_exists = exist(fullfile(myWritePath{1},filename_txt));


if txt_exists
    fileID = fopen(fullfile(myWritePath{1}, filename_txt),'r');
    while ~feof(fileID)
        line = fgetl(fileID);
        if strcmp('Value of the interpolation parameter :', line)
            line = fgetl(fileID);
            k = str2double(line);
            flag = 1;
        end
        if strcmp('Radius ratio :', line)
            line = fgetl(fileID);
            radius_ratio = str2double(line);
            flag = 2;
        end
    end
    fclose(fileID) ;

else
    flag = 0;
    k = 1;
    radius_ratio = 0.18;
    mkdir(myWritePath{1});
    fileID = fopen(fullfile(myWritePath{1}, filename_txt),'w');
    fprintf(fileID,[...
    'Value of the interpolation parameter :\n%d\n' ...
    'Radius ratio :\n%d\n'], ...
    k, ...
    radius_ratio);
    fclose(fileID) ;
end
    

