function [] =  checkPulsewaveParamsFromTxt(path)
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
    disp("Parameter file already exists")

else
    disp("Parameter file does not exist, writing in process")
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
    

