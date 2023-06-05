function [] =  checkPulsewaveParamsFromTxt(path)

myPath = split(path,'\');
myWritePath = {myPath{1},myPath{2},'txt'};
myWritePath  = join(myWritePath,'\');
filename_txt = '_PulsewaveParams.txt';
filename_txt = strcat(myPath{2},filename_txt);
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
    

