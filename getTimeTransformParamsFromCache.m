function [flag,type, F1, F2, minPCA, maxPCA] = getTimeTransformParamsFromCache(one_cycle_dir)
% find .mat with cache
myFolders = split(one_cycle_dir,'\');
myFolders = myFolders(1:size(myFolders,1)-1);
myFolders{size(myFolders,1)-1} = '..';


myReadPath = myFolders;
myReadPath{size(myReadPath,1)} = 'mat';
myReadPath = join(myReadPath,'\');

myFolders{size(myFolders,1)+1} = myFolders{size(myFolders,1)-3};

myNewPath = join(myFolders,'\');
myNewPath = [myNewPath{1},'.mat'];
[~,filename_mat,~] = fileparts(myNewPath);
% filename_mat = filename(1:length(filename)-3) ;
cache_exists = exist(fullfile(myReadPath{1},strcat(filename_mat,'.mat')));
%
if cache_exists % .mat with cache from holowaves is present, timeline can be computed
    disp('reading cache parameters');    
    load(fullfile(myReadPath{1},strcat(filename_mat,'.mat')),'cache') ; 
    stride = cache.batch_stride ;
    type = cache.time_transform.type;
    F1 = cache.time_transform.f1;
    F2 = cache.time_transform.f2;
    minPCA = cache.time_transform.min_PCA;
    maxPCA = cache.time_transform.max_PCA;
    flag = 1;
    
    disp('done.')
elseif exist(myNewPath)
    disp('reading cache parameters'); 
    load(myNewPath,'cache') ; 
    stride = cache.batch_stride ;
    type = cache.time_transform.type;
    F1 = cache.time_transform.f1;
    F2 = cache.time_transform.f2;
    minPCA = cache.time_transform.min_PCA;
    maxPCA = cache.time_transform.max_PCA;
    flag = 1;

else
    flag = 0;
    type = 'None';
    F1 = 0;
    F2 = 0;
    minPCA = 0;
    maxPCA = 0;
end
end

