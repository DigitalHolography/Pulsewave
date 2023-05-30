function [flag,stride, Fs] = getTimelineParamsFromCache(one_cycle_dir)
% find .mat with cache
myFolders = split(one_cycle_dir,'\');
myFolders{size(myFolders,1)+1} = '..';
myFolders{size(myFolders,1)+1} = '..';
myFolders{size(myFolders,1)+1} = myFolders{size(myFolders,1)-4};

myReadPath = {myFolders{1},myFolders{2},'mat'};
myReadPath = join(myReadPath,'\');


%FIXME : aller dans le rep mat
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
    Fs = cache.Fs;
    flag = 1;
    disp('done.')
else
    flag = 0;
    stride = 0;
    Fs = 0;
end
end

