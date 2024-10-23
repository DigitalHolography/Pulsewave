% Test and evaluation script

paths = readlines("C:\Users\Vladikavkaz\Documents\data_test_list.txt");

figure(435)
title("Segmentation")

figure(321)
title("Blood Volume Rate")

for ind =1:length(paths)
    split_path = strsplit(paths(ind), '\');
    main_foldername = split_path{end};
    PW_folder_name = strcat(main_foldername, '_PW');
    pw_path = fullfile(paths(ind),'pulsewave');
    list_dir = dir(pw_path);
    idx = 0;
    for i=1:length(list_dir)
        if contains(list_dir(i).name, PW_folder_name)
            match = regexp(list_dir(i).name, '\d+$', 'match');
            if ~isempty(match) && str2double(match{1}) >= idx
                idx = str2double(match{1}); %suffix
            end
        end
    end
    last_PW_folder_name = sprintf('%s_%d', PW_folder_name, idx);
    
    segmentation_paths{ind} = fullfile(pw_path,last_PW_folder_name,'png','mask',[main_foldername,'_arteryVeinSegmentation.png']);
    bvr_paths{ind} = fullfile(pw_path,last_PW_folder_name,'png','bloodVolumeRate',[main_foldername,'_bloodVolumeRateallradxtime.png']);
end

figure(435)
montage(segmentation_paths);
figure(321)
%montage(bvr_paths);