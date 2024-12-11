function [outputArg1,outputArg2] = ShowOutputs(PW_paths,output_dir)
% This function show multiple outputs from the foldermanagement drawerlist
figure(435)
title("Segmentation")

figure(321)
title("Blood Volume Rate")

for ind =1:length(PW_paths)
    split_path = strsplit(PW_paths{ind}, '\');
    main_foldername = split_path{end};
    PW_folder_name = strcat(main_foldername, '_PW');
    pw_path = fullfile(PW_paths{ind},'pulsewave');
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
    
    if exist(fullfile(pw_path,last_PW_folder_name,'png','mask',[main_foldername,'_arteryVeinSegmentation.png']))
        segmentation_paths{ind} = fullfile(pw_path,last_PW_folder_name,'png','mask',[main_foldername,'_arteryVeinSegmentation.png']);
    end
    if exist(fullfile(pw_path,last_PW_folder_name,'png','volumeRate',[main_foldername,'_volumeRateallradxtime.png']))
        bvr_paths{ind} = fullfile(pw_path,last_PW_folder_name,'png','volumeRate',[main_foldername,'_volumeRateallradxtime.png']);
    end

end

figure(435)
montage(segmentation_paths);
exportgraphics(gca,fullfile(output_dir,'segmentations.png'));
figure(321)
montage(bvr_paths);
exportgraphics(gca,fullfile(output_dir,'bloodVolumeRate.png'));

end