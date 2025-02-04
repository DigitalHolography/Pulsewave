function ShowOutputs(PW_paths,output_dir)
% This function show multiple outputs from the foldermanagement drawerlist

N = length(PW_paths);

for ind =1:N
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
    
    if exist(fullfile(pw_path,last_PW_folder_name,'png','mask',[main_foldername,'_Segmentation.png']))
        segmentation_paths{ind} = fullfile(pw_path,last_PW_folder_name,'png','mask',[main_foldername,'_Segmentation.png']);
    end
    if exist(fullfile(pw_path,last_PW_folder_name,'png','volumeRate',[main_foldername,'_volumeRate_allrad_Artery_time.png']))
        bvr_paths{ind} = fullfile(pw_path,last_PW_folder_name,'png','volumeRate',[main_foldername,'_volumeRate_allrad_Artery_time.png']);
    end
    if exist(fullfile(pw_path,last_PW_folder_name,'png','pulseAnalysis',[main_foldername,'_1_Arteries_fRMS.png']))
        Arteries_fRMS_paths{ind} = fullfile(pw_path,last_PW_folder_name,'png','pulseAnalysis',[main_foldername,'_1_Arteries_fRMS.png']);
    end
    if exist(fullfile(pw_path,last_PW_folder_name,'png','pulseAnalysis',[main_foldername,'_ARI_velocity.png']))
        ARI_velocity_paths{ind} = fullfile(pw_path,last_PW_folder_name,'png','pulseAnalysis',[main_foldername,'_ARI_velocity.png']);
    end
    if exist(fullfile(pw_path,last_PW_folder_name,'png','bloodFlowVelocity',[main_foldername,'_histogramVelocityArteries.png']))
        histo_art_velocity_paths{ind} = fullfile(pw_path,last_PW_folder_name,'png','bloodFlowVelocity',[main_foldername,'_histogramVelocityArteries.png']);
    end
    if exist(fullfile(pw_path,last_PW_folder_name,'png','pulseAnalysis',[main_foldername,'_2_Arteries_velocity.png']))
        art_velocity_paths{ind} = fullfile(pw_path,last_PW_folder_name,'png','pulseAnalysis',[main_foldername,'_2_Arteries_velocity.png']);
    end
    if exist(fullfile(pw_path,last_PW_folder_name,'png','volumeRate',[main_foldername,'_ARI_BVR.png']))
        ARI_BVR_paths{ind} = fullfile(pw_path,last_PW_folder_name,'png','volumeRate',[main_foldername,'_ARI_BVR.png']);
    end
    if exist(fullfile(pw_path,last_PW_folder_name,'png','volumeRate',[main_foldername,'_strokeAndTotalVolume.png']))
        Stroke_total_volume{ind} = fullfile(pw_path,last_PW_folder_name,'png','volumeRate',[main_foldername,'_strokeAndTotalVolume.png']);
    end
    
end

[l, L] = bestMontageLayout(n);

figure(320)
montage(segmentation_paths, Size = [l L]);
exportgraphics(gca,fullfile(output_dir,'segmentations.png'));
figure(321)
montage(Arteries_fRMS_paths, Size = [l L]);
exportgraphics(gca,fullfile(output_dir,'ArteriesfRMS.png'));
figure(322)
montage(ARI_velocity_paths, Size = [l L]);
exportgraphics(gca,fullfile(output_dir,'ARIvelocity.png'));
figure(323)
montage(bvr_paths, Size = [l L]);
exportgraphics(gca,fullfile(output_dir,'bloodVolumeRate.png'));
figure(324)
montage(histo_art_velocity_paths, Size = [l L]);
exportgraphics(gca,fullfile(output_dir,'histogramVelocityArteries.png'));
figure(325)
montage(art_velocity_paths, Size = [l L]);
exportgraphics(gca,fullfile(output_dir,'ArteriesVelocity.png'));
figure(326)
montage(ARI_BVR_paths, Size = [l L]);
exportgraphics(gca,fullfile(output_dir,'ARIvelocityBVR.png'));
figure(327)
montage(Stroke_total_volume, Size = [l L]);
exportgraphics(gca,fullfile(output_dir,'_strokeAndTotalVolume.png'));

end