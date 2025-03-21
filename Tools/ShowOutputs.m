function ShowOutputs(paths, output_dir)
% This function show multiple outputs from the foldermanagement drawerlist

N = length(paths);

for path_idx = 1:N
    split_path = strsplit(paths{path_idx}, '\');
    main_foldername = split_path{end};
    folder_name = strcat(main_foldername, '_EF');
    ef_path = fullfile(paths{path_idx}, 'eyeflow');
    list_dir = dir(ef_path);
    idx = 0;

    for i = 1:length(list_dir)

        if contains(list_dir(i).name, folder_name)
            match = regexp(list_dir(i).name, '\d+$', 'match');

            if ~isempty(match) && str2double(match{1}) >= idx
                idx = str2double(match{1}); %suffix
            end

        end

    end

    last_folder_name = sprintf('%s_%d', folder_name, idx);

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'mask', [main_foldername, '_vesselMap.png']))
        segmentation_paths{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'mask', [main_foldername, '_vesselMap.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'volumeRate', [main_foldername, '_allrad_Artery_time.png']))
        bvr_paths{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'volumeRate', [main_foldername, '_allrad_Artery_time.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'pulseAnalysis', [main_foldername, '_1_Arteries_fRMS_graph.png']))
        Arteries_fRMS_paths{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'pulseAnalysis', [main_foldername, '_1_Arteries_fRMS_graph.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'pulseAnalysis', [main_foldername, '_ARI_velocity_graph.png']))
        ARI_velocity_paths{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'pulseAnalysis', [main_foldername, '_ARI_velocity_graph.png']);
    end

    if isfile(fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_histogramVelocityArteries.png']))
        histo_art_velocity_paths{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'bloodFlowVelocity', [main_foldername, '_histogramVelocityArteries.png']);
    end

    %     if isfile(fullfile(ef_path,last_folder_name,'png','pulseAnalysis',[main_foldername,'_2_Arteries_velocity.png']))
    %         art_velocity_paths{path_idx} = fullfile(ef_path,last_folder_name,'png','pulseAnalysis',[main_foldername,'_2_Arteries_velocity.png']);
    %     end
    %     if isfile(fullfile(ef_path,last_folder_name,'png','volumeRate',[main_foldername,'_ARI_BVR.png']))
    %         ARI_BVR_paths{path_idx} = fullfile(ef_path,last_folder_name,'png','volumeRate',[main_foldername,'_ARI_BVR.png']);
    %     end
    if isfile(fullfile(ef_path, last_folder_name, 'png', 'volumeRate', [main_foldername, '_strokeAndTotalVolume.png']))
        Stroke_total_volume{path_idx} = fullfile(ef_path, last_folder_name, 'png', 'volumeRate', [main_foldername, '_strokeAndTotalVolume.png']);
    end

end

[l, L] = bestMontageLayout(N);

figure(320)
montage(segmentation_paths, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'segmentations.png'));
figure(321)
montage(Arteries_fRMS_paths, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'ArteriesfRMS.png'));
figure(322)
montage(ARI_velocity_paths, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'ARIvelocity.png'));
figure(323)
montage(bvr_paths, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'bloodVolumeRate.png'));
figure(324)
montage(histo_art_velocity_paths, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, 'histogramVelocityArteries.png'));
% figure(325)
% montage(art_velocity_paths, Size = [l L]);
% exportgraphics(gca,fullfile(output_dir,'ArteriesVelocity.png'));
% figure(326)
% montage(ARI_BVR_paths, Size = [l L]);
% exportgraphics(gca,fullfile(output_dir,'ARIvelocityBVR.png'));
figure(327)
montage(Stroke_total_volume, Size = [l L]);
exportgraphics(gca, fullfile(output_dir, '_strokeAndTotalVolume.png'));

end
