paths = readlines("D:\HoloDopplerFolders\folders.txt");

%% ensure set default parameters and no forced mask

for ind = 1:length(paths)
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

    copyfile('DefaultPulsewaveParams.json',fullfile(pw_path,'json','InputPulsewaveParams.json'));

    if isfile(fullfile(pw_path,'mask','forceMaskArtery.png'))
        movefile(fullfile(pw_path,'mask','forceMaskArtery.png'),fullfile(pw_path,'mask','oldForceMaskArtery.png'));
    end
end


%% launch

for ind = 1:length(paths)
    path = paths(ind);
    if isfolder(path)
        path = strcat(path, '\');
    end
    OcClass = OneCycleClass(path);
    
    OcClass = OcClass.registerVideo();
    OcClass = OcClass.cropAllVideo();
    OcClass = OcClass.MomentNormalize();
    OcClass = OcClass.VideoResize();
    OcClass = OcClass.Interpolate();

    OcClass.flag_SH_analysis = 0;
    OcClass.flag_PulseWave_analysis = 0;
    OcClass.flag_velocity_analysis = 0;
    OcClass.flag_ExtendedPulseWave_analysis = 0;
    OcClass.flag_bloodVolumeRate_analysis = 0;
    OcClass.flag_bloodVelocityProfile_analysis = 0;

    OcClass.onePulse(255);
end

%% Show

Show_multiple_outputs;

