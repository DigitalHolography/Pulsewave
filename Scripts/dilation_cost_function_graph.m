paths = readlines("X:\250117\AUZ_L\aaa.txt");

%% ensure set default parameters and no forced mask

Ndilation = 9;

ratio = 2*512;

for ind = 1:length(paths)
    pw_path = fullfile(paths(ind),'pulsewave');
    jsonFiles = dir(fullfile(pw_path,'json', '*.json'));

    for i = 1:length(jsonFiles)
        filePath = fullfile(fullfile(pw_path,'json'), jsonFiles(i).name);
        delete(filePath);  
    end

    copyfile(fullfile('Parameters','DefaultPulsewaveParams.json'),fullfile(pw_path,'json','InputPulsewaveParams.json'));


    for j = 0:Ndilation
        
        fileID = fopen(fullfile(pw_path,'json','InputPulsewaveParams.json'), 'r');
        jsonData = fread(fileID, inf, 'uint8')';
        fclose(fileID);
        jsonData = char(jsonData);
        decodedData = jsondecode(jsonData);
        
        decodedData.CreationOfMasks.ForceVesselWidth = j;
        
        jsonStr = jsonencode(decodedData,"PrettyPrint",true);
        fileID = fopen(fullfile(pw_path,'json',sprintf('InputPulsewaveParams_%d.json',j)), 'w');
        fprintf(fileID, '%s', jsonStr);
        fclose(fileID);
    end
end


%% launch

for ind = 1:length(paths)
    path = paths(ind);
    if isfolder(path)
        path = strcat(path, '\');
    end
    OcClass = OneCycleClass(path);
    
    OcClass = OcClass.preprocessData();
    
    OcClass.flag_SH_analysis = 0;
    OcClass.flag_Segmentation = 1;
    OcClass.flag_PulseWave_analysis = 1;
    OcClass.flag_velocity_analysis = 0;
    OcClass.flag_ExtendedPulseWave_analysis = 0;
    OcClass.flag_bloodVolumeRate_analysis = 0;
    OcClass.flag_bloodVelocityProfile_analysis = 0;

    for i = 1:length(OcClass.PW_params_names)
        OcClass.PW_param_name = OcClass.PW_params_names{i};
        OcClass.onePulse();
    end
end

%% Show
figure(9958);
hold on;
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
    
    diffRMS = {};
    for j = 0:Ndilation
        
        fileID = fopen(fullfile(pw_path,sprintf('%s_%d', PW_folder_name, idx-(Ndilation-j)),'txt',sprintf("%s_advanced_outputs.txt",sprintf('%s', PW_folder_name))), 'r');
        if fileID>=0
            tline = fgetl(fileID);
            while ischar(tline)
                if contains(tline, 'Mean fRMS difference artery ')
                    tmp = regexp(tline, '\d+(\.\d+)?', 'match');
                    diffRMS{j+1} = str2double(tmp{1});
                end
                tline = fgetl(fileID);
            end
            fclose(fileID);
        else
            disp('error opening file')
        end
        
    end
    
    
    plot((2*(1:Ndilation)+1)/ratio ,cell2mat(diffRMS(2:end)), 'LineWidth', 2);
    %findpeaks(cell2mat(diffRMS(2:end)));
end
ylabel('Mean difference f_{RMS}');
xlabel('width on max dimension ratio');


