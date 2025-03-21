function topvel2csv(t, top_vel, top_vel_std, name)

ToolBox = getGlobalToolBox;

%csv output of the widths
T = table();
numR = length(top_vel); % number of radii

if not(isempty(t)) % if time is a variable
    T.time = t';
end

for rIdx = 1:numR
    numSection = size(top_vel{rIdx}, 1);

    for sectionIdx = 1:numSection
        T.(sprintf('Max_Vel_R%d_S%d_%s', rIdx, sectionIdx, name)) = squeeze(squeeze(top_vel{rIdx}(sectionIdx, :)))';
        T.(sprintf('STD_Max_Vel_R%d_S%d_%s', rIdx, sectionIdx, name)) = squeeze(squeeze(top_vel_std{rIdx}(sectionIdx, :)))';
    end

end

writetable(T, fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'MaxVelocityTable', '_', name, '.csv')));

end
