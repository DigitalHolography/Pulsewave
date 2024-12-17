function plot2txtForAllRadSection(t, vr_avg_cell, vr_std_cell, vr_avg_mat, vr_std_mat, type)

numR = size(vr_avg_cell, 2);

for rIdx = 1:numR
    numSection = size(vr_avg_cell{rIdx}, 1);
    for sectionIdx = 1:numSection
        plot2txt(t, vr_avg_cell{rIdx}(sectionIdx, :), sprintf('AVGVolumeRate_R%d_S%d_%s', rIdx, sectionIdx, type))
        plot2txt(t, vr_std_cell{rIdx}(sectionIdx, :), sprintf('STDVolumeRate_R%d_S%d_%s', rIdx, sectionIdx, type))
    end
    plot2txt(t, sum(vr_avg_mat(rIdx, :, :), 2), sprintf('AVGVolumeRate_%d_Total_%s', rIdx, type))
    plot2txt(t, sum(vr_std_mat(rIdx, :, :), 2), sprintf('STDVolumeRate_%d_Total_%s', rIdx, type))
end

end