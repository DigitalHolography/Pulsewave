function plot2csvForAllRadSection(t, vr_avg_cell, vr_std_cell, vr_avg_mat, vr_std_mat, type)

numR = size(vr_avg_cell, 2);


T = table();
T.time = t;
for rIdx = 1:numR
    numSection = size(vr_avg_cell{rIdx}, 1);
    for sectionIdx = 1:numSection
        T.(sprintf('AVGVolumeRate_R%d_S%d_%s', rIdx, sectionIdx, type)) = vr_avg_cell{rIdx}(sectionIdx, :);
        T.(sprintf('STDVolumeRate_R%d_S%d_%s', rIdx, sectionIdx, type)) = vr_std_cell{rIdx}(sectionIdx, :);
    end
end

for rIdx = 1:numR
    T.(sprintf('AVGVolumeRate_%d_Total_%s', rIdx, type)) = sum(vr_avg_mat(rIdx, :, :), 2);
    T.(sprintf('STDVolumeRate_%d_Total_%s', rIdx, type)) = sum(vr_std_mat(rIdx, :, :), 2);
end
end