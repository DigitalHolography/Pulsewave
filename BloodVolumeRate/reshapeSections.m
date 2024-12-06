function [area_mat, width_std_mat, vr_avg_mat, vr_std_mat] = reshapeSections(numSections, area_r, width_std_r, vr_avg_r, vr_std_r)

numSectionsMax = max(numSections);
numCircles = size(numSections, 2);
numFrames = size(vr_avg_r{1,1},2);

area_mat = zeros(numCircles, numSectionsMax);
width_std_mat = zeros(numCircles, numSectionsMax);
vr_avg_mat = zeros(numCircles, numSectionsMax, numFrames);
vr_std_mat = zeros(numCircles, numSectionsMax, numFrames);

for circleIdx = 1:numCircles
    numSection = numSections(circleIdx);
    if numSection < numSectionsMax
        for sectionIdx = numSection:numSectionsMax
            area_mat(circleIdx, sectionIdx) = 0;
            width_std_mat(circleIdx, sectionIdx) = 0;
            vr_avg_mat(circleIdx, sectionIdx, :) = zeros(1, 1, numFrames);
            vr_std_mat(circleIdx, sectionIdx, :) = zeros(1, 1, numFrames);
        end
    end
    for sectionIdx = 1:numSection
        area_mat(circleIdx, sectionIdx) = area_r{circleIdx}(sectionIdx);
        width_std_mat(circleIdx, sectionIdx) = width_std_r{circleIdx}(sectionIdx);
        vr_avg_mat(circleIdx, sectionIdx, :) = vr_avg_r{circleIdx}(sectionIdx, :);
        vr_std_mat(circleIdx, sectionIdx, :) = vr_std_r{circleIdx}(sectionIdx, :);
    end
end
end