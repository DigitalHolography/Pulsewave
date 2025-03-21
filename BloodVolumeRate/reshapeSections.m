function [A_mat, Q_mat, Q_std_mat] = reshapeSections(numFrames, numSections, A_r, Q_r, Q_std_r)

numSectionsMax = max(numSections);
numCircles = size(numSections, 2);

A_mat = zeros(numCircles, numSectionsMax);
Q_mat = zeros(numCircles, numSectionsMax, numFrames);
Q_std_mat = zeros(numCircles, numSectionsMax, numFrames);

for cIdx = 1:numCircles
    numSection = numSections(cIdx);

    if (numSection < numSectionsMax) && (numSection ~= 0)

        for sectionIdx = numSection:numSectionsMax
            A_mat(cIdx, sectionIdx) = nan;
            Q_mat(cIdx, sectionIdx, :) = nan(1, 1, numFrames);
            Q_std_mat(cIdx, sectionIdx, :) = nan(1, 1, numFrames);
        end

    end

    if numSection ~= 0

        for sectionIdx = 1:numSection
            A_mat(cIdx, sectionIdx) = A_r{cIdx}(sectionIdx);
            Q_mat(cIdx, sectionIdx, :) = Q_r{cIdx}(sectionIdx, :);
            Q_std_mat(cIdx, sectionIdx, :) = Q_std_r{cIdx}(sectionIdx, :);
        end

    end

end

end
