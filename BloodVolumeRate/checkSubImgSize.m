function [subImgSize] = checkSubImgSize(sub_images_r)

subImgSize = [];
i = 1;

while isempty(subImgSize)

    if ~isempty(sub_images_r{1, i})
        subImgSize = size(sub_images_r{1, i}{1, 1});
    end

    i = i + 1;
end

end
