function [mask_output] = removeDisconnected(mask_input, mask_input2, circle, name, ToolBox)

mask_output = mask_input & bwareafilt(mask_input2 | circle, 1, 8);
saveImage(mask_output + circle * 0.5, ToolBox, sprintf('%s_clear.png', name), isStep = true)
end
