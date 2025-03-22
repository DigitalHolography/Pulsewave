function obj = VideoRegistering(obj)
tic
% Registers the video using intensity based registration
params = Parameters_json(obj.directory, obj.param_name);

if ~params.json.Preprocess.Register.Flag
    return % do nothing if not required
end

video = obj.M0_ff_video;
numX = size(video, 1);
numY = size(video, 2);

disc_ratio = 0.35; % parametrize this coef if needed

disc = diskMask(numX, numY, disc_ratio);
video_reg = video .* disc - disc .* sum(video .* disc, [1, 2]) / nnz(disc); % minus the mean in the disc of each frame
video_reg = reshape(video_reg, size(video, 1), size(video, 2), 1, size(video, 3)); % insert a dimension to match reegistration functions

video_reg = video_reg ./ (max(abs(video_reg), [], [1, 2])); % rescaling each frame but keeps mean at zero

refStart = params.json.Preprocess.Register.RefStart;
refEnd = params.json.Preprocess.Register.RefEnd;

image_ref = mean(video_reg(:, :, refStart:refEnd), 3); % ref image is from 10 to 20 for example
[~, shifts, scores] = register_video_from_reference(video_reg, image_ref);

figure(16), plot(scores / mean(scores), 'k'), ylim([0 max(0.1, 1.2 * max(scores) / mean(scores))]), title('Registration Correlation Score (1 is good) (u.a.)');
saveas(16, fullfile(obj.directory, 'eyeflow', sprintf("%s_%s", obj.filenames, 'RegistrationCorrelationScore.png')));
figure(17), plot(shifts(1, :)), hold on, plot(shifts(2, :))
imwrite(17, fullfile(obj.directory, 'eyeflow', sprintf("%s_%s", obj.filenames, 'RegistrationShiftsXY.png')));

obj.M0_ff_video = register_video_from_shifts(video, shifts);

obj.M0_data_video = register_video_from_shifts(obj.M0_data_video, shifts);
obj.M1_data_video = register_video_from_shifts(obj.M1_data_video, shifts);
obj.M2_data_video = register_video_from_shifts(obj.M2_data_video, shifts);

toc
close 16 17;

end
