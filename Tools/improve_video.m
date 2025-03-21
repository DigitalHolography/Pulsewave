function [video] = improve_video(video_input, contrast_enhancement_tol, temporal_filter_sigma, contrast_inversion, options)
% some post processing
% commonly done for rendering hologram videos
%
% video: a 4D array containing a raw video to save to a file
%
% contrast_enhancement_tol: a parameter to adjust video contrast ([] if not wanted)
% temporal_filter_sigma: magnitude of temporal gaussian filter ([] if no filter is wanted)
% contrast_inversion: if true, contrast of the video will be inverted
% export_raw: if true, the video is also exported as a raw file in the raw directory
% export_avg_img: if true, save the temporal average of the video as a png
%                 file in the 'png' directory

arguments
    video_input
    contrast_enhancement_tol
    temporal_filter_sigma
    contrast_inversion
    options.NoIntensity = false
    options.cornerNorm = false
end

video = video_input;

%% temporal filter
if ~isempty(temporal_filter_sigma)
    sigma = [0.0001 0.0001 temporal_filter_sigma];

    for c = 1:size(video, 3)
        video(:, :, c, :) = imgaussfilt3(squeeze(video(:, :, c, :)), sigma);
    end

end

%% fix intensity flashes

if ~options.NoIntensity
    video = video - mean(mean(video, 2), 1);
end

%% corner normalizations

if options.cornerNorm > 0
    Nx = size(video, 1);
    Ny = size(video, 2);
    [X, Y] = meshgrid(linspace(-Nx / 2, Nx / 2, Nx), linspace(-Ny / 2, Ny / 2, Ny));
    disc_ratio = options.cornerNorm;
    disc = X .^ 2 + Y .^ 2 < (disc_ratio * min(Nx, Ny) / 2) ^ 2;
    video = video ./ mean(video .* ~disc, [1, 2]);
end

%% contrast enhancement
if ~isempty(contrast_enhancement_tol)
    video = mat2gray(video);
end

%% contrast inversion
if contrast_inversion
    video = -1.0 * video;
end

end
