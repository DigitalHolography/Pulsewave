function obj = readRaw(obj)

% if ~exist(strcat(obj.filenames, '_moment0.raw')) | ~exist(strcat(obj.filenames, '_moment1.raw')) | ~exist(strcat(obj.filenames, '_moment2.raw'))
%     error(' No raw moment files found. Please check folder path. Filenames should end with (_moment0.raw, _moment1.raw, _moment2.raw)) .')
% end
dir_path_avi = fullfile(obj.directory, 'avi');
NameRefAviFile = strcat(obj.filenames, '_M0');
RefAviFilePath = fullfile(dir_path_avi, NameRefAviFile);
ext = '.avi';
fprintf('- reading : %s\n', RefAviFilePath);

V = VideoReader(fullfile(dir_path_avi, [NameRefAviFile, ext]));
M0_disp_video = zeros(V.Height, V.Width, V.NumFrames);

for n = 1:V.NumFrames
    M0_disp_video(:, :, n) = rgb2gray(read(V, n));
end

obj.M0_ff_video = M0_disp_video;

clear V M0_disp_video

refvideosize = size(obj.M0_ff_video);
dir_path_raw = fullfile(obj.directory, 'raw');
ext = '.raw';

try
    % Import Moment 0

    NameRawFile = strcat(obj.filenames, '_moment0');
    fprintf('- reading : %s\n', fullfile(dir_path_raw, [NameRawFile, ext]));

    fileID = fopen(fullfile(dir_path_raw, [NameRawFile, ext]));
    M0_data_video = fread(fileID, 'float32');
    fclose(fileID);
    obj.M0_data_video = reshape(M0_data_video, refvideosize);
    clear M0_data_video

    % Import Moment 1

    NameRawFile = strcat(obj.filenames, '_moment1');
    fprintf('- reading : %s\n', fullfile(dir_path_raw, [NameRawFile, ext]));

    fileID = fopen(fullfile(dir_path_raw, [NameRawFile, ext]));
    M1_data_video = fread(fileID, 'float32');
    fclose(fileID);
    obj.M1_data_video = reshape(M1_data_video, refvideosize);
    clear M1_data_video

    % Import Moment 2

    NameRawFile = strcat(obj.filenames, '_moment2');
    fprintf('- reading : %s\n', fullfile(dir_path_raw, [NameRawFile, ext]));

    fileID = fopen(fullfile(dir_path_raw, [NameRawFile, ext]));
    M2_data_video = fread(fileID, 'float32');
    fclose(fileID);
    obj.M2_data_video = reshape(M2_data_video, refvideosize);
    clear M2_data_video

catch ME
    disp(['ID: ' ME.identifier])
    %
    rethrow(ME)

end

end
