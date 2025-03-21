function [videoM0, videoM1, videoM2] = readMoments(filename)
% read moments stored in a .holo file format output by holovibes

header_mmap = memmapfile(filename, 'Format', ...
    {'uint8', 4, 'magic_number'; ...
     'uint16', 1, 'version'; ...
     'uint16', 1, 'bit_depth'; ...
     'uint32', 1, 'width'; ...
     'uint32', 1, 'height'; ...
     'uint32', 1, 'num_frames'; ...
     'uint64', 1, 'total_size'; ...
     'uint8', 1, 'endianness'; ...
% padding - skip
 }, 'Repeat', 1);

version = header_mmap.Data.version;
num_frames = header_mmap.Data.num_frames;
frame_width = header_mmap.Data.width;
frame_height = header_mmap.Data.height;
data_size = header_mmap.Data.total_size;
bit_depth = header_mmap.Data.bit_depth;
endianness = header_mmap.Data.endianness;

if ~isequal(header_mmap.Data.magic_number', unicode2native('HOLO'))
    error('Bad holo file.');
end

%% parse footer
footer_skip = 64 + uint64(frame_width * frame_height) * uint64(num_frames) * uint64(bit_depth / 8);
s = dir(filename);
footer_size = s.bytes - footer_skip;

fd = fopen(filename, 'r');
fseek(fd, 64, 'bof'); % skip the header
frame_size = frame_width * frame_height * 32; % as it is float32 images

videoM0 = zeros(frame_width, frame_height, floor(num_frames / 3), 'single');
videoM1 = zeros(frame_width, frame_height, floor(num_frames / 3), 'single');
videoM2 = zeros(frame_width, frame_height, floor(num_frames / 3), 'single');

if endianness == 0
    endian = 'l';
elseif endianness == 1
    endian = 'b';
end

max_idx = 1;

if mod(num_frames, 3) == 0
    max_idx = floor(num_frames / 3);
else
    max_idx = floor(num_frames / 3 - 1);
end

for i = 1:max_idx
    videoM0(:, :, i) = fftshift(reshape(fread(fd, frame_width * frame_height, 'single=>single', endian), [frame_width, frame_height]));
    videoM1(:, :, i) = fftshift(reshape(fread(fd, frame_width * frame_height, 'single=>single', endian), [frame_width, frame_height]));
    videoM2(:, :, i) = fftshift(reshape(fread(fd, frame_width * frame_height, 'single=>single', endian), [frame_width, frame_height]));
end

% implay(rescale(videoM0));
% implay(rescale(videoM1));
% implay(rescale(videoM2));
end
