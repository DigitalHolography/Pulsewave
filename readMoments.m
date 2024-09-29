function [videoM0,videoM1,videoM2] = readMoments(filename)
% read moments stored in a .holo file format output by holovibes

header_mmap = memmapfile(filename, 'Format',...
            {'uint8', 4, 'magic_number';...
             'uint16', 1, 'version';...
             'uint16', 1, 'bit_depth';...
             'uint32', 1, 'width';...
             'uint32', 1, 'height';...
             'uint32', 1, 'num_frames';...
             'uint64', 1, 'total_size';...
             'uint8', 1,  'endianness';...
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
footer_skip = 64 + uint64(obj.frame_width * obj.frame_height) * uint64(obj.num_frames) * uint64(obj.bit_depth/8);
s=dir(filename);
footer_size = s.bytes - footer_skip;

% read old footer
 if  footer_skip >= s.bytes
        obj.footer.compute_settings.image_rendering.lambda = 8.5200e-07';
        obj.footer.info.pixel_pitch.x = 12;
        obj.footer.info.pixel_pitch.y = 12;
        obj.footer.compute_settings.image_rendering.propagation_distance = 0.4000;
else
    % open file and seek footer
    fd = fopen(filename, 'r');
    fseek(fd, footer_skip, 'bof');
    footer_unparsed = fread(fd, footer_size, '*char');
    footer_parsed = jsondecode(convertCharsToStrings(footer_unparsed));
    obj.footer = footer_parsed;
    fclose(fd);
end

fd = fopen(obj.filename, 'r');
fseek(fd, 64, 'bof'); % skip the header
frame_size = frame_width * obj.frame_height * 32; % as it is float32 images

if mod(num_frames,3) ~=0
    error('Bad moments file.')
end
videoM0 = zeros(frame_width, frame_height, floor(num_frames/3), 'single');
videoM1 = zeros(frame_width, frame_height, floor(num_frames/3), 'single');
videoM2 = zeros(frame_width, frame_height, floor(num_frames/3), 'single');

if endianness == 0
    endian= 'l';
elseif endianness == 1
    endian= 'b';
end

for i=1:num_frames

    videoM0(:,:,i) = reshape(fread(fd,frame_width * frame_height, 'single=>single', endian),[frame_width, frame_height]);
    videoM1(:,:,i) = reshape(fread(fd,frame_width * frame_height, 'single=>single', endian),[frame_width, frame_height]);
    videoM2(:,:,i) = reshape(fread(fd,frame_width * frame_height, 'single=>single', endian),[frame_width, frame_height]);
end
end