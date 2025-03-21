function readMomentsFooter(directory)
% Reads the rendering parameters of a .holo processed file and save them
% inside a json file
filename = strcat(directory, '.holo');
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
fseek(fd, footer_skip, 'bof');
footer_unparsed = fread(fd, footer_size, '*char');
footer_parsed = jsondecode(convertCharsToStrings(footer_unparsed));
fclose(fd);
json_txt = jsonencode(footer_parsed, PrettyPrint = true);
fid = fopen(fullfile(directory, 'Holovibes_rendering_parameters.json'), 'w+');
fprintf(fid, json_txt);
fclose(fid);

end
