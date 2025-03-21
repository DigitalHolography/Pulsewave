function writeVideoOnDisc(data, filename, type)
% Writing function for videos to make use of parallel capacities with parfeval
%  filename - path to file and file name without extension
%  type     - nothing for avi, 'MPEG-4' for mp4

% Check input data
if isempty(data)
    error('Input data is empty.');
end

% Ensure data is in the correct format (4D: height × width × channels × frames)
if ndims(data) == 3 && size(data, 3) ~= 3 % Grayscale video
    data = reshape(data, size(data, 1), size(data, 2), 1, size(data, 3));
end

% Ensure data is of type uint8
if ~isa(data, 'uint8')
    data = uint8(rescale(data) * 255); % Convert to uint8
end

% Create VideoWriter object
if nargin > 2
    w = VideoWriter(filename, type);
else
    w = VideoWriter(filename);
end

% Open VideoWriter and write frames
open(w);

for jj = 1:size(data, 4)
    writeVideo(w, squeeze(data(:, :, :, jj)));
end

close(w);
end
