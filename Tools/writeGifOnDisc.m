function writeGifOnDisc(data, filename, timePeriodMin, numFramesFixed, opt)
% Writing function for gif to make use of parallel capacities
% with parfeval
%  filename - - path to file and file name with '.gif' extension

arguments
    data
    filename
    timePeriodMin = NaN
    numFramesFixed = NaN
    opt.ToolBox = []
end

if size(size(data)) == [1, 3]
    data = reshape(data, size(data, 1), size(data, 2), 1, size(data, 3));
end

numFrames = size(data, 4);

gifWriter = GifWriter(filename, numFrames, ...
    timePeriodMin, numFramesFixed, "ToolBox", opt.ToolBox);

for frameIdx = 1:numFrames
    gifWriter.write(data(:, :, :, frameIdx), frameIdx);
end

gifWriter.generate();
gifWriter.delete();
end
