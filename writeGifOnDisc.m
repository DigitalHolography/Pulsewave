function writeGifOnDisc(data,filename,timePeriod)
% Writing function for gif to make use of parallel capacities 
% with parfeval
%  filename - - path to file and file name with '.gif' extension
N_frame = size(data,4);
gifWriter = GifWriter(filename, timePeriod, 0.04, N_frame);

for frameIdx = 1:N_frame
    gifWriter.write(data(:, :, :, frameIdx), frameIdx);
end

gifWriter.generate();
gifWriter.delete();
end