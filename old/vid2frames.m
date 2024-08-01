function I=vid2frames(videoObject,NumberFrames,m,n)
% vid2frames converts a video into a matrix 
% USE: I=vid2vrames(videoObject,NumberFrames,m,n)
%
% INPUT:
%   videoObject: The video as seen as an object by Matlab
%   NumberFrames: Number of frames in  the video (can choose integer or whole)
%   m,n: Size of the images in the video
%
% OUPUT:
%   I: final matrix which contains all the Frames (m x n x NumberFrames) 
%
% Authors: 05/2021 Salim BENNANI & Malek AZAIZ 

I=zeros(m,n,NumberFrames);

for i=1:NumberFrames
    Frame = read(videoObject, i);
    I(:,:,i)=rgb2gray(Frame);
end

end