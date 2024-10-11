function pulseVelocity(M,maskArtery,ToolBox,path)
%
[maskLongArtery,L,adjMatrix] = getLongestArteryBranch(maskArtery,M,ToolBox,path);
close all
%PWV = pulseWaveVelocity(M,maskLongArtery,ToolBox,path);
end