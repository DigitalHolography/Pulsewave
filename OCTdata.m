classdef OCTdata
    
    %     OCTdata is in PulseWave just to evitate the following Warnings :
    
    %     Warning: Cannot load an object of class 'OCTdata':
    %     Its class cannot be found.
    %     Warning: Cannot load an object of class 'OCTdata':
    %     Its class cannot be found.
    
    %     When the variable 'cache' is called from the function
    %     'getTimelineParamsFromCache.m' in 'pulseAnalysis.m'
    
    properties (Access = public)
        stack
        projection_xz
        projection_xy
        range_y
        range_z
    end
    
    methods
        function obj = OCTdata()
            obj.range_y = 1:10;
            obj.range_z = 1:10;
            obj.stack = zeros(512,512,512);
        end
    end
end