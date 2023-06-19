classdef Parameters
    % Class for storing parameters from an input text file
    
    properties
        path
        k
        radius_ratio
        resistivity_gaussFiltSigma
        resistivity_satAmp
        resistivity_minTolVal
        resistivity_maxTolVal
        resistivity_gamma
        oneCycle_Ninterp
        oneCycle_outNoiseThreshold
        oneCycle_dataReliabilityThreshold
        arteryMask_vesselness_sigma
        arteryMask_vesselness_beta
        arteryMask_vesselness_bin_threshold
        arteryMask_magicwand_nb_of_area_vessels
        arteryMask_magicwand_nb_of_area_artery
        arteryMask_ArteryCorrThreshold
        centralRetinal_arteryThreshold
        centralRetinal_veinThreshold
        centralRetinal_backgndThreshold
        vesselMask_stdBinTreshold
        cropSection_scaleFactorWidth
        cropSection_scaleFactorSize
        cropSection_maskThreshold
        cropSection_pixelSize
        elasWave_nDomFreq
        elasWave_pixelSize
        elasWave_gaussFiltPadding
        elasWave_butterFiltOrder
        elasWave_butterFiltBand
        systoleThreshold
        flatField_gwRatio
        flatField_border
        flatField_borderDMap
        flatField_borderPulseAnal
        flowRate_gaussFiltPadding
        flowRate_yelorRangeHSV
        flowRate_cydblueRangeHSV
        flowRate_minTolVal
        flowRate_maxTolVal
        flowRate_sliceHalfThickness
        pupilRadius
        iris2retinaDist
        theta
        opticalIndex
        lambda
        pulseAnal_dataReliabilityFactor
        pulseAnal_peakHeightThreshold
        pulseAnal_outNoiseThreshold
        pulseAnal_blurScaleFactor
        pulseAnal_frameMinDiastole
        pulseAnal_framePeakSystole
        pulseAnal_exp
        trWavelength_f0
        trWavelength_r0
        video2vessels_radiusRatio
        video2vessels_gaussFiltPadding
        video2vessels_gaussFiltFactor
        viscosity_Ninterp
        viscosity_interpParam
        viscosity_interpSizeFactor
        viscosity_projimThreshold
        viscosity_fctnFactor
        viscosity_fitMatrix
        viscosity_listParamA
        viscosity_listParamB
        nbSides
    end
    
    methods

       function obj = Parameters(dir_path)
            obj.path = dir_path;
            persistent PulsewaveParams
            if isempty(PulsewaveParams) 
                PulsewaveParams = obj.GetParameters();
            end
           obj = PulsewaveParams;
       end

        function obj = GetParameters(obj)
            % Constructor method     
            %[~, filename, ~] = fileparts(obj.path);
            filename_txt = 'InputPulsewaveParams.txt';
            dir_path_txt = fullfile(obj.path,'txt');
            txtPath = fullfile(dir_path_txt, filename_txt);
            
            if exist(txtPath, 'file')
                fileContent = fileread(txtPath);

                % Recherche de chaque paramètre
                obj.k = obj.extractValue(fileContent, 'Value of the interpolation parameter :');
                obj.radius_ratio = obj.extractValue(fileContent, 'Radius ratio :');
                obj.resistivity_gaussFiltSigma = obj.extractValue(fileContent, 'Gaussian filter sigma :');
                obj.resistivity_satAmp = obj.extractValue(fileContent, 'Amplitude saturation  ponderation :');
                obj.resistivity_minTolVal = obj.extractValue(fileContent, 'Tolerante value minimum :');
                obj.resistivity_maxTolVal = obj.extractValue(fileContent, 'Tolerante value maximum :');
                obj.resistivity_gamma = obj.extractValue(fileContent, 'Gamma :');
                obj.oneCycle_Ninterp = obj.extractValue(fileContent, 'Number of interpolation :');
                obj.oneCycle_outNoiseThreshold = obj.extractValue(fileContent, 'Out of noise treshold :');
                obj.oneCycle_dataReliabilityThreshold = obj.extractValue(fileContent, 'Data reliability index treshold :');
                obj.arteryMask_vesselness_sigma = obj.extractValue(fileContent, 'Vesselnes parameter sigma :');
                obj.arteryMask_vesselness_beta = obj.extractValue(fileContent, 'Vesselnes parameter beta :');
                obj.arteryMask_vesselness_bin_threshold = obj.extractValue(fileContent, 'Vesselness Binarization treshold :');
                obj.arteryMask_magicwand_nb_of_area_vessels = obj.extractValue(fileContent, 'MagicWand Number of segmented area detected for vessels :');
                obj.arteryMask_magicwand_nb_of_area_artery = obj.extractValue(fileContent, 'MagicWand Number of segmented area detected for artery :');
                obj.arteryMask_ArteryCorrThreshold = obj.extractValue(fileContent, 'Artery correlation threshold :');
                obj.centralRetinal_arteryThreshold = obj.extractValue(fileContent, 'Central retina artery treshold :');
                obj.centralRetinal_veinThreshold = obj.extractValue(fileContent, 'Central retina vein treshold :');
                obj.centralRetinal_backgndThreshold = obj.extractValue(fileContent, 'Central retina background treshold :');
                obj.vesselMask_stdBinTreshold = obj.extractValue(fileContent, 'Standard binarization treshold :');
                obj.cropSection_scaleFactorWidth = obj.extractValue(fileContent, 'Scale factor width :');
                obj.cropSection_scaleFactorSize = obj.extractValue(fileContent, 'Scale factor size :');
                obj.cropSection_maskThreshold = obj.extractValue(fileContent, 'Mask slice treshold :');
                obj.cropSection_pixelSize = obj.extractValue(fileContent, 'Pixel size ? :');
                obj.elasWave_nDomFreq = obj.extractValue(fileContent, 'Number of dominant frequency :');
                obj.elasWave_pixelSize = obj.extractValue(fileContent, 'Pixel size :');
                obj.elasWave_gaussFiltPadding = obj.extractValue(fileContent, 'Gaussian filter padding :');
                obj.elasWave_butterFiltOrder = obj.extractValue(fileContent, 'Butterworth filter order :');
                obj.elasWave_butterFiltBand = obj.extractValue(fileContent, 'Butterworth filter thickness band :');
                obj.systoleThreshold = obj.extractValue(fileContent, 'Systole threshold :');
                obj.flatField_gwRatio = obj.extractValue(fileContent, 'GW ratio :');
                obj.flatField_border = obj.extractValue(fileContent, 'Border :');
                obj.flatField_borderDMap = obj.extractValue(fileContent, 'Border dMAP :');
                obj.flatField_borderPulseAnal = obj.extractValue(fileContent, 'Border PulseAnalisys :');
                obj.flowRate_gaussFiltPadding = obj.extractValue(fileContent, 'Gaussian filter padding :');
                obj.flowRate_yelorRangeHSV = obj.extractValue(fileContent, 'Yellow-orange range HSV :');
                obj.flowRate_cydblueRangeHSV = obj.extractValue(fileContent, 'Cyan-dark blue range HSV :');
                obj.flowRate_minTolVal = obj.extractValue(fileContent, 'Tolerante value minimum :');
                obj.flowRate_maxTolVal = obj.extractValue(fileContent, 'Tolerante value maximum :');
                obj.flowRate_sliceHalfThickness = obj.extractValue(fileContent, 'Slice half thickness :');
                obj.pupilRadius = obj.extractValue(fileContent, 'Radius pupil :');
                obj.iris2retinaDist = obj.extractValue(fileContent, 'Distance iris-retina :');
                obj.theta = obj.extractValue(fileContent, 'Theta :');
                obj.opticalIndex = obj.extractValue(fileContent, 'Optical index :');
                obj.lambda = obj.extractValue(fileContent, 'Lambda :');
                obj.pulseAnal_dataReliabilityFactor = obj.extractValue(fileContent, 'Data reliatibility factor ? :');
                obj.pulseAnal_peakHeightThreshold = obj.extractValue(fileContent, 'Peak height treshold :');
                obj.pulseAnal_outNoiseThreshold = obj.extractValue(fileContent, 'Out of noise treshold :');
                obj.pulseAnal_blurScaleFactor = obj.extractValue(fileContent, 'Blur scale factor :');
                obj.pulseAnal_frameMinDiastole = obj.extractValue(fileContent, 'Frame percentage before minimum of diastole :');
                obj.pulseAnal_framePeakSystole = obj.extractValue(fileContent, 'Frame percentage around peak systole :');
                obj.pulseAnal_exp = obj.extractValue(fileContent, 'exponentiel 1 ? :');
                obj.trWavelength_f0 = obj.extractValue(fileContent, 'F0 :');
                obj.trWavelength_r0 = obj.extractValue(fileContent, 'R0 :');
                obj.video2vessels_radiusRatio = obj.extractValue(fileContent, 'Radius ratio factor :');
                obj.video2vessels_gaussFiltPadding = obj.extractValue(fileContent, 'Gaussian filter padding :');
                obj.video2vessels_gaussFiltFactor = obj.extractValue(fileContent, 'Gaussian filter factor :');
                obj.viscosity_Ninterp = obj.extractValue(fileContent, 'Number of interpolation :');
                obj.viscosity_interpParam = obj.extractValue(fileContent, 'Interpolation parameter :');
                obj.viscosity_interpSizeFactor = obj.extractValue(fileContent, 'Interpolation size factor :');
                obj.viscosity_projimThreshold = obj.extractValue(fileContent, 'Image projection treshold :');
                obj.viscosity_fctnFactor = obj.extractValue(fileContent, 'Viscosity function factor :');
                obj.viscosity_fitMatrix = obj.extractValue(fileContent, 'Fit function parameter (matrix) :');
                obj.viscosity_listParamA = obj.extractValue(fileContent, 'Viscosity list parameter a :');
                obj.viscosity_listParamB = obj.extractValue(fileContent, 'Viscosity list parameter b :');
                obj.nbSides = obj.extractValue(fileContent, 'Number of sides :');

            else
                error('The text file could not be found.');
            end
        end
        
        function value = extractValue(~, fileContent, tag)
            % Extracts the value associated with the given tag from the file content
            
            startIndex = strfind(fileContent, tag);
            
            if isempty(startIndex)
                error('Tag not found in the file content.');
            end
            
            % Finding the value 
            startIndex = startIndex + length(tag) + 1;
            endIndex = strfind(fileContent(startIndex:end), newline);
            
            if isempty(endIndex)
                % If the tag is the last element in the file
                endIndex = length(fileContent);
            else
                endIndex = startIndex + endIndex(1) - 2; % Exclude endline character
            end
            
            valueStr = fileContent(startIndex:endIndex);
            value = str2double(valueStr);
            
            if isnan(value)
                error(['Invalid value for tag: ' tag]);
            end
        end
    


end
end
