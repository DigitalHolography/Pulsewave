classdef Parameters_json < handle
    % Class for storing parameters from an input json file

    properties
        path
        name
        params
        registerVideoFlag
        refAvgStart
        refAvgEnd
        videoStartFrameIndex
        videoEndFrameIndex
        NormTempMode
        alphaConvolveNorm
        frameWidth
        frameHeight
        videoLength
        k
        removeOutliers
        radius_ratio
        radius_gap
        gauss_filt_size_for_barycenter
        oneCycleNinterp
        oneCycle_outNoiseThreshold
        oneCycle_dataReliabilityThreshold
        local_background_width
        masks_vesselness_sigma
        masks_vesselness_beta
        masks_crop_radius
        masks_vascular_threshold
        masks_vascular_classes
        masks_arterial_threshold
        masks_arterial_classes
        masks_venous_threshold
        masks_venous_classes
        masks_minSize
        masks_diaphragmRadius
        masks_imclose_radius
        masks_min_width
        masks_imdilateFinal
        CRACRV_Threshold
        cropSection_scaleFactorWidth
        cropSection_scaleFactorSize
        cropSection_maskThreshold
        cropSection_pixelSize
        elasWave_nDomFreq
        elasWave_pixelSize
        elasWave_gaussFiltPadding
        elasWave_butterFiltOrder
        elasWave_butterFiltBand
        exportVideos
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
        normFlag
        normPowerCalibrationSlope
        normPoweryIntercept
        normRefBloodFlow
        normActualCoefficient
        pupilRadius
        iris2retinaDist
        theta
        opticalIndex
        lambda
        phi
        pulseAnal_dataReliabilityFactor
        pulseAnal_peakHeightThreshold
        pulseAnal_outNoiseThreshold
        pulseAnal_blurScaleFactor
        pulseAnal_frameMinDiastole
        pulseAnal_framePeakSystole
        pulseAnal_exp
        trWavelength_f0
        trWavelength_r0
        veins_analysis
        velocitySmallRadiusRatio
        velocityBigRadiusRatio
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
        DiffFirstCalculationsFlag
        AllCirclesFlag
        nbCircles
        forcewidth
        forcebarycenter
        timePeriodMin
    end

    methods

        function obj = Parameters_json(dir_path,filename)
            obj.path = dir_path;
            obj.name = filename;

            obj = obj.GetParameters();
        end

        function obj = GetParameters(obj)
            % Constructor method
            %[~, filename, ~] = fileparts(obj.path);
            filename_json =  obj.name;
            dir_path_json = fullfile(obj.path, 'pulsewave', 'json');
            jsonPath = fullfile(dir_path_json, filename_json);

            if exist(jsonPath, 'file')
                jsonData = fileread(jsonPath);
                parsedData = jsondecode(jsonData);

                obj.params = parsedData;
                % Recherche de chaque paramètre
                obj.registerVideoFlag = parsedData.Video.Register;
                obj.refAvgStart = parsedData.Video.RefStart;
                obj.refAvgEnd = parsedData.Video.RefEnd;

                obj.videoStartFrameIndex = parsedData.Video.StartFrame;
                obj.videoEndFrameIndex = parsedData.Video.EndFrame;

                obj.NormTempMode = parsedData.NormalizeByTimeAverage;
                obj.alphaConvolveNorm = parsedData.MomentNormalizeConvolutionParameter;

                obj.frameWidth = parsedData.ResizeVideo.FrameWidth;
                obj.frameHeight = parsedData.ResizeVideo.FrameHeight;
                obj.videoLength = parsedData.ResizeVideo.VideoLength;

                obj.k = parsedData.ValueOfTheInterpolationParameter;

                obj.removeOutliers = parsedData.RemoveOutliersOption;

                obj.veins_analysis = parsedData.VeinsAnalysis;
                obj.exportVideos = parsedData.ExportVideos;


                obj.flatField_gwRatio = parsedData.FlatFieldCorrection.GWRatio;
                obj.flatField_border = parsedData.FlatFieldCorrection.Border;
                obj.flatField_borderDMap = parsedData.FlatFieldCorrection.BorderDMAP;
                obj.flatField_borderPulseAnal = parsedData.FlatFieldCorrection.BorderPulseAnalysis;

                obj.masks_vesselness_sigma = parsedData.CreationOfMasks.VesselnesParameterSigma;
                obj.masks_vesselness_beta = parsedData.CreationOfMasks.VesselnesParameterBeta;
                obj.gauss_filt_size_for_barycenter = parsedData.CreationOfMasks.GaussianFilterSizeForBarycenter;
                obj.masks_crop_radius = parsedData.CreationOfMasks.CropCoroidRadius;
                obj.masks_vascular_threshold = parsedData.CreationOfMasks.VascularCorrelationMapThreshold;
                obj.masks_vascular_classes = parsedData.CreationOfMasks.VascularClasses;
                obj.masks_arterial_threshold = parsedData.CreationOfMasks.ArterialCorrelationMapThreshold;
                obj.masks_arterial_classes = parsedData.CreationOfMasks.ArterialClasses;
                obj.masks_venous_threshold = parsedData.CreationOfMasks.VenousCorrelationMapThreshold;
                obj.masks_venous_classes = parsedData.CreationOfMasks.VenousClasses;
                obj.masks_minSize = parsedData.CreationOfMasks.MinimumSeedAreaSize;
                obj.masks_diaphragmRadius = parsedData.CreationOfMasks.DiaphragmRadius;
                obj.masks_imclose_radius = parsedData.CreationOfMasks.ImcloseRadius;
                obj.masks_min_width = parsedData.CreationOfMasks.MinimumVesselWidth;
                obj.masks_imdilateFinal = parsedData.CreationOfMasks.FinalDilation;
                obj.CRACRV_Threshold = parsedData.CreationOfMasks.CentralRetinaVeinArteryThreshold;

                obj.systoleThreshold = parsedData.SystoleDetection.SystoleThresholdRatioOfMaximum;

                obj.oneCycleNinterp = parsedData.CreationOfOneCycle.InterpolationPoints;
                obj.oneCycle_outNoiseThreshold = parsedData.CreationOfOneCycle.OutOfNoiseThreshold;
                obj.oneCycle_dataReliabilityThreshold = parsedData.CreationOfOneCycle.DataReliabilityIndexThreshold;

                obj.DiffFirstCalculationsFlag = parsedData.Velocity.DiffFirstCalculationsFlag;
                obj.local_background_width = parsedData.Velocity.LocalBackgroundWidth;

                obj.velocitySmallRadiusRatio = parsedData.SizeOfField.SmallRadiusRatio;
                obj.velocityBigRadiusRatio = parsedData.SizeOfField.BigRadiusRatio;

                obj.nbCircles = parsedData.BloodVolumeRate.NumberOfCircles;
                obj.radius_ratio = parsedData.BloodVolumeRate.RadiusRatio;
                obj.radius_gap = parsedData.BloodVolumeRate.RadiusGap;
                obj.cropSection_scaleFactorWidth = parsedData.BloodVolumeRate.ScaleFactorWidth;
                obj.cropSection_scaleFactorSize = parsedData.BloodVolumeRate.ScaleFactorSize;
                obj.cropSection_maskThreshold = parsedData.BloodVolumeRate.MaskSliceThreshold;
                obj.cropSection_pixelSize = parsedData.BloodVolumeRate.PixelSize;
                obj.flowRate_sliceHalfThickness = parsedData.BloodVolumeRate.SliceHalfThickness;

                obj.pupilRadius = parsedData.PulseAnalysis.RadiusPupil;
                obj.iris2retinaDist = parsedData.PulseAnalysis.DistanceIrisRetina;
                obj.theta = parsedData.PulseAnalysis.Theta;
                obj.opticalIndex = parsedData.PulseAnalysis.OpticalIndex;
                obj.lambda = parsedData.PulseAnalysis.Lambda;
                obj.phi = parsedData.PulseAnalysis.Phi;
                obj.pulseAnal_dataReliabilityFactor = parsedData.PulseAnalysis.DataReliatibilityFactor; %?
                obj.pulseAnal_peakHeightThreshold = parsedData.PulseAnalysis.PeakHeightThreshold;
                obj.pulseAnal_outNoiseThreshold = parsedData.PulseAnalysis.OutOfNoiseThreshold;
                obj.pulseAnal_blurScaleFactor = parsedData.PulseAnalysis.BlurScaleFactor;
                obj.pulseAnal_frameMinDiastole = parsedData.PulseAnalysis.FramePercentageBeforeMinimumOfDiastole;
                obj.pulseAnal_framePeakSystole = parsedData.PulseAnalysis.FramePercentageAroundPeakSystole;
                obj.pulseAnal_exp = parsedData.PulseAnalysis.Exponentiel;

                obj.video2vessels_radiusRatio = parsedData.VesselsVideo.RadiusRatioFactor;
                obj.video2vessels_gaussFiltPadding = parsedData.VesselsVideo.GaussianFilterPadding;
                obj.video2vessels_gaussFiltFactor = parsedData.VesselsVideo.GaussianFilterFactor;

                obj.viscosity_Ninterp = parsedData.Viscosity.NumberOfInterpolation;
                obj.viscosity_interpParam = parsedData.Viscosity.InterpolationParameter;
                obj.viscosity_interpSizeFactor = parsedData.Viscosity.InterpolationSizeFactor;
                obj.viscosity_projimThreshold = parsedData.Viscosity.ImageProjectionThreshold;
                obj.viscosity_fctnFactor = parsedData.Viscosity.ViscosityFunctionFactor;
                obj.viscosity_fitMatrix = parsedData.Viscosity.FitFunctionParameterMatrix;
                obj.viscosity_listParamA = parsedData.Viscosity.ViscosityListParameterA;
                obj.viscosity_listParamB = parsedData.Viscosity.ViscosityListParameterB;

                obj.elasWave_nDomFreq = parsedData.ElasticWaves.NumberOfDominantFrequency;
                obj.elasWave_pixelSize = parsedData.ElasticWaves.PixelSize;
                obj.elasWave_gaussFiltPadding = parsedData.ElasticWaves.GaussianFilterPadding;
                obj.elasWave_butterFiltOrder = parsedData.ElasticWaves.ButterworthFilterOrder;
                obj.elasWave_butterFiltBand = parsedData.ElasticWaves.ButterworthFilterThicknessBand;

                obj.nbSides = parsedData.Other.NumberOfSides;
                obj.AllCirclesFlag = parsedData.Other.AllCircles;
                obj.forcewidth = parsedData.Other.ForceWidthInPixels;
                obj.forcebarycenter = parsedData.Other.ForceBarycenter;
                obj.timePeriodMin = parsedData.Other.MinimumGifPeriod;

            else
                error('The json file could not be found.');
            end

        end

        function WriteParametersToJson(obj, outputPath)

            % Convert the structure into a JSON string
            jsonString = jsonencode(obj.params, "PrettyPrint", true);

            % Write the JSON string to the specified output path
            fileID = fopen(outputPath, 'w');
            if fileID == -1
                error('Could not open the file for writing: %s', outputPath);
            end
            fwrite(fileID, jsonString, 'char');
            fclose(fileID);

            disp('New JSON file created successfully.');
        end

    end

end
