classdef Parameters_json < handle
    % Class for storing parameters from an input json file

    properties
        path
        registerVideoFlag
        refAvgStart
        refAvgEnd
        videoStartFrameIndex
        videoEndFrameIndex
        frameWidth
        frameHeight
        videoLength
        k
        radius_ratio
        radius_gap
        gauss_filt_size_for_barycentre
        oneCycle_outNoiseThreshold
        oneCycle_dataReliabilityThreshold
        local_background_width
        masks_vesselness_sigma
        masks_vesselness_beta
        masks_radius
        masks_minSize
        masks_diaphragmRadius
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
        entirePulseAnalysis
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
    end

    methods

        function obj = Parameters_json(dir_path)
            obj.path = dir_path;

            obj = obj.GetParameters();
        end

        function obj = GetParameters(obj)
            % Constructor method
            %[~, filename, ~] = fileparts(obj.path);
            filename_json = 'InputPulsewaveParams.json';
            dir_path_json = fullfile(obj.path, 'pulsewave', 'json');
            jsonPath = fullfile(dir_path_json, filename_json);

            if exist(jsonPath, 'file')
                jsonData = fileread(jsonPath);
                parsedData = jsondecode(jsonData);

                % Recherche de chaque paramètre
                obj.registerVideoFlag = parsedData.Video.Register;
                obj.refAvgStart = parsedData.Video.RefStart;
                obj.refAvgEnd = parsedData.Video.RefEnd;


                obj.videoStartFrameIndex = parsedData.Video.StartFrame;
                obj.videoEndFrameIndex = parsedData.Video.EndFrame;

                obj.frameWidth = parsedData.ResizeVideo.FrameWidth;
                obj.frameHeight = parsedData.ResizeVideo.FrameHeight;
                obj.videoLength = parsedData.ResizeVideo.VideoLength;

                obj.k = parsedData.ValueOfTheInterpolationParameter;
                obj.radius_ratio = parsedData.RadiusRatio;
                obj.radius_gap = parsedData.RadiusGap;
                obj.gauss_filt_size_for_barycentre = parsedData.GaussianFilterSizeForBarycentre;

                obj.veins_analysis = parsedData.VeinsAnalysis;
                obj.exportVideos = parsedData.ExportVideos;
                obj.entirePulseAnalysis = parsedData.EntirePulseAnalysis;

                obj.oneCycle_outNoiseThreshold = parsedData.CreationOfOneCycle.OutOfNoiseThreshold;
                obj.oneCycle_dataReliabilityThreshold = parsedData.CreationOfOneCycle.DataReliabilityIndexThreshold;
                obj.local_background_width = parsedData.CreationOfOneCycle.LocalBackgroundWidth;

                obj.masks_vesselness_sigma = parsedData.CreationOfMasks.VesselnesParameterSigma;
                obj.masks_vesselness_beta = parsedData.CreationOfMasks.VesselnesParameterBeta;
                obj.masks_radius = parsedData.CreationOfMasks.CropCoroidRadius;
                obj.masks_minSize = parsedData.CreationOfMasks.MinimumSeedAreaSize;
                obj.masks_diaphragmRadius = parsedData.CreationOfMasks.DiaphragmRadius;
                obj.CRACRV_Threshold = parsedData.CreationOfMasks.CentralRetinaVeinArteryThreshold;

                obj.cropSection_scaleFactorWidth = parsedData.CropSectionAnalysis.ScaleFactorWidth;
                obj.cropSection_scaleFactorSize = parsedData.CropSectionAnalysis.ScaleFactorSize;
                obj.cropSection_maskThreshold = parsedData.CropSectionAnalysis.MaskSliceThreshold;
                obj.cropSection_pixelSize = parsedData.CropSectionAnalysis.PixelSize; %??

                obj.elasWave_nDomFreq = parsedData.ElasticWaves.NumberOfDominantFrequency;
                obj.elasWave_pixelSize = parsedData.ElasticWaves.PixelSize;
                obj.elasWave_gaussFiltPadding = parsedData.ElasticWaves.GaussianFilterPadding;
                obj.elasWave_butterFiltOrder = parsedData.ElasticWaves.ButterworthFilterOrder;
                obj.elasWave_butterFiltBand = parsedData.ElasticWaves.ButterworthFilterThicknessBand;

                obj.systoleThreshold = parsedData.Systole.SystoleThreshold;

                obj.flatField_gwRatio = parsedData.FlatFieldCorrection.GWRatio;
                obj.flatField_border = parsedData.FlatFieldCorrection.Border;
                obj.flatField_borderDMap = parsedData.FlatFieldCorrection.BorderDMAP;
                obj.flatField_borderPulseAnal = parsedData.FlatFieldCorrection.BorderPulseAnalysis;

                obj.flowRate_gaussFiltPadding = parsedData.FlowRate.GaussianFilterPadding;
                obj.flowRate_yelorRangeHSV = parsedData.FlowRate.YellowOrangeRangeHSV;
                obj.flowRate_cydblueRangeHSV = parsedData.FlowRate.CyanDarkBlueRangeHSV;
                obj.flowRate_minTolVal = parsedData.FlowRate.TolerantValueMinimum;
                obj.flowRate_maxTolVal = parsedData.FlowRate.TolerantValueMaximum;
                obj.flowRate_sliceHalfThickness = parsedData.FlowRate.SliceHalfThickness;

                obj.normFlag = parsedData.OpticalPowerNormalization.NormalizationFlag;
                obj.normPowerCalibrationSlope = parsedData.OpticalPowerNormalization.PowerCalibrationCurveSlope;
                obj.normPoweryIntercept = parsedData.OpticalPowerNormalization.PowerCalibrationYIntercept;
                obj.normRefBloodFlow = parsedData.OpticalPowerNormalization.ReferenceTotalRetinalBloodFlow;

                obj.pupilRadius = parsedData.PulseAnalysis.RadiusPupil;
                obj.iris2retinaDist = parsedData.PulseAnalysis.DistanceIrisRetina;
                obj.theta = parsedData.PulseAnalysis.Theta;
                obj.opticalIndex = parsedData.PulseAnalysis.OpticalIndex;
                obj.lambda = parsedData.PulseAnalysis.Lambda;
                obj.phi = parsedData.PulseAnalysis.Phi; % Mie scattering angle
                obj.pulseAnal_dataReliabilityFactor = parsedData.PulseAnalysis.DataReliatibilityFactor; %?
                obj.pulseAnal_peakHeightThreshold = parsedData.PulseAnalysis.PeakHeightThreshold;
                obj.pulseAnal_outNoiseThreshold = parsedData.PulseAnalysis.OutOfNoiseThreshold;
                obj.pulseAnal_blurScaleFactor = parsedData.PulseAnalysis.BlurScaleFactor;
                obj.pulseAnal_frameMinDiastole = parsedData.PulseAnalysis.FramePercentageBeforeMinimumOfDiastole;
                obj.pulseAnal_framePeakSystole = parsedData.PulseAnalysis.FramePercentageAroundPeakSystole;
                obj.pulseAnal_exp = parsedData.PulseAnalysis.Exponentiel; %??

                obj.trWavelength_f0 = parsedData.Wavelength.F0;
                obj.trWavelength_r0 = parsedData.Wavelength.R0;

                obj.velocitySmallRadiusRatio = parsedData.Velocity.SmallRadiusRatio;
                obj.velocityBigRadiusRatio = parsedData.Velocity.BigRadiusRatio;

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

                obj.nbSides = parsedData.Other.NumberOfSides;
                obj.DiffFirstCalculationsFlag = parsedData.Other.DiffFirstCalculationsFlag;
                obj.AllCirclesFlag = parsedData.Other.AllCircles;
                obj.nbCircles = parsedData.Other.NumberOfCircles ;
                obj.forcewidth = parsedData.Other.ForceWidthInPixels;
                obj.forcebarycenter = parsedData.Other.ForceBarycenter;


            else
                error('The json file could not be found.');
            end

        end

    end

end
