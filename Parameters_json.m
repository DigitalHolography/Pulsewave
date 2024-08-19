classdef Parameters_json
    % Class for storing parameters from an input json file

    properties
        path
        registerVideoFlag
        videoStartFrameIndex
        videoEndFrameIndex
        frameWidth
        frameHeight
        videoLength
        k
        radius_ratio
        radius_gap
        gauss_filt_size_for_barycentre
        resistivity_gaussFiltSigma
        resistivity_satAmp
        resistivity_minTolVal
        resistivity_maxTolVal
        resistivity_gamma
        oneCycle_Ninterp
        oneCycle_outNoiseThreshold
        oneCycle_dataReliabilityThreshold
        local_background_width
        arteryMask_vesselness_sigma
        arteryMask_vesselness_beta
        arteryMask_vesselness_bin_threshold
        arteryMask_magicwand_nb_of_area_vessels
        arteryMask_magicwand_nb_of_area_artery
        arteryMask_ArteryCorrThreshold
        arteryMask_CorrelationMatrixThreshold
        arteryMask_vesselnessContrastNumSlides
        centralRetinal_arteryThreshold
        centralRetinal_veinThreshold
        centralRetinal_backgndThreshold
        vesselMask_BinTreshold
        masks_radius
        masks_radius_treshold
        masks_minSize
        masks_cleaningCoroid
        masks_showIntermediateFigures
        RG_FloorThreshold
        RG_veinConditionThreshold
        RG_ArteryConditionThreshold
        RG_ArterySeedsThreshold
        CRACRV_Threshold
        RG_alpha
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
        velocity_smallRadiusRatio
        velocity_bigRadiusRatio
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
    end

    methods

        function obj = Parameters_json(dir_path)
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
            filename_json = 'InputPulsewaveParams.json';
            dir_path_json = fullfile(obj.path, 'pulsewave', 'json');
            jsonPath = fullfile(dir_path_json, filename_json);

            if exist(jsonPath, 'file')
                jsonData = fileread(jsonPath);
                parsedData = jsondecode(jsonData);

                % Recherche de chaque paramètre
                obj.registerVideoFlag = parsedData.Video.Register;

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

                obj.resistivity_gaussFiltSigma = parsedData.ResistivityIndex.GaussianFilterSigma;
                obj.resistivity_satAmp = parsedData.ResistivityIndex.AmplitudeSaturationPonderation;
                obj.resistivity_minTolVal = parsedData.ResistivityIndex.TolerantValueMinimum;
                obj.resistivity_maxTolVal = parsedData.ResistivityIndex.TolerantValueMaximum;
                obj.resistivity_gamma = parsedData.ResistivityIndex.Gamma;

                obj.oneCycle_Ninterp = parsedData.CreationOfOneCycle.NumberOfInterpolation;
                obj.oneCycle_outNoiseThreshold = parsedData.CreationOfOneCycle.OutOfNoiseThreshold;
                obj.oneCycle_dataReliabilityThreshold = parsedData.CreationOfOneCycle.DataReliabilityIndexThreshold;
                obj.local_background_width = parsedData.CreationOfOneCycle.LocalBackgroundWidth;

                obj.arteryMask_vesselness_sigma = parsedData.CreationOfMasks.VesselnesParameterSigma;
                obj.arteryMask_vesselness_beta = parsedData.CreationOfMasks.VesselnesParameterBeta;
                obj.arteryMask_vesselness_bin_threshold = parsedData.CreationOfMasks.VesselnessBinarizationThreshold;
                obj.arteryMask_magicwand_nb_of_area_vessels = parsedData.CreationOfMasks.MagicWandNumberOfSegmentedAreaDetectedForVessels;
                obj.arteryMask_magicwand_nb_of_area_artery = parsedData.CreationOfMasks.MagicWandNumberOfSegmentedAreaDetectedForArtery;
                obj.arteryMask_ArteryCorrThreshold = parsedData.CreationOfMasks.ArteryCorrelationThreshold;
                obj.arteryMask_CorrelationMatrixThreshold = parsedData.CreationOfMasks.ArteryCorrelationMatrixThreshold;
                obj.arteryMask_vesselnessContrastNumSlides = parsedData.CreationOfMasks.ContrastNumSlides;

                obj.centralRetinal_arteryThreshold = parsedData.CentralRetinaMask.CentralRetinaArteryThreshold;
                obj.centralRetinal_veinThreshold = parsedData.CentralRetinaMask.CentralRetinaVeinThreshold;
                obj.centralRetinal_backgndThreshold = parsedData.CentralRetinaMask.CentralRetinaBackgroundThreshold;
                obj.vesselMask_BinTreshold = parsedData.CentralRetinaMask.StandardBinarizationThreshold;
                obj.masks_radius = parsedData.CentralRetinaMask.CropCoroidRadius;
                obj.masks_radius_treshold = parsedData.CentralRetinaMask.TresholdRadiusValue;
                obj.masks_minSize = parsedData.CentralRetinaMask.MinimumSeedAreaSize;
                obj.masks_cleaningCoroid = parsedData.CentralRetinaMask.CleaningCoroidOrNot;
                obj.masks_showIntermediateFigures = parsedData.CentralRetinaMask.ShowingIntermediateFiguresInTheProcess;
                obj.RG_FloorThreshold = parsedData.CentralRetinaMask.RegionGrowingFloorThreshold;
                obj.RG_veinConditionThreshold = parsedData.CentralRetinaMask.RegionGrowingVeinConditionThreshold;
                obj.RG_ArteryConditionThreshold = parsedData.CentralRetinaMask.RegionGrowingArteryConditionThreshold;
                obj.RG_ArterySeedsThreshold = parsedData.CentralRetinaMask.RegionGrowingArterySeedsThreshold;
                obj.CRACRV_Threshold = parsedData.CentralRetinaMask.CentralRetinaVeinArteryThreshold;
                obj.RG_alpha = parsedData.CentralRetinaMask.RegionGrowingThreshold;

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

                obj.velocity_smallRadiusRatio = parsedData.Velocity.SmallRadiusRatio;
                obj.velocity_bigRadiusRatio = parsedData.Velocity.BigRadiusRatio;

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

            else
                error('The json file could not be found.');
            end

        end

    end

end
