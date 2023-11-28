function [json_data] = txt2json_param(fileContent)

parsedData = {};

parsedData.Video.StartFrame = extractValue(fileContent, 'Start frame :');
parsedData.Video.EndFrame = extractValue(fileContent, 'End frame :');

parsedData.ValueOfTheInterpolationParameter = extractValue(fileContent, 'Value of the interpolation parameter :');
parsedData.RadiusRatio = extractValue(fileContent, 'Radius ratio :');
parsedData.RadiusGap = extractValue(fileContent, 'Radius gap :');

parsedData.ResistivityIndex.GaussianFilterSigma = extractValue(fileContent, 'Gaussian filter sigma :');
parsedData.ResistivityIndex.AmplitudeSaturationPonderation = extractValue(fileContent, 'Amplitude saturation  ponderation :');
parsedData.ResistivityIndex.TolerantValueMinimum = extractValue(fileContent, 'Tolerante value minimum :');
parsedData.ResistivityIndex.TolerantValueMaximum = extractValue(fileContent, 'Tolerante value maximum :');
parsedData.ResistivityIndex.Gamma = extractValue(fileContent, 'Gamma :');

parsedData.CreationOfOneCycle.NumberOfInterpolation = extractValue(fileContent, 'Number of interpolation :');
parsedData.CreationOfOneCycle.OutOfNoiseThreshold = extractValue(fileContent, 'Out of noise treshold :');
parsedData.CreationOfOneCycle.DataReliabilityIndexThreshold = extractValue(fileContent, 'Data reliability index treshold :');
parsedData.CreationOfOneCycle.LocalBackgroundWidth = extractValue(fileContent, 'Local background width :');

parsedData.CreationOfMasks.VesselnesParameterSigma = extractValue(fileContent, 'Vesselnes parameter sigma :');
parsedData.CreationOfMasks.VesselnesParameterBeta = extractValue(fileContent, 'Vesselnes parameter beta :');
parsedData.CreationOfMasks.VesselnessBinarizationThreshold = extractValue(fileContent, 'Vesselness Binarization treshold :');
parsedData.CreationOfMasks.MagicWandNumberOfSegmentedAreaDetectedForVessels = extractValue(fileContent, 'MagicWand Number of segmented area detected for vessels :');
parsedData.CreationOfMasks.MagicWandNumberOfSegmentedAreaDetectedForArtery = extractValue(fileContent, 'MagicWand Number of segmented area detected for artery :');
parsedData.CreationOfMasks.ArteryCorrelationThreshold = extractValue(fileContent, 'Artery correlation threshold :');

parsedData.CentralRetinaMask.CentralRetinaArteryThreshold = extractValue(fileContent, 'Central retina artery treshold :');
parsedData.CentralRetinaMask.CentralRetinaVeinThreshold = extractValue(fileContent, 'Central retina vein treshold :');
parsedData.CentralRetinaMask.CentralRetinaBackgroundThreshold = extractValue(fileContent, 'Central retina background treshold :');
parsedData.CentralRetinaMask.StandardBinarizationThreshold = extractValue(fileContent, 'Standard binarization treshold :');
parsedData.CentralRetinaMask.CropCoroidRadius = extractValue(fileContent, 'Crop coroid radius :');
parsedData.CentralRetinaMask.MinimumSeedAreaSize = extractValue(fileContent, 'Minimum seed area size :');
parsedData.CentralRetinaMask.CleaningCoroidOrNot = extractValue(fileContent, 'Cleaning Coroid or not :');
parsedData.CentralRetinaMask.ShowingIntermediateFiguresInTheProcess = extractValue(fileContent, 'Showing intermediate figures in the process :');
parsedData.CentralRetinaMask.RegionGrowingFloorThreshold = extractValue(fileContent, 'Region growing floor threshold :');
parsedData.CentralRetinaMask.RegionGrowingVeinConditionThreshold = extractValue(fileContent, 'Region growing vein condition threshold :');
parsedData.CentralRetinaMask.RegionGrowingArteryConditionThreshold = extractValue(fileContent, 'Region growing artery condition threshold :');
parsedData.CentralRetinaMask.CentralRetinaVeinArteryThreshold = extractValue(fileContent, 'Central Retina Vein/Artery threshold :');
parsedData.CentralRetinaMask.RegionGrowingThreshold = extractValue(fileContent, 'Region growing threshold :');

parsedData.CropSectionAnalysis.ScaleFactorWidth = extractValue(fileContent, 'Scale factor width :');
parsedData.CropSectionAnalysis.ScaleFactorSize = extractValue(fileContent, 'Scale factor size :');
parsedData.CropSectionAnalysis.MaskSliceThreshold = extractValue(fileContent, 'Mask slice treshold :');
parsedData.CropSectionAnalysis.PixelSize = extractValue(fileContent, 'Pixel size ? :');

parsedData.ElasticWaves.NumberOfDominantFrequency = extractValue(fileContent, 'Number of dominant frequency :');
parsedData.ElasticWaves.PixelSize = extractValue(fileContent, 'Pixel size :');
parsedData.ElasticWaves.GaussianFilterPadding = extractValue(fileContent, 'Gaussian filter padding :');
parsedData.ElasticWaves.ButterworthFilterOrder = extractValue(fileContent, 'Butterworth filter order :');
parsedData.ElasticWaves.ButterworthFilterThicknessBand = extractValue(fileContent, 'Butterworth filter thickness band :');

parsedData.Systole.SystoleThreshold = extractValue(fileContent, 'Systole threshold :');

parsedData.FlatFieldCorrection.GWRatio = extractValue(fileContent, 'GW ratio :');
parsedData.FlatFieldCorrection.Border = extractValue(fileContent, 'Border :');
parsedData.FlatFieldCorrection.BorderDMAP = extractValue(fileContent, 'Border dMAP :');
parsedData.FlatFieldCorrection.BorderPulseAnalysis = extractValue(fileContent, 'Border PulseAnalisys :');

parsedData.FlowRate.GaussianFilterPadding = extractValue(fileContent, 'Gaussian filter padding :');
parsedData.FlowRate.YellowOrangeRangeHSV = extractValue(fileContent, 'Yellow-orange range HSV :');
parsedData.FlowRate.CyanDarkBlueRangeHSV = extractValue(fileContent, 'Cyan-dark blue range HSV :');
parsedData.FlowRate.TolerantValueMinimum = extractValue(fileContent, 'Tolerante value minimum :');
parsedData.FlowRate.TolerantValueMaximum = extractValue(fileContent, 'Tolerante value maximum :');
parsedData.FlowRate.SliceHalfThickness = extractValue(fileContent, 'Slice half thickness :');

parsedData.PulseAnalysis.RadiusPupil = extractValue(fileContent, 'Radius pupil :');
parsedData.PulseAnalysis.DistanceIrisRetina = extractValue(fileContent, 'Distance iris-retina :');
parsedData.PulseAnalysis.Theta = extractValue(fileContent, 'Theta :');
parsedData.PulseAnalysis.OpticalIndex = extractValue(fileContent, 'Optical index :');
parsedData.PulseAnalysis.Lambda = extractValue(fileContent, 'Lambda :');
parsedData.PulseAnalysis.DataReliatibilityFactor = extractValue(fileContent, 'Data reliatibility factor ? :');
parsedData.PulseAnalysis.PeakHeightThreshold = extractValue(fileContent, 'Peak height treshold :');
parsedData.PulseAnalysis.OutOfNoiseThreshold = extractValue(fileContent, 'Out of noise treshold :');
parsedData.PulseAnalysis.BlurScaleFactor = extractValue(fileContent, 'Blur scale factor :');
parsedData.PulseAnalysis.FramePercentageBeforeMinimumOfDiastole = extractValue(fileContent, 'Frame percentage before minimum of diastole :');
parsedData.PulseAnalysis.FramePercentageAroundPeakSystole = extractValue(fileContent, 'Frame percentage around peak systole :');
parsedData.PulseAnalysis.Exponentiel = extractValue(fileContent, 'exponentiel 1 ? :');

parsedData.Wavelength.F0 = extractValue(fileContent, 'F0 :');
parsedData.Wavelength.R0 = extractValue(fileContent, 'R0 :');

parsedData.VesselsVideo.RadiusRatioFactor = extractValue(fileContent, 'Radius ratio factor :');
parsedData.VesselsVideo.GaussianFilterPadding = extractValue(fileContent, 'Gaussian filter padding :');
parsedData.VesselsVideo.GaussianFilterFactor = extractValue(fileContent, 'Gaussian filter factor :');

parsedData.Viscosity.NumberOfInterpolation = extractValue(fileContent, 'Number of interpolation :');
parsedData.Viscosity.InterpolationParameter = extractValue(fileContent, 'Interpolation parameter :');
parsedData.Viscosity.InterpolationSizeFactor = extractValue(fileContent, 'Interpolation size factor :');
parsedData.Viscosity.ImageProjectionThreshold = extractValue(fileContent, 'Image projection treshold :');
parsedData.Viscosity.ViscosityFunctionFactor = extractValue(fileContent, 'Viscosity function factor :');
parsedData.Viscosity.FitFunctionParameterMatrix = extractValue(fileContent, 'Fit function parameter (matrix) :');
parsedData.Viscosity.ViscosityListParameterA = extractValue(fileContent, 'Viscosity list parameter a :');
parsedData.Viscosity.ViscosityListParameterB = extractValue(fileContent, 'Viscosity list parameter b :');

parsedData.Other.NumberOfSides = extractValue(fileContent, 'Number of sides :');


fields = fieldnames(parsedData);
json_data = rmfield(parsedData, fields(structfun(@isempty, parsedData)));

end