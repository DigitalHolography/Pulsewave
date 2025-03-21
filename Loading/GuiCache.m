classdef GuiCache
% A parameter cache to store parameters from GUI
% so that during computations, values are fetched from cache instead of GUI
% so that the user can modify the values in the GUI without messing
% everything up
properties (Access = public)

    nb_cpu_cores double
    batch_size double
    spatialTransformation char
    z double
    z_retina double
    z_iris double
    z_switch char
    time_transform struct % object with : type of transformation, f1, f2
    blur double
    imageChoice char
    wavelength double
    Fs double
    pix_width double
    pix_height double
    parallelism char

    aberration_compensation logical
    notes cell
    position_in_file double
    output_videos char
    SVD logical

    % video rendering

    ref_batch_size double
    batch_stride double
    DX double
    DY double

    % video rendering logical checkboxes
    low_memory logical
    save_raw logical
    rephasing logical
    registration logical
    registration_via_phase logical
    iterative_registration logical
    temporal_filter_flag logical
    temporal_filter double
    registration_disc logical
    registration_disc_ratio double

    % color image parameters
    color_f1 double
    color_f2 double
    color_f3 double
    low_frequency logical

    % dark field parameters
    xystride double
    num_unit_cells_x double
    r1 double

    % iterative optimization parameters
    iterative_aberration_compensation double
    optimization_zernike_ranks double
    mask_num_iter double
    zernikes_tol double
    max_constraint double

    % shack-hartmann parameters
    shack_hartmann_aberration_compensation double
    shack_hartmann_zernike_ranks double
    image_subapertures_size_ratio double
    num_subapertures_positions double
    subaperture_margin double
    SubAp_PCA double
    minSubAp_PCA double
    maxSubAp_PCA double
    zernike_projection double
    shack_hartmann_ref_image double

end

methods (Access = public)

    function obj = GuiCache(app)
        obj.nb_cpu_cores = app.numworkersSpinner.Value;
        obj.batch_size = app.batchsizeEditField.Value;
        obj.ref_batch_size = app.refbatchsizeEditField.Value;
        obj.batch_stride = app.batchstrideEditField.Value;
        obj.spatialTransformation = app.spatialTransformationDropDown.Value; % Fresnel or angular spectrum...&obj.
        obj.z_switch = app.Switch.Value;
        obj.z_retina = app.zretinaEditField.Value;
        obj.z_iris = app.zirisEditField.Value;
        obj.time_transform = app.time_transform;
        obj.blur = app.blurEditField.Value;
        obj.imageChoice = strrep(app.ImageChoiceDropDown.Value, ' ', '_');
        obj.wavelength = app.wavelengthEditField.Value;
        obj.Fs = app.Fs;
        obj.pix_width = app.pix_width;
        obj.pix_height = app.pix_height;
        obj.registration = app.imageregistrationCheckBox.Value;
        obj.temporal_filter_flag = app.temporalfilterCheckBox.Value;
        obj.temporal_filter = app.temporalfilterEditField.Value;
        obj.parallelism = app.ParallelismDropDown.Value;
        obj.registration_via_phase = app.phaseregistrationCheckBox.Value;
        obj.aberration_compensation = app.IterativeoptimizationCheckBox.Value;
        obj.SVD = app.SVDCheckBox.Value;

        obj.notes = app.NotesTextArea.Value;
        obj.DX = app.DX;
        obj.DY = app.DY;
        obj.position_in_file = app.positioninfileSlider.Value;
        obj.output_videos = (strrep(app.outputvideoDropDown.Value, ' ', '_'));
        obj.low_memory = app.lowmemoryCheckBox.Value;
        obj.rephasing = app.rephasingCheckBox.Value;
        obj.save_raw = app.saverawvideosCheckBox.Value;
        obj.registration_disc = app.registrationdiscCheckBox.Value;
        obj.registration_disc_ratio = app.regDiscRatioEditField.Value;

        obj.color_f1 = app.compositef1EditField.Value;
        obj.color_f2 = app.compositef2EditField.Value;
        obj.color_f3 = app.compositef3EditField.Value;
        obj.low_frequency = app.lowfrequencyCheckBox.Value;

        % dark field parameters
        obj.xystride = app.xystrideEditField.Value;
        obj.num_unit_cells_x = app.unitcellsinlatticeEditField.Value;
        obj.r1 = app.r1EditField.Value;

        % iterative optimization parameters
        obj.iterative_aberration_compensation = app.IterativeoptimizationCheckBox.Value;
        obj.optimization_zernike_ranks = app.optimizationzernikeranksEditField.Value;
        obj.mask_num_iter = app.masknumiterEditField.Value;
        obj.zernikes_tol = app.zernikestolEditField.Value;
        obj.max_constraint = app.maxconstraintEditField.Value;
        obj.iterative_registration = app.iterativeregistrationCheckBox.Value;

        % shack-hartmann parameters
        obj.shack_hartmann_aberration_compensation = app.ShackHartmannCheckBox.Value;
        obj.shack_hartmann_zernike_ranks = app.shackhartmannzernikeranksEditField.Value;
        obj.image_subapertures_size_ratio = app.imagesubapsizeratioEditField.Value;
        obj.num_subapertures_positions = app.subapnumpositionsEditField.Value;
        obj.subaperture_margin = app.subaperturemarginEditField.Value;
        obj.SubAp_PCA = app.SubAp_PCACheckBox.Value;
        obj.minSubAp_PCA = app.minSubAp_PCAEditField.Value;
        obj.maxSubAp_PCA = app.maxSubAp_PCAEditField.Value;
        obj.zernike_projection = app.ZernikeProjectionCheckBox.Value;
        obj.shack_hartmann_ref_image = app.referenceimageDropDown.Value;

        % bufferize (and lock during computation) current paremeter values from front end
    end

    function load2Gui(obj, app)
        % set gui parameters from cache

        app.batchsizeEditField.Value = loadGUIVariable(obj.batch_size);
        app.max_PCAEditField.Limits = loadGUIVariable([0 double(app.batchsizeEditField.Value)]);
        app.spatialTransformationDropDown.Value = loadGUIVariable(obj.spatialTransformation);
        app.refbatchsizeEditField.Value = loadGUIVariable(obj.ref_batch_size);
        app.batchstrideEditField.Value = loadGUIVariable(obj.batch_stride);
        app.Switch.Value = loadGUIVariable(obj.z_switch);
        app.zretinaEditField.Value = loadGUIVariable(obj.z_retina);
        app.zirisEditField.Value = loadGUIVariable(obj.z_iris);
        app.time_transform = loadGUIVariable(obj.time_transform);
        app.blurEditField.Value = loadGUIVariable(obj.blur);
        app.ImageChoiceDropDown.Value = loadGUIVariable(strrep(obj.imageChoice, '_', ' '));
        app.timetransformDropDown.Value = loadGUIVariable(obj.time_transform.type);
        app.SVDCheckBox.Value = loadGUIVariable(obj.SVD);
        app.f1EditField.Value = loadGUIVariable(obj.time_transform.f1);
        app.f2EditField.Value = loadGUIVariable(obj.time_transform.f2);
        app.min_PCAEditField.Value = loadGUIVariable(obj.time_transform.min_PCA);
        app.max_PCAEditField.Value = loadGUIVariable(obj.time_transform.max_PCA);
        app.imageregistrationCheckBox.Value = loadGUIVariable(obj.registration);
        app.iterativeregistrationCheckBox.Value = loadGUIVariable(obj.iterative_registration);
        app.wavelengthEditField.Value = loadGUIVariable(obj.wavelength);
        app.outputvideoDropDown.Value = loadGUIVariable(strrep(obj.output_videos, '_', ' '));
        app.lowmemoryCheckBox.Value = loadGUIVariable(obj.low_memory);
        app.rephasingCheckBox.Value = loadGUIVariable(obj.rephasing);
        app.NotesTextArea.Value = loadGUIVariable(obj.notes);
        app.positioninfileSlider.Value = loadGUIVariable(obj.position_in_file);
        app.EditField.Value = loadGUIVariable(obj.position_in_file);
        app.IterativeoptimizationCheckBox.Value = loadGUIVariable(obj.iterative_aberration_compensation);
        app.masknumiterEditField.Value = loadGUIVariable(obj.mask_num_iter);
        app.ShackHartmannCheckBox.Value = loadGUIVariable(obj.shack_hartmann_aberration_compensation);
        app.ZernikeProjectionCheckBox.Value = loadGUIVariable(obj.zernike_projection);
        app.referenceimageDropDown.Value = loadGUIVariable(obj.shack_hartmann_ref_image);
        app.shackhartmannzernikeranksEditField.Value = loadGUIVariable(obj.shack_hartmann_zernike_ranks);
        app.imagesubapsizeratioEditField.Value = loadGUIVariable(obj.image_subapertures_size_ratio);
        app.subapnumpositionsEditField.Value = loadGUIVariable(obj.num_subapertures_positions);
        app.subaperturemarginEditField.Value = loadGUIVariable(obj.subaperture_margin);
        app.SubAp_PCACheckBox.Value = loadGUIVariable(obj.SubAp_PCA);
        app.minSubAp_PCAEditField.Value = loadGUIVariable(obj.minSubAp_PCA);
        app.maxSubAp_PCAEditField.Value = loadGUIVariable(obj.maxSubAp_PCA);
        app.compositef1EditField.Value = loadGUIVariable(obj.color_f1);
        app.compositef2EditField.Value = loadGUIVariable(obj.color_f2);
        app.compositef3EditField.Value = loadGUIVariable(obj.color_f3);
        app.lowfrequencyCheckBox.Value = loadGUIVariable(obj.low_frequency);
        app.saverawvideosCheckBox.Value = loadGUIVariable(obj.save_raw);

    end

end

end
