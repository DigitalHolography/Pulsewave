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
            obj.minSubAp_PCA = app.minSubAp_PCAEditField_2.Value;
            obj.maxSubAp_PCA = app.maxSubAp_PCAEditField_2.Value;
            obj.zernike_projection = app.ZernikeProjectionCheckBox.Value;
            obj.shack_hartmann_ref_image = app.referenceimageDropDown.Value;

            % bufferize (and lock during computation) current paremeter values from front end
        end

        function load2Gui(obj, app)
            % set gui parameters from cache
            try
                app.batchsizeEditField.Value = obj.batch_size;
                app.max_PCAEditField.Limits = [0 double(app.batchsizeEditField.Value)];
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            % try
            %     app.numworkersSpinner.Value = obj.nb_cpu_cores;
            % end

            try
                app.spatialTransformationDropDown.Value = obj.spatialTransformation;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.refbatchsizeEditField.Value = obj.ref_batch_size;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.batchstrideEditField.Value = obj.batch_stride;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.Switch.Value = obj.z_switch;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.zretinaEditField.Value = obj.z_retina;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.zirisEditField.Value = obj.z_iris;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.time_transform = obj.time_transform;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.blurEditField.Value = obj.blur;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.ImageChoiceDropDown.Value = strrep(obj.imageChoice, '_', ' ');
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.timetransformDropDown.Value = obj.time_transform.type;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.SVDCheckBox.Value = obj.SVD;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.f1EditField.Value = obj.time_transform.f1;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.f2EditField.Value = obj.time_transform.f2;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.min_PCAEditField.Value = obj.time_transform.min_PCA;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.max_PCAEditField.Value = obj.time_transform.max_PCA;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.imageregistrationCheckBox.Value = obj.registration;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.iterativeregistrationCheckBox.Value = obj.iterative_registration;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.wavelengthEditField.Value = obj.wavelength;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.outputvideoDropDown.Value = strrep(obj.output_videos, '_', ' ');
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.lowmemoryCheckBox.Value = obj.low_memory;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.rephasingCheckBox.Value = obj.rephasing;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.NotesTextArea.Value = obj.notes;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.positioninfileSlider.Value = obj.position_in_file;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.EditField.Value = obj.position_in_file;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.IterativeoptimizationCheckBox.Value = obj.iterative_aberration_compensation;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.masknumiterEditField.Value = obj.mask_num_iter;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            %             app.zernikestolEditField.Value = obj.low_order_zernikes;
            %             app.maxconstraintEditField.Value = obj.low_order_max_constraint;
            try
                app.ShackHartmannCheckBox.Value = obj.shack_hartmann_aberration_compensation;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.ZernikeProjectionCheckBox.Value = obj.zernike_projection;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.referenceimageDropDown.Value = obj.shack_hartmann_ref_image;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.shackhartmannzernikeranksEditField.Value = obj.shack_hartmann_zernike_ranks;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.imagesubapsizeratioEditField.Value = obj.image_subapertures_size_ratio;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.subapnumpositionsEditField.Value = obj.num_subapertures_positions;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.subaperturemarginEditField.Value = obj.subaperture_margin;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.SubAp_PCACheckBox.Value = obj.SubAp_PCA;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.minSubAp_PCAEditField_2.Value = obj.minSubAp_PCA;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.maxSubAp_PCAEditField_2.Value = obj.maxSubAp_PCA;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.compositef1EditField.Value = obj.color_f1;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.compositef2EditField.Value = obj.color_f2;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.compositef3EditField.Value = obj.color_f3;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.lowfrequencyCheckBox.Value = obj.low_frequency;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

            try
                app.saverawvideosCheckBox.Value = obj.save_raw;
            catch ME
                disp('Error Message:')
                disp(ME.message)
                for i = 1:numel(ME.stack)
                    ME.stack(i)
                end
            end

        end

    end

end
