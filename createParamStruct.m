function param = createParamStruct()

% arteries mask
field_1 = 'video_std_binarization_threshold'; value_1 = 0.2;
field_4 = 'mask_artery_choroid_correlation_threshold'; value_4 = 0.2;
field_5 = 'magic_wand_num_of_patterns_arteries'; value_5 = 4; %% additional default parameters

% peak detection
field_6 = 'find_peaks_in_pulse_threshold'; value_6 = 0.7;

% vessels mask
field_14 = 'vessels_mask_binarization_threshold'; value_14 = 0.4;

% veins mask
field_13 = 'magic_wand_num_of_patterns_veins'; value_13 = 6;

% flat field parameters
field_2 = 'flat_field_gw_ratio'; value_2 = 0.7;
field_3 = 'flat_field_border'; value_3 = 0.25;
field_15 = 'flat_field_border_dMap'; value_15 = 0.33;
field_7 = 'flat_field_border_pulseAnalysis'; value_7 = 0;

% physical parameters
field_8 = 'pupil_radius'; value_8 = 2.0; %mm
field_9 = 'iris_to_retina_distance'; value_9 = 20.0; %mm
field_10 = 'theta'; value_10 = 0.05; %rad
field_11 = 'optical_index'; value_11  = 1.35;
field_12 = 'lambda'; value_12 = 852e-9;


% additional parameters in index reliability index

%it takes a lot of time, would it be possible to unify the structure : have
%only the image field or only the additional images?

param = struct(field_1, value_1, field_2, value_2, field_3, value_3, field_4, value_4, field_5, value_5, field_6, value_6, ...
    field_7, value_7, field_8, value_8, field_9, value_9, field_10, value_10, field_11, value_11, field_12, value_12, ...
    field_13, value_13, field_14, value_14, field_15, value_15);

end