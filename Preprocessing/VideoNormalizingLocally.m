function obj = VideoNormalizingLocally(obj)

    PW_params = Parameters_json(obj.directory);

    obj.f_RMS_video = sqrt(double(obj.M2_data_video) ./ double(obj.M0_data_video));
    obj.f_AVG_video = double(obj.M1_data_video) ./ double(obj.M0_data_video);
    obj.M0_disp_video = flat_field_correction(obj.M0_disp_video, ceil(PW_params.flatField_gwRatio * size(obj.M0_disp_video, 1)), PW_params.flatField_border);

end