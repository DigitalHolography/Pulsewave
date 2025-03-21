function obj = VideoNormalizing(obj)

params = Parameters_json(obj.directory, obj.param_name);
M0_data_mean = mean(double(obj.M0_data_video), [1 2]);

obj.f_RMS_video = sqrt(double(obj.M2_data_video) ./ M0_data_mean);
obj.f_AVG_video = double(obj.M1_data_video) ./ M0_data_mean;
obj.M0_ff_video = flat_field_correction(obj.M0_ff_video, ceil(params.flatField_gwRatio * size(obj.M0_ff_video, 1)), params.flatField_border);

end
