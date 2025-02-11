function obj = VideoNormalizingLocally(obj)

PW_params = Parameters_json(obj.directory,obj.PW_param_name);

[N, M, F] = size(obj.M0_data_video);
alpha = PW_params.alphaConvolveNorm;
D = (M + N) / 2;

if alpha == 1
    % behaves as if conv_size = alpha*(2*D-1) just faster;
    M0_data_convoluated = double(mean(obj.M0_data_video, [1, 2]));
elseif alpha == 0
    % forces the pixel M0 normalisation;
    M0_data_convoluated = double(obj.M0_data_video);
else
    conv_size = round(alpha * (2 * D - 1));
    M0_data_convoluated = zeros([N, M, F]);
    conv_kern = ones(conv_size);
    
    parfor i = 1:F
        M0_data_convoluated(:, :, i) = conv2(double(obj.M0_data_video(:, :, i)), conv_kern, 'same');
    end
    
    S = sum(obj.M0_data_video, [1, 2]);
    S2 = sum(M0_data_convoluated, [1, 2]);

    
    imwrite(rescale(mean(M0_data_convoluated,3)),fullfile(obj.directory,'pulsewave', sprintf("%s_alpha=%s_%s", obj.filenames, num2str(alpha), 'M0_Convolution_Norm.png')), 'png');
    
    M0_data_convoluated = M0_data_convoluated .* S ./ S2; % normalizing to get the average with alpha = 0;
end

if PW_params.NormTempMode
        M0_data_convoluated = mean(M0_data_convoluated, 3);
end

obj.f_RMS_video = sqrt(double(obj.M2_data_video) ./ M0_data_convoluated);
obj.f_AVG_video = double(obj.M1_data_video) ./ M0_data_convoluated;
obj.M0_ff_video = flat_field_correction(obj.M0_ff_video, ceil(PW_params.flatField_gwRatio * size(obj.M0_ff_video, 1)), PW_params.flatField_border);

end
