function mask = getHolonetprediction(M0, net)

if nargin < 2 || isempty(net)
    net = importONNXNetwork("C:\Users\Mikhalkino\Downloads\onnx_model.onnx");
end

[Nx, Ny] = size(M0);

M0 = imresize(rescale(M0), [400, 400]);

input = repmat(255 * rescale(M0), 1, 1, 3);

output = predict(net, input);

mask = sigmoid(output(:, :, 2) / 200) > 0.6;

mask = imresize(mask, [Nx, Ny]);
end
