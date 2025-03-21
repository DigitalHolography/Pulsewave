function [cmap] = cmapLAB(n, color_n, pos_n)

arguments
    n = 256
end

arguments (Repeating)
    color_n
    pos_n
end

pos_n{1} = 0;
pos_n{end} = 1;

numArgs = size(color_n, 2);

L = zeros(n, 1);
a = zeros(n, 1);
b = zeros(n, 1);

for i = 1:numArgs - 1

    n1 = round(n * pos_n{i} + 1);
    n2 = round(n * pos_n{i + 1});

    color1 = color_n{i};
    color2 = color_n{i + 1};

    lab1 = rgb2lab(color1);
    lab2 = rgb2lab(color2);

    dL = lab2(1) - lab1(1);
    da = lab2(2) - lab1(2);
    db = lab2(3) - lab1(3);

    x = linspace(0, 1, n2 - n1 + 1);

    L(n1:n2) = lab1(1) + dL * x;
    a(n1:n2) = lab1(2) + da * x;
    b(n1:n2) = lab1(3) + db * x;

end

cmap = zeros(n, 3);

for idx = 1:n
    cmap(idx, :) = lab2rgb([L(idx) a(idx) b(idx)]);
end

cmap(cmap > 1) = 1;
cmap(cmap < 0) = 0;
cmap(isnan(cmap)) = 0;

end
