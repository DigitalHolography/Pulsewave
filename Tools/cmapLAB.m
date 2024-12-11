function [cmap] = cmapLAB(color1, color2, n, colorn, posn)

arguments
    color1
    color2
    n = 256
end

arguments (Repeating)
    colorn
    posn
end

numArgs = size(colorn, 1);

for i = 1:numArgs

lab1 = rgb2lab(color1);
lab2 = rgb2lab(color2);

dL = lab2(1) - lab1(1);
da = lab2(2) - lab1(2);
db = lab2(3) - lab1(3);

x = linspace(0, 1, n);

L = lab1(1) + dL * x;
a = lab1(2) + da * x;
b = lab1(3) + db * x;

cmap = zeros(n, 3);

for idx = 1:n
    cmap(idx, :) = lab2rgb([L(idx) a(idx) b(idx)]);
end

end