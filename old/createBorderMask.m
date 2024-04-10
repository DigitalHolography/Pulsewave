function borderMask = createBorderMask(image, borderAmount)

if borderAmount == 0
    a = 1;
    b = size(image,1);
    c = 1;
    d = size(image,2);
else
    a = ceil(size(image,1)*borderAmount);
    b = floor(size(image,1)*(1-borderAmount));
    c = ceil(size(image,2)*borderAmount);
    d = floor(size(image,2)*(1-borderAmount));
end

borderMask = zeros(size(image));
borderMask(a:b,c:d) = 1;

end

