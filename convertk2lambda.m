function [x, profile] = convertk2lambda(profile, pixel_size)

x = zeros(length(profile), 1);

for i = 1 : length(profile)
    x(i) = (length(profile) - i)*pixel_size;
end

end