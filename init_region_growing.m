function [mask] = init_region_growing(img,nb_div)


[N,M] = size(img);

mask = zeros(N,M);

for i=1:floor(N/nb_div):N-floor(N/nb_div)
    for j=1:floor(M/nb_div):M-floor(M/nb_div)
        
        window = img(i:i+floor(N/nb_div), j:j+floor(M/nb_div));

        high_int = window > (0.5*max(window,[],'all'));

        high_int_one_region = bwareafilt(high_int,1);

        mask(i:i+floor(N/nb_div), j:j+floor(M/nb_div)) = high_int_one_region;


    end
end


end

