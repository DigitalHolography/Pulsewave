function matrix = video2vessels(v_RMS, maskArtery, maskCRA)

nb_sides = 120;

matrix = zeros(100, nb_sides/10);

% why radius ratio if this is effectively an absolute number?
radius_ratio = round(0.27* size(v_RMS,1));
%FIXME : anamorphic image
blurred_mask = imgaussfilt(double(mean(v_RMS,3).*double(maskCRA)),round(size(maskCRA,1)/4),'Padding',0);
[~,x_center] = findpeaks(sum(blurred_mask,1));
[~,y_center] = findpeaks(sum(blurred_mask,2));


for rr = 0 : 9
    radius_ratio = radius_ratio + 1;
    polygon = nsidedpoly(nb_sides, 'Center', [x_center, y_center], 'Radius', radius_ratio);
    points_x = polygon.Vertices(:,1);
    points_x(end + 1) = points_x(1);
    points_y = polygon.Vertices(:,2);
    points_y(end + 1) = points_y(1);
    % figure(121)
    % for ii = 1:nb_sides
    %     l   = line([points_x(ii), points_x(ii + 1)], [points_y(ii), points_y(ii + 1)]);
    %     l.Color = 'red';
    %     l.LineWidth = 2;
    % end

    %Vertices, Edges
    [cx, cy, ~] = improfile(maskArtery, points_x, points_y);

    jj = 0;
    % Delete all points which are not in the maskArtery
    for ii=1:size(cx, 1)
        ry = round(cy(ii));
        rx = round(cx(ii));
        if (ry > 0 && ry <= size(v_RMS, 1) && rx > 0 && rx <= size(v_RMS, 2))
            jj = jj + 1;
            cy(jj) = ry;
            cx(jj) = rx;
        end
    end
    if (jj == 0) %If no points, no analysis.
        return;
    end

    signal = find_vessel_signal(maskArtery, v_RMS, cx, cy, jj);
    matrix(rr + 1, 1:length(signal)) = signal;
end


end