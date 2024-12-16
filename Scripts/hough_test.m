% Read the image
A = imread("X:\241216\Vienna\HoloDoppler_Vienna_1\241209_capture_HD_0\png\241209_capture_HD_0_M0.png");

% Detect circles in the image
[centers, radii] = imfindcircles(A, [40, 100], 'ObjectPolarity', 'dark');

% Display the image
imshow(A);
hold on;

% Plot the detected circles in red
viscircles(centers, radii, 'Color', 'r');

% Release the hold on the current figure
hold off;
