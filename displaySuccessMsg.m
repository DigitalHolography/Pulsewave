function [success] = displaySuccessMsg(n)
            arguments
                n = [] % optional arg. empty array by default
            end
kaomoji = {'( ͡° ͜ʖ ͡°)', '⊙﹏⊙ ', '°‿‿°','(⊃｡•́‿•̀｡)⊃'};
if ~isempty(n)
    n = mod(n,size(kaomoji,2));
else % n is empty
    n = ceil(0.5+rand*(size(kaomoji,2)-1));
end
disp(' ');
disp('all done.');
disp(' ');
disp(['   ' kaomoji{n}]);
disp(' ');
disp('for robust rendering : ');
disp('1-flat-field correction, 2-background substraction');
disp('make disc mask : ');
disp('1- image field, to avoid bad segmentation, 2-fourier, to enable qty estimation of velocity');
disp(' ');
success = 1;
end

