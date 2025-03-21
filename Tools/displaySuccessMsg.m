function [success] = displaySuccessMsg(n)

arguments
    n = [] % optional arg. empty array by default
end

kaomoji = {'( ͡° ͜ʖ ͡°)', '⊙﹏⊙ ', '°‿‿°', '(⊃｡•́‿•̀｡)⊃'};

if ~isempty(n)
    n = mod(n, size(kaomoji, 2));
else % n is empty
    n = ceil(0.5 + rand * (size(kaomoji, 2) - 1));
end

disp(' ');
disp('all done.');
disp(' ');
disp(['   ' kaomoji{n}]);
disp(' ');
success = 1;
end
