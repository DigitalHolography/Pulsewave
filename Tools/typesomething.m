function typesomething(str, robot)
%typesomething
% function to type a string with a matlab robot
% made for azerty keyboards
% all charaecters are included for path : / \ (not so obvious)

if nargin < 2

    robot = java.awt.Robot;
end

str = char(str);
% Loop through each character in the string
for i = 1:length(str)
    charValue = str(i);

    % Convert character to the appropriate Java key code
    keyCode = getKeyCode(charValue);

    if isempty(keyCode)
        continue; % Skip non recognised characters
    end

    HoldShift = ischar(charValue) && (isstrprop(charValue, 'upper') || isstrprop(charValue, 'digit') || charValue == '.');

    if HoldShift
        robot.keyPress(java.awt.event.KeyEvent.VK_SHIFT); % Press Shift
    end

    HoldAltGr = ischar(charValue) && charValue == '\';

    if HoldAltGr
        robot.keyPress(17); robot.keyPress(18); % Press Ctrl AND Alt
    end

    robot.keyPress(keyCode);
    robot.keyRelease(keyCode);

    if HoldShift
        robot.keyRelease(java.awt.event.KeyEvent.VK_SHIFT); % Release Shift
    end

    if HoldAltGr
        robot.keyRelease(17); robot.keyRelease(18); % Press Ctrl AND Alt
    end

end

function keyCode = getKeyCode(character)
    % Map characters to key codes
    switch character
        case 'a', keyCode = java.awt.event.KeyEvent.VK_A;
        case 'b', keyCode = java.awt.event.KeyEvent.VK_B;
        case 'c', keyCode = java.awt.event.KeyEvent.VK_C;
        case 'd', keyCode = java.awt.event.KeyEvent.VK_D;
        case 'e', keyCode = java.awt.event.KeyEvent.VK_E;
        case 'f', keyCode = java.awt.event.KeyEvent.VK_F;
        case 'g', keyCode = java.awt.event.KeyEvent.VK_G;
        case 'h', keyCode = java.awt.event.KeyEvent.VK_H;
        case 'i', keyCode = java.awt.event.KeyEvent.VK_I;
        case 'j', keyCode = java.awt.event.KeyEvent.VK_J;
        case 'k', keyCode = java.awt.event.KeyEvent.VK_K;
        case 'l', keyCode = java.awt.event.KeyEvent.VK_L;
        case 'm', keyCode = java.awt.event.KeyEvent.VK_M;
        case 'n', keyCode = java.awt.event.KeyEvent.VK_N;
        case 'o', keyCode = java.awt.event.KeyEvent.VK_O;
        case 'p', keyCode = java.awt.event.KeyEvent.VK_P;
        case 'q', keyCode = java.awt.event.KeyEvent.VK_Q;
        case 'r', keyCode = java.awt.event.KeyEvent.VK_R;
        case 's', keyCode = java.awt.event.KeyEvent.VK_S;
        case 't', keyCode = java.awt.event.KeyEvent.VK_T;
        case 'u', keyCode = java.awt.event.KeyEvent.VK_U;
        case 'v', keyCode = java.awt.event.KeyEvent.VK_V;
        case 'w', keyCode = java.awt.event.KeyEvent.VK_W;
        case 'x', keyCode = java.awt.event.KeyEvent.VK_X;
        case 'y', keyCode = java.awt.event.KeyEvent.VK_Y;
        case 'z', keyCode = java.awt.event.KeyEvent.VK_Z;

        case 'A', keyCode = java.awt.event.KeyEvent.VK_A;
        case 'B', keyCode = java.awt.event.KeyEvent.VK_B;
        case 'C', keyCode = java.awt.event.KeyEvent.VK_C;
        case 'D', keyCode = java.awt.event.KeyEvent.VK_D;
        case 'E', keyCode = java.awt.event.KeyEvent.VK_E;
        case 'F', keyCode = java.awt.event.KeyEvent.VK_F;
        case 'G', keyCode = java.awt.event.KeyEvent.VK_G;
        case 'H', keyCode = java.awt.event.KeyEvent.VK_H;
        case 'I', keyCode = java.awt.event.KeyEvent.VK_I;
        case 'J', keyCode = java.awt.event.KeyEvent.VK_J;
        case 'K', keyCode = java.awt.event.KeyEvent.VK_K;
        case 'L', keyCode = java.awt.event.KeyEvent.VK_L;
        case 'M', keyCode = java.awt.event.KeyEvent.VK_M;
        case 'N', keyCode = java.awt.event.KeyEvent.VK_N;
        case 'O', keyCode = java.awt.event.KeyEvent.VK_O;
        case 'P', keyCode = java.awt.event.KeyEvent.VK_P;
        case 'Q', keyCode = java.awt.event.KeyEvent.VK_Q;
        case 'R', keyCode = java.awt.event.KeyEvent.VK_R;
        case 'S', keyCode = java.awt.event.KeyEvent.VK_S;
        case 'T', keyCode = java.awt.event.KeyEvent.VK_T;
        case 'U', keyCode = java.awt.event.KeyEvent.VK_U;
        case 'V', keyCode = java.awt.event.KeyEvent.VK_V;
        case 'W', keyCode = java.awt.event.KeyEvent.VK_W;
        case 'X', keyCode = java.awt.event.KeyEvent.VK_X;
        case 'Y', keyCode = java.awt.event.KeyEvent.VK_Y;
        case 'Z', keyCode = java.awt.event.KeyEvent.VK_Z;

        case '0', keyCode = java.awt.event.KeyEvent.VK_0;
        case '1', keyCode = java.awt.event.KeyEvent.VK_1;
        case '2', keyCode = java.awt.event.KeyEvent.VK_2;
        case '3', keyCode = java.awt.event.KeyEvent.VK_3;
        case '4', keyCode = java.awt.event.KeyEvent.VK_4;
        case '5', keyCode = java.awt.event.KeyEvent.VK_5;
        case '6', keyCode = java.awt.event.KeyEvent.VK_6;
        case '7', keyCode = java.awt.event.KeyEvent.VK_7;
        case '8', keyCode = java.awt.event.KeyEvent.VK_8;
        case '9', keyCode = java.awt.event.KeyEvent.VK_9;

        case ' ', keyCode = java.awt.event.KeyEvent.VK_SPACE;

        case '\', keyCode = java.awt.event.KeyEvent.VK_8;

        case '/', keyCode = 111;

        case ':', keyCode = java.awt.event.KeyEvent.VK_COLON;

        case '_', keyCode = java.awt.event.KeyEvent.VK_8;

        case '.', keyCode = java.awt.event.KeyEvent.VK_SEMICOLON;

        case ']', keyCode = java.awt.event.KeyEvent.VK_SEMICOLON;

        otherwise , keyCode = []; % Skip characters that don't have keycodes
    end

end

end
