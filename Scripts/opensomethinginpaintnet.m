robot = java.awt.Robot;

INPUT_PATH = "X:\250117\AUZ_L\250117_AUZ0752_L_HD_0\pulsewave\mask\M0.png";
INPUT_PATH_Layer = "X:\250220\250220_BOM0753_L_1_HD_1\pulsewave\mask\DiaSysRGB.png";

st = system(['start "" "', 'C:\Program Files\paint.net\paintdotnet.exe', '"']); % launches pdn in non blocking mode
robot.delay(1500); % wait 1.5s for paint net to open
displaysize = get(0, 'ScreenSize'); % for a script indifferent of the display size 3=width 4=height

sratio = displaysize(3:4) ./ [1920 1200];

set(0,'PointerLocation',[20 1165] .* sratio); % Go on "File"
robot.mousePress  (java.awt.event.InputEvent.BUTTON1_MASK); % Click
robot.mouseRelease(java.awt.event.InputEvent.BUTTON1_MASK);

set(0,'PointerLocation',[20 1115] .* sratio);
robot.mousePress  (java.awt.event.InputEvent.BUTTON1_MASK); % Click "Open"
robot.mouseRelease(java.awt.event.InputEvent.BUTTON1_MASK);


pause(1);

typesomething(INPUT_PATH);
set(0,'PointerLocation',[1088 464] .* sratio);
robot.mousePress(java.awt.event.InputEvent.BUTTON1_MASK); % Click "Open"
robot.mouseRelease(java.awt.event.InputEvent.BUTTON1_MASK);

pause(1);
sz = size(imread(INPUT_PATH));

if sz(1)~=sz(2)
    robot.keyPress(17); robot.keyPress(java.awt.event.KeyEvent.VK_R); % ctrlR resize
    robot.keyRelease(17); robot.keyRelease(java.awt.event.KeyEvent.VK_R); % ctrlR resize
    typesomething(string(num2str(max(sz(1),sz(2)))));
    pause(0.05);
    robot.keyPress(java.awt.event.KeyEvent.VK_TAB);
    robot.keyRelease(java.awt.event.KeyEvent.VK_TAB);
    pause(0.05);
    typesomething(string(num2str(max(sz(1),sz(2)))));
    robot.keyPress(java.awt.event.KeyEvent.VK_ENTER);
    robot.keyRelease(java.awt.event.KeyEvent.VK_ENTER);
end

set(0,'PointerLocation',[187	1164] .* sratio);
robot.mousePress(java.awt.event.InputEvent.BUTTON1_MASK); % Click "Layers"
robot.mouseRelease(java.awt.event.InputEvent.BUTTON1_MASK);

pause(1);
sz1 = size(imread(INPUT_PATH_Layer));
if sz1(1:2) ~=sz(1:2)
    img = imresize(imread(INPUT_PATH_Layer),[max(sz(1),sz(2)) max(sz(1),sz(2))]);
    imwrite(img, INPUT_PATH_Layer, 'png');
end

set(0,'PointerLocation',[196	1013] .* sratio); 
robot.mousePress(java.awt.event.InputEvent.BUTTON1_MASK); % Click "Import from file"
robot.mouseRelease(java.awt.event.InputEvent.BUTTON1_MASK);

pause(0.5);
typesomething(INPUT_PATH_Layer);
robot.keyPress(java.awt.event.KeyEvent.VK_ENTER);
robot.keyRelease(java.awt.event.KeyEvent.VK_ENTER);

pause(0.5);
robot.keyPress(java.awt.event.KeyEvent.VK_F4);
robot.keyRelease(java.awt.event.KeyEvent.VK_F4);

pause(0.5);
robot.keyPress(java.awt.event.KeyEvent.VK_TAB);
robot.keyRelease(java.awt.event.KeyEvent.VK_TAB);

robot.keyPress(java.awt.event.KeyEvent.VK_TAB);
robot.keyRelease(java.awt.event.KeyEvent.VK_TAB);

typesomething("55"); % set opacity to 55

robot.keyPress(java.awt.event.KeyEvent.VK_ENTER);
robot.keyRelease(java.awt.event.KeyEvent.VK_ENTER);

robot.keyPress(17);robot.keyPress(java.awt.event.KeyEvent.VK_SHIFT);robot.keyPress(java.awt.event.KeyEvent.VK_N);
robot.keyRelease(17);robot.keyRelease(java.awt.event.KeyEvent.VK_SHIFT);robot.keyRelease(java.awt.event.KeyEvent.VK_N);

set(0,'PointerLocation',[955	593] .* sratio); 
robot.keyPress(java.awt.event.KeyEvent.VK_X); % inverse color default
robot.keyRelease(java.awt.event.KeyEvent.VK_X);
robot.keyPress(java.awt.event.KeyEvent.VK_F); % select fill
robot.keyRelease(java.awt.event.KeyEvent.VK_F);
robot.mousePress(java.awt.event.InputEvent.BUTTON1_MASK); % Click to fill layer
robot.mouseRelease(java.awt.event.InputEvent.BUTTON1_MASK);

pause(0.5);
robot.keyPress(java.awt.event.KeyEvent.VK_F4);
robot.keyRelease(java.awt.event.KeyEvent.VK_F4);

pause(0.5);
robot.keyPress(java.awt.event.KeyEvent.VK_TAB);
robot.keyRelease(java.awt.event.KeyEvent.VK_TAB);

robot.keyPress(java.awt.event.KeyEvent.VK_TAB);
robot.keyRelease(java.awt.event.KeyEvent.VK_TAB);

typesomething("55"); % set opacity to 55

robot.keyPress(java.awt.event.KeyEvent.VK_ENTER);
robot.keyRelease(java.awt.event.KeyEvent.VK_ENTER);

robot.keyPress(java.awt.event.KeyEvent.VK_B); % select pen
robot.keyRelease(java.awt.event.KeyEvent.VK_B);




% save a pdn file
robot.keyPress(17);robot.keyPress(java.awt.event.KeyEvent.VK_S);
robot.keyRelease(17);robot.keyRelease(java.awt.event.KeyEvent.VK_S);

robot.keyPress(java.awt.event.KeyEvent.VK_ENTER);
robot.keyRelease(java.awt.event.KeyEvent.VK_ENTER);
