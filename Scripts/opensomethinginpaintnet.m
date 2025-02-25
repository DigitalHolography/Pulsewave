robot = java.awt.Robot;

st = system(['start "" "', 'C:\Program Files\paint.net\paintdotnet.exe', '"']); % launches pdn in non blocking mode
robot.delay(1500); % wait 1.5s for paint net to open
displaysize = get(0, 'ScreenSize'); % for a script indifferent of the display size 3=width 4=height

sratio = displaysize(3:4) ./ [1920 1200];

set(0,'PointerLocation',[20 1165] .* sratio); % Go on "File"
robot.mousePress  (java.awt.event.InputEvent.BUTTON1_MASK);
robot.mouseRelease(java.awt.event.InputEvent.BUTTON1_MASK);




