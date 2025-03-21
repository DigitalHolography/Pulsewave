function displaySplashScreen()
% Displays the splash screen during the app launch pausing for a short time
% and displaying git info
s = SplashScreen('Splashscreen', 'eyeflow_logo.png', 'ProgressBar', 'on', 'ProgressPosition', 5, 'ProgressRatio', 0.4);

t = timer;
t.StartDelay = 3;
t.TimerFcn = @(~, ~)delete(s);
start(t);

end
