function displaySplashScreen()
% Displays the splash screen during the app launch pausing for a short time
% and displaying git info
s = SplashScreen( 'Splashscreen', 'pulsewave_logo_temp.png','ProgressBar', 'on','ProgressPosition', 5,'ProgressRatio', 0.4 );
s.addText( 30, 50, 'Git info :', 'FontSize', 15, 'Color', [0 0 0.0] , 'Shadow','off')

gitBranchCommand = 'git symbolic-ref --short HEAD';
[statusBranch, resultBranch] = system(gitBranchCommand);

if statusBranch == 0
    resultBranch = strtrim(resultBranch);
    MessBranch = 'Current branch : ';
else
    MessBranch = 'Error getting current branch name. Git command output: ';
end

gitHashCommand = 'git rev-parse HEAD';
[statusHash, resultHash] = system(gitHashCommand);

if statusHash == 0 %hash command was successful
    resultHash = strtrim(resultHash);
    MessHash = 'Latest Commit Hash : ';
else
    MessHash = 'Error getting latest commit hash. Git command output: ';
end

% display git info
s.addText( 30, 80, [MessBranch, resultBranch], 'FontSize', 10, 'Color', [0 0 0] , 'Shadow','off')
s.addText( 30, 100, [MessHash, resultHash], 'FontSize', 10, 'Color', [0 0 0], 'Shadow','off' )

t = timer;
t.StartDelay = 3;
t.TimerFcn = @(~,~)delete(s);
start(t);

end