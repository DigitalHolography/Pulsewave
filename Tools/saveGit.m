function saveGit()
% SAVING GIT VERSION
% In the txt file in the folder : "log"

% Get the current branch name
gitBranchCommand = 'git symbolic-ref --short HEAD';
[statusBranch, resultBranch] = system(gitBranchCommand);

if statusBranch == 0
    resultBranch = strtrim(resultBranch);
    MessBranch = 'Current branch : %s \r';
else
    vers = readlines('version.txt');
    MessBranch = sprintf('PulseWave GitHub version %s\r', char(vers));
end

% Get the latest commit hash
gitHashCommand = 'git rev-parse HEAD';
[statusHash, resultHash] = system(gitHashCommand);

if statusHash == 0 %hash command was successful
    resultHash = strtrim(resultHash);
    MessHash = 'Latest Commit Hash : %s \r';
else
    MessHash = '';
end

% Get the most recent tag
gitTagCommand = 'git describe --tags';
[statusTag, resultTag] = system(gitTagCommand);

if statusTag == 0 %tag command was successful
    resultTag = strtrim(resultTag);
    MessTag = 'Most recent tag : %s \r';
else
    MessTag = '';
end

% Print the results
fprintf('==========================================\rGIT VERSION :\r');
fprintf(MessBranch, resultBranch);
fprintf(MessHash, resultHash);
fprintf(MessTag, resultTag);
fprintf('==========================================\r');

end
