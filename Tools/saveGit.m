function saveGit()
% SAVING GIT VERSION
% In the txt file in the folder : "log"

gitBranchCommand = 'git symbolic-ref --short HEAD';
[statusBranch, resultBranch] = system(gitBranchCommand);

if statusBranch == 0
    resultBranch = strtrim(resultBranch);
    MessBranch = 'Current branch : %s \r';
else

    vers = readlines('version.txt');
    MessBranch = ['PulseWave GitHub version ', char(vers)];
end

gitHashCommand = 'git rev-parse HEAD';
[statusHash, resultHash] = system(gitHashCommand);

if statusHash == 0 %hash command was successful
    resultHash = strtrim(resultHash);
    MessHash = 'Latest Commit Hash : %s \r';
else
    MessHash = '';
end

fprintf('==========================================\rGIT VERSION :\r');
fprintf(MessBranch, resultBranch);
fprintf(MessHash, resultHash);
fprintf('==========================================\r\n ');

end