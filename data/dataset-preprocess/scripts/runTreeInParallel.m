tic
M = 12;
parfor (i = 72:77, M)
    treeGenerationParallel(i,500)
end
toc

% DONE:
% 10 trees
% i = 1:10, sample = 1
% 482.981059 seconds.
% 100 trees
% i = 11:20, sample = 10
% 3817.830117 seconds.
% 1000 trees
% i = 21:30, sample = 100
% Elapsed time is 36547.275090 seconds.
% maximum workers 6
% 1500 trees
% i = 31:45, sample = 100
% Elapsed time is 53052.976094 seconds.

% i = 46:71, sample = 100
% Elapsed time is 1287.060170 seconds.

% i = 72:77, sample = 500
% Elapsed time is 185.511309 seconds.