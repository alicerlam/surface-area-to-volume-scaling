function [listOfFolderNames, listOfFileNames, topLevelFolder] = find_files(ext) 
%% Finds all files in a folder; adapted from https://matlab.fandom.com/wiki/FAQ#How_can_I_process_a_sequence_of_files.3F
clc;   
format long g;
format compact;

% Defines a starting folder.
start_path = fullfile(matlabroot, '\toolbox');
if ~exist(start_path, 'dir')
	start_path = matlabroot;
end

% Asks user to confirm the folder/change folder
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end
addpath(topLevelFolder);
fprintf('The top level folder is "%s".\n', topLevelFolder);

% Specifies the file pattern.
filePattern = sprintf(strcat('%s/**/*', ext), topLevelFolder);
allFileInfo = dir(filePattern);
isFolder = [allFileInfo.isdir]; % Logical list of what item is a folder or not.
allFileInfo(isFolder) = [];

% Gets a cell array of strings. 
listOfFolderNames = unique({allFileInfo.folder});
numberOfFolders = length(listOfFolderNames);
fprintf('The total number of folders to look in is %d.\n', numberOfFolders);
listOfFileNames = {allFileInfo.name};
totalNumberOfFiles = length(listOfFileNames);
fprintf('The total number of files in those %d folders is %d.\n', numberOfFolders, totalNumberOfFiles);

% Processes all files in those folders.
totalNumberOfFiles = length(allFileInfo);
if totalNumberOfFiles >= 1
	for k = 1 : totalNumberOfFiles
		thisFolder = allFileInfo(k).folder;
		thisBaseFileName = allFileInfo(k).name;
		fullFileName = fullfile(thisFolder, thisBaseFileName);
		[~, baseNameNoExt, ~] = fileparts(thisBaseFileName);
		fprintf('%s\n', baseNameNoExt);
	end
else
	fprintf('Folder %s has no files in it.\n', thisFolder);
end
fprintf('\nDone looking in all %d folders!\nFound %d files in the %d folders.\n', numberOfFolders, totalNumberOfFiles, numberOfFolders);
end 