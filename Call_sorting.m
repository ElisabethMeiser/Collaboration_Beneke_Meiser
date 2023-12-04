% Specify the folder path
folderPath = '/Volumes/E_Meiser/LeishBASE';

% Other parameters for LeishiGuideRanking
path_in = '/Volumes/E_Meiser/LeishBASE';
path_out = '/Volumes/E_Meiser/LeishBASE/Sorting120423';
rankByGeneID = false ;

% List all files in the folder
files = dir(fullfile(folderPath, '*LeishBASEedit_v1_output.primer.txt')); % Replace 'your_file_extension' with the actual file extension

% Iterate over each file
for i = 1:length(files)
    % Get the current filename
    filename = files(i).name;
    
    % Construct the full path to the file
    fullPath = fullfile(folderPath, filename);
    
    % Call your LeishiGuideRanking function for each file
    LeishiGuideRanking(path_in, path_out, fullPath, rankByGeneID);
end
