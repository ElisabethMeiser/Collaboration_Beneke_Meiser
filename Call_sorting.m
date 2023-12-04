%LeishiGuideRanking(path_in, path_out, filename, rankByGeneID)

%LeishiGuideRanking('/Users/eli/Desktop/TOM_Leish_Guide/', '/Users/eli/Desktop/TOM_Leish_Guide/TestSorting', 'TriTrypDB-59_LmexicanaMHOMGT2001U1103_Genome_LeishBASEedit_v1_output.primer.txt', false)

% Specify the folder path
folderPath = '/Users/eli/Desktop/TOM_Leish_Guide/GuideIDTable';

% Other parameters for LeishiGuideRanking
pathIn = '/Users/eli/Desktop/TOM_Leish_Guide/';
pathOut = '/Users/eli/Desktop/TOM_Leish_Guide/Sorting120423';
rankByGeneID = false ;

% List all files in the folder
files = dir(fullfile(folderPath, '*.txt')); % Replace 'your_file_extension' with the actual file extension

% Iterate over each file
for i = 1:length(files)
    % Get the current filename
    currentFile = files(i).name;
    
    % Construct the full path to the file
    fullPath = fullfile(folderPath, currentFile);
    
    % Call your LeishiGuideRanking function for each file
    LeishiGuideRanking(pathIn, pathOut, fullPath, rankByGeneID);
end
