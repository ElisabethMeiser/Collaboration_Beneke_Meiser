%% LeishiGuideRanking - Process TOMS CSV data and perform data ranking.
%
% This function loads TOMS CSV files, processes the data, and performs data ranking
% based on specified criteria. The ranked data can be saved as a CSV file.
%
% Usage:
%   LeishiGuideRanking(path_in, path_out, filename, rankByGeneID)
%   LeishiGuideRanking('/Users/eli/Desktop/TOM_Leish_Guide/', '/Users/eli/Desktop/TOM_Leish_Guide/TestSorting', 'TriTrypDB-59_LmexicanaMHOMGT2001U1103_Genome_LeishBASEedit_v1_output.primer.txt', true)

%
% Input:
%   - path_in: The path to the folder containing TOMS CSV data files.
%   - path_out: The path where the sorted data will be saved.
%   - filename: The name of the input CSV data file.
%   - rankByGeneID: A boolean flag to indicate whether to rank data by
%     GeneID. Not ranking is very quick, ranking is done with a loop and
%     takes some additional time.
%
% Output:
%   - The final sorted data table is saved as a .csv file in 'ResultFolder'
%     under the given 'filename'.
%
% Created by Elisabeth Meiser on 23rd September 2013.
%

function LeishiGuideRanking(path_in, path_out, filename, rankByGeneID)

ResultFolder = fullfile(path_out, 'GuideIDTable');

% Create the directory if it doesn't exist
if ~exist(ResultFolder, 'dir')
    mkdir(ResultFolder);
end


%% Step 1: Data Loading from a CSV file into a MATLAB table
GuideID = readtable(fullfile(path_in, filename),'PreserveVariableNames',true);

%% Step 2: Data Ranking with grouping by GeneID

% Initialize new columns to store results
GuideID.NumCytosine = zeros(size(GuideID.GuideSeq));
GuideID.NumPossibleStopCodons = zeros(size(GuideID.GuideSeq));
GuideID.RelativeProbability = zeros(size(GuideID.GuideSeq));
GuideID.C4to8Score = zeros(size(GuideID.GuideSeq));


% Step 1: Calculate NumCytosine in the 'Target' column using a 
% vectorized approach and 
% Step 2: determine the number of possible stop codons in the 
% 'GuideCoordinates' column.

for i = 1:numel(GuideID.GuideTarget)
    inputStringTarget = GuideID.GuideTarget{i};
    pattern = '\[(.*?)\]';
    matches = regexp(inputStringTarget, pattern, 'tokens');
    textWithinBrackets = [matches{:}];

    % Count 'C' occurrences in all textWithinBrackets
    numCsWithinBrackets = sum(cellfun(@(x) sum(x == 'C'), textWithinBrackets));
    GuideID.NumCytosine(i) = numCsWithinBrackets;

    % Count '->' occurrences in all textWithinBrackets to calculate NumPossibleStopCodons
    numStopCodons = sum(cellfun(@(x) sum(x == '>'), textWithinBrackets));
    GuideID.NumPossibleStopCodons(i) = numStopCodons;
end

% Step 3: Calculate the relative probability to generate a stop codon from any C.
GuideID.RelativeProbability = GuideID.NumPossibleStopCodons ./ GuideID.NumCytosine;

% Step 4: Check if 'C' is located in positions 4 to 8 also override toms
% column, which was not working
for i = 1:size(GuideID, 1)
    c_position = strfind(GuideID.GuideSeq{i}, 'C');
    if ~isempty(c_position) && any(c_position >= 4 & c_position <= 8)
        GuideID.C4to8Score(i) = 1;
    end
end

%Rename Toms Varible, because the nomenclature is trouble
if ismember('EditWindow4to8', GuideID.Properties.VariableNames)
    % Column 'EditWindow4to8' already exists, no need to rename it
    GuideID.EditWindow4to8 = GuideID.C4to8Score;
else
    % Rename 'EditWindow_4-8' to 'EditWindow4to8' if it doesn't exist
    GuideID = renamevars(GuideID, 'EditWindow_4-8', 'EditWindow4to8');
    GuideID.EditWindow4to8 = GuideID.C4to8Score;
end

% Step 4_b: Check if 'C' is located in positions 4 to 10 also override toms
% column, which was not working
for i = 1:size(GuideID, 1)
    c_position = strfind(GuideID.GuideSeq{i}, 'C');
    if ~isempty(c_position) && any(c_position >= 4 & c_position <= 10)
        GuideID.C4to10Score(i) = 1;
    end
end

%Rename Toms Varible, because the nomenclature is trouble
if ismember('EditWindow4to10', GuideID.Properties.VariableNames)
    % Column 'EditWindow4to10' already exists, no need to rename it
    GuideID.EditWindow4to10 = GuideID.C4to10Score;
else
    % Rename 'EditWindow_4-10' to 'EditWindow4to10' if it doesn't exist
    GuideID = renamevars(GuideID, 'EditWindow_4-10', 'EditWindow4to10');
    GuideID.EditWindow4to10 = GuideID.C4to10Score;
end

%% Calculate the final scores for all sequences in the table
scores = ((GuideID.RelativeProbability * 4) + GuideID.C4to8Score + GuideID.("TargetWithinFirst_20%_of_CDS")+ GuideID.("TargetWithinFirst_40%_of_CDS"));

%normalize scores by maximal possible score = 100
normalized_scores = scores*100/7; %num of contributions

% Add the normalized_scores as a new column 'C_Scores' to the table
GuideID.C_Scores = 100-normalized_scores; % 0 is the best guide, 100 is 
% the worst guide, counter intuitive, but this way Tom can
% sort the identidiers in a acending way...


%% Create a Guide ID called 'NewIdentifier' that contains all important
% scoring parmeters.


% Initialize a cell array to store the 'NewIdentifier' values
NewIdentifier = cell(size(GuideID, 1), 1);

GuideID.C_Scores_round = round(GuideID.C_Scores); % Only for naming

for i = 1:size(GuideID, 1)
    % Extract 'GeneID' value for the current row
    geneID = GuideID.GeneID{i};

    % Extract 'Strand' value for the current row
    strand = GuideID.Strand{i};

    % Extract 'C_Scores' value for the current row
    cScore = GuideID.C_Scores_round(i);

    % Extract guide coordinates
    guideCoordinatesValue = GuideID.GuideCoordinates{i};
    

    % Determine the strand part of the 'NewIdentifier'
    strandPart = '';
    if strcmp(strand, 'Plus')
        strandPart = '1';
    elseif strcmp(strand, 'Minus')
        strandPart = '0';
    end

    % Create the 'NewIdentifier' by concatenating the parts
    NewIdentifier{i} = [geneID '_'  num2str(cScore) '_' guideCoordinatesValue '_' strandPart];
end

% Add the 'NewIdentifier' to the 'GuideID' table
GuideID.NewIdentifier = NewIdentifier;


if rankByGeneID

    % Create a grouping variable based on 'GeneID'
    [uniqueGeneIDs, ~, groupID] = unique(GuideID.GeneID);

    % Initialize an empty table for the sorted result
    sortedGuideID = table();

    % Loop through each group defined by 'groupID' %Caution, this is slow, but
    % I havent found a good alternative yet

    for i = 1:max(groupID)
        groupIdx = groupID == i; % Find rows in the current group
        group = GuideID(groupIdx, :); % Extract the group
        sortedGroup = sortrows(group, 'C_Scores', 'descend'); % Sort the group by 'C_Scores' in descending order
        sortedGuideID = [sortedGuideID; sortedGroup]; % Append the sorted group to the result
    end



    %else
    % Step 2: Data Ranking without grouping
    %     sortedGuideID = customRankingFunction(GuideID);
    % Save the sortedGuideID as a CSV file

% Define the CSV file name (you can use strrep to replace the file extension)
csvFileName = fullfile(ResultFolder, strrep(filename, '.csv', '_sorted.csv'));

% Write the table to the CSV file
writetable(sortedGuideID, csvFileName);

else

% Save the GuideID as a CSV file
% Define the CSV file name (you can use strrep to replace the file extension)
csvFileName = fullfile(ResultFolder, strrep(filename, '.csv', '_unsort.csv'));

% Write the table to the CSV file
writetable(GuideID, csvFileName);

end

end
