# Collaboration_Tom_Beneke LeishiGuideRanking Function

## Description
The `LeishiGuideRanking` function is designed to process TOMS CSV files and perform data ranking based on specified criteria. It can be used to load TOMS CSV data, calculate various scoring parameters, and rank the data. The ranked data can be saved as a CSV file for further analysis.

## Features
- Loads TOMS CSV data files into a MATLAB table.
- Calculates various scoring parameters for the data.
- Ranks the data based on the calculated scores.
- Allows ranking by GeneID if desired.
- Saves the final sorted data as a CSV file for easy access.

## Usage
To use this function, follow these steps:

1. Provide the path to the folder containing TOMS CSV data files using `path_in`.
2. Specify the path where the sorted data will be saved using `path_out`.
3. Provide the name of the input CSV data file using `filename`.
4. Set the `rankByGeneID` flag to `true` if you want to rank data by GeneID.

Example:
```matlab
LeishiGuideRanking('/path/to/input_data', '/path/to/output_data', 'data_file.csv', true);
