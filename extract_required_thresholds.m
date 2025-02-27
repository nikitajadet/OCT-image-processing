%----------------------------------------------------------------------------
%NAME
%   extract_required_thresholds.m
%PURPOSE
%   Extract thresholds from stimulus locations to be interpolated
%DESCRPTION
%   This program extracts thresholds (16 horizontal, 16 vertical)
%   required for interpolation at specified points within foveal OCT image.
%   The thresholds are extracted from a larger file containing thresholds
%   measured across the entire central visual field.
%AUTHOR
%   Nikita Thomas
%----------------------------------------------------------------------------
%============================================================================
%               INITIALISE
%============================================================================
%               Clear MATLAB workspace and variables
clear; clc;
%               Select required .csv file containing thresholds measured
%               across the entire central visual field and convert to matrix
[vf_filename, vf_filepath] = uigetfile({'*.csv', 'CSV files'}, 'Select ');
vf_file = fullfile(vf_filepath,vf_filename);
fid = fopen(vf_file,'r');
importCsv = textscan(fid,'%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q','Delimiter',',','Headerlines',1);
importCsv = [importCsv{:}];
vf_matrix = cellfun(@str2double,importCsv,'UniformOutput',false);
vf_matrix = cell2mat(vf_matrix);

%               Define x and y coordinates to search along horizontal axis
x_coordinates_horiz = [-7, -7, -5, -5, -3, -3, -1, -1, 1, 1, 3, 3, 5, 5, 7, 7]';
y_coordinates_horiz = [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1]';
%               Collate all horizontal coordinates together
horiz_coordinates = [x_coordinates_horiz y_coordinates_horiz];
%               Define x and y coordinates to search along vertical axis
x_coordinates_vert = [-1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1]';
y_coordinates_vert = [7, 7, 5, 5, 3, 3, 1, 1, -1, -1, -3, -3, -5, -5, -7, -7]';
%               Collate all vertical coordinates together
vert_coordinates = [x_coordinates_vert y_coordinates_vert];
%               Create a struct to house all visual field data from extracted
%               locations
Horizontal = struct('visualfield_data',[]);
Vertical = struct('visualfield_data',[]);

%============================================================================
%               EXTRACT VISUAL FIELD DATA FROM HORIZONTAL STIMULUS LOCATIONS
%============================================================================
for i = 1:length(horiz_coordinates)
    x = find(vf_matrix(:,1) == horiz_coordinates(i,1));
    y = find(vf_matrix(:,2) == horiz_coordinates(i,2));
    for j = 1:length(x)
        present(j) = ismember(x(j),y);
        if present(j) == true
            indexposition = x(j);
        else
            Horizontal(i).visualfield_data(1:17) = nan;
        end
    end
    %               Extract all visual field data for the extracted location
    Horizontal(i).visualfield_data(1:17) = vf_matrix(indexposition,1:17);
    %               Put threshold into seperate matrix (14th column
    %               of vf_matrix)
    horizontal_thresholds_all_locations(i) = Horizontal(i).visualfield_data(14);
end

%============================================================================
%               EXTRACT VISUAL FIELD DATA FROM HORIZONTAL STIMULUS LOCATIONS
%============================================================================
for i = 1:length(vert_coordinates)
    x = find(vf_matrix(:,1) == vert_coordinates(i,1));
    y = find(vf_matrix(:,2) == vert_coordinates(i,2));
    for j = 1:length(x)
        present(j) = ismember(x(j),y);
        if present(j) == true
            indexposition = x(j);
        else
            Vertical(i).visualfield_data(1:17) = nan;
        end
    end
    %               Extract all visual field data for the extracted location
    Vertical(i).visualfield_data(1:17) = vf_matrix(indexposition,1:17);
    %               Put threshold into seperate matrix (14th column
    %               of vf_matrix)
    vertical_thresholds_all_locations(i) = Vertical(i).visualfield_data(14);
end

%               Interpolate thresholds at specified OCT foveal locations
%               with 'interpolate_thresholds.m' function
[interpolated_horizontalthresholds, interpolated_verticalthresholds] = interpolate_thresholds(horizontal_thresholds_all_locations,vertical_thresholds_all_locations);

%============================================================================
%               FINISH PROGRAM
%============================================================================
%               Print completion notification to the command window
disp('--------------------------------------------');
disp('Program complete');
