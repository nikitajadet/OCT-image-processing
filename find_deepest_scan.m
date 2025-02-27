%----------------------------------------------------------------------------
%NAME
%   find_deepest_scan.m
%PURPOSE
%   Find the B-scan containing the foveal image with the smallest foveal pit 
%   height (deepest foveal pit)
%DESCRPTION
%   This program isolates 25-pixel windows across the ILM within a B-scan
%   and draws a line of best fit for each window. Following this, a line
%   normal to the line  of best fit is drawn to the photoreceptor layer.
%   The program then finds the smallest of these normal lines and
%   identifies this as the foveal pit height. The program then looks
%   through all B-scans to find the image with the smallest foveal pit
%   height overall.
%AUTHOR
%   Nikita Thomas
%----------------------------------------------------------------------------
%============================================================================
%               INITIALISE
%============================================================================
%               Clear MATLAB workspace and variables
clear; clc;
%               Select required xml file that contains segmentation data
%               (for all B-scans segmented)
[surfacesfilename, surfacesfilepath] = uigetfile({'*.xml', 'XML files'}, 'Select surfaces file');
%               Read in xml file 
surfaces_data = xml2struct(surfacesfilename); 
%               Isolate ILM and IS/OS segmented data and convert into matrix
bscans_ILM = cell2mat(surfaces_data.surfaces.surface{1}.bscan(1:end));
bscans_ISOS = cell2mat(surfaces_data.surfaces.surface{11}.bscan(1:end));
%               Create a struct to house scan data from each individual B-scan
pre_allocate_scanproperties = 1:length(bscans_ILM);
ScanProperties(pre_allocate_scanproperties) = struct('bscans1',[],'bscans11',[],'ILM',[],'ISOS',[],'yILM',[],'xILM',[],'ILM_InterX_format',[],'min_distance',[]); 

%============================================================================
%               FILL STRUCT
%============================================================================
%               Loop to house segmented data into created struct fields for
%               each individual B-scan
for j = 1:length(bscans_ILM)
    ScanProperties(j).bscans1 = cell2mat(bscans_ILM(j).y);
    ScanProperties(j).ILM = nan(length(ScanProperties(j).bscans1),1);
    for ii = 1:length(ScanProperties(j).bscans1)
        ScanProperties(j).ILM(ii) = str2double(ScanProperties(j).bscans1(ii).Text);
    end
    
    ScanProperties(j).bscans11 = cell2mat(bscans_ISOS(j).y);
    ScanProperties(j).ISOS = nan(length(ScanProperties(j).bscans11),1);
    for ii = 1:length(ScanProperties(j).bscans11)
        ScanProperties(j).ISOS(ii) = str2double(ScanProperties(j).bscans11(ii).Text);
    end
end

%============================================================================
%               PROCESS DATA TO FIND FOVEAL PIT HEIGHT
%============================================================================
%               Create a new struct for processing data and pre-allocate to 
%               fill fields with 25-pixel windows
pre_allocate_scanparameters = 1:length(ScanProperties(j).xILM)-24;
ScanParameters(pre_allocate_scanparameters) = struct('pixelwindow_ILM',[],'x_coordinates_ILM',[],'ILM_coefficients',[],'normalline_calculation',[],'m1_to_extend_normal',[],'m2_to_extend_normal',[],'extend_normalline_x',[],'extend_normalline_y',[],'find_300_points_x',[],'extended_normalline',[],'pixel_InterX_format',[],'ISOSpixel_normal_crosses',[],'distance_ILM_to_ISOS',[]);
%               Loop to isolate each pixel window, obtain the coefficients of
%               the line of best fit, plot normal line to ISOS layer, and
%               find smallest normal line (foveal pit height)
for j = 1:length(bscans_ILM)
    for i = 1:length(ScanProperties(j).xILM)-24
        ScanProperties(j).ScanParameters(i).pixelwindow_ILM = ScanProperties(j).yILM(i:i+24)'; 
        ScanProperties(j).ScanParameters(i).x_coordinates_ISOS = ScanProperties(j).xILM(i:i+24); 
        ScanProperties(j).ScanParameters(i).ILM_coefficients = polyfit(ScanProperties(j).ScanParameters(i).x_coordinates_ILM,ScanProperties(j).ScanParameters(i).pixelwindow_ILM,1); 
        ScanProperties(j).ScanParameters(i).normalline_calculation = -1/ScanProperties(j).ScanParameters(i).ILM_coefficients(1); 
        %       Extend normal line in order to reach ISOS layer
        ScanProperties(j).ScanParameters(i).m1_to_extend_normal = mean(ScanProperties(j).ScanParameters(i).x_coordinates_ILM)+0.5 - mean(ScanProperties(j).ScanParameters(i).x_coordinates_ILM); 
        ScanProperties(j).ScanParameters(i).m2_to_extend_normal = mean(ScanProperties(j).ScanParameters(i).pixelwindow_ILM)+0.5*ScanProperties(j).ScanParameters(i).normalline_calculation - mean(ScanProperties(j).ScanParameters(i).pixelwindow_ILM); 
        ScanProperties(j).ScanParameters(i).extend_normalline_x = [mean(ScanProperties(j).ScanParameters(i).x_coordinates_ILM) mean(ScanProperties(j).ScanParameters(i).x_coordinates_ILM)+0.5]+300*[-ScanProperties(j).ScanParameters(i).m1_to_extend_normal ScanProperties(j).ScanParameters(i).m1_to_extend_normal]; ScanProperties(j).ScanParameters(i).extend_normalline_y = [mean(ScanProperties(j).ScanParameters(i).pixelwindow_ILM) mean(ScanProperties(j).ScanParameters(i).pixelwindow_ILM)+0.5*ScanProperties(j).ScanParameters(i).normalline_calculation]+300*[-ScanProperties(j).ScanParameters(i).m2_to_extend_normal ScanProperties(j).ScanParameters(i).m2_to_extend_normal];
        ScanProperties(j).ScanParameters(i).find_300_points_x = ScanProperties(j).ScanParameters(i).extend_normalline_x(1):ScanProperties(j).ScanParameters(i).extend_normalline_x(2);
        ScanProperties(j).ScanParameters(i).extended_normalline = interp1(ScanProperties(j).ScanParameters(i).extend_normalline_x,ScanProperties(j).ScanParameters(i).extend_normalline_y,ScanProperties(j).ScanParameters(i).find_300_points_x);
        %       Group x and y normal coordinates into a compatible matrix 
        %       format to be later used with function 'InterX'
        ScanProperties(j).ScanParameters(i).pixel_InterX_format = [ScanProperties(j).ScanParameters(i).find_300_points_x;ScanProperties(j).ScanParameters(i).extended_normalline]; 
        %       Find the ILM pixel that the normal line crosses 
        %       (InterX format = x coordinate;y coordinate)
        ScanProperties(j).ScanParameters(i).ISOSpixel_normal_crosses = InterX(ScanProperties(j).ScanParameters(i).pixel_InterX_format,ScanProperties(j).ILM_InterX_format);
        if ~any(isempty(ScanProperties(j).ScanParameters(i).ILMpixel_normal_crosses))
            ScanProperties(j).ScanParameters(i).distance_ILM_to_ISOS = median(ScanProperties(j).ScanParameters(i).pixelwindow_ILM) - ScanProperties(j).ScanParameters(i).ILMpixel_normal_crosses(2);
        else
            ScanProperties(j).ScanParameters(i).distance_ILM_to_ISOS = [];
        end
    end
    %           Find the B-scan with the smallest foveal pit height overall
    ScanProperties(j).min_distance = min([ScanProperties(j).ScanParameters.distance_ILM_to_ISOS]);
    [min_distance_allscans, min_scanNumber] = min([ScanProperties.min_distance]);
end

%============================================================================
%               PLOT DATA TO VISUALISE NORMAL LINES FOR IMAGE WITH SMALLEST
%               FOVEAL PIT HEIGHT
%============================================================================
%               Plot data from image with smallest foveal pit height
plot(ScanProperties(min_scanNumber).ILM);
plot(ScanProperties(min_scanNumber).ISOS);
for i = 1:length(ScanProperties(min_scanNumber).xISOS)-24
    plot([median(ScanProperties(min_scanNumber).ScanParameters(i).x_coordinates_ISOS) median(ScanProperties(min_scanNumber).ScanParameters(i).x_coordinates_ISOS)+0.5]+300*[-ScanProperties(min_scanNumber).ScanParameters(i).m1_to_extend_normal ScanProperties(min_scanNumber).ScanParameters(i).m1_to_extend_normal], [median(ScanProperties(min_scanNumber).ScanParameters(i).pixelwindow_ISOS) median(ScanProperties(min_scanNumber).ScanParameters(i).pixelwindow_ISOS)+0.5*ScanProperties(min_scanNumber).ScanParameters(i).normalline_calculation]+300*[-ScanProperties(min_scanNumber).ScanParameters(i).m2_to_extend_normal ScanProperties(min_scanNumber).ScanParameters(i).m2_to_extend_normal]);
end

%============================================================================
%               DISPLAY RESULTS
%============================================================================
%              Print the results 
disp(['Smallest foveal pit height out of all B-scan images: ', num2str(min_distance_allscans)])
disp(['Image with the smallest foveal pit height: no.',num2str(min_scanNumber)])
%              Print completion notification to the command window
disp('--------------------------------------------');
disp('Program complete');
