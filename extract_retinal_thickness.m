%----------------------------------------------------------------------------
%NAME
%   extract_retinal_thickness.m
%PURPOSE
%   Calculate retinal layer thicknesses (RNFL, GCL, total receptor thickness) 
%   from specified locations within a foveal OCT image 
%DESCRPTION
%   This program extracts segmentation data from all retinal layers and
%   then calculates the RNFL, GCL, and total receptor thicknesses at 
%   specified locations within a foveal OCT image. 
%AUTHOR
%   Nikita Thomas
%----------------------------------------------------------------------------
%============================================================================
%               INITIALISE
%============================================================================
%               Clear MATLAB workspace and variables
clear; clc;
%               Select required xml file containing segmentation data for the
%               image with the deepest foveal pit (horizontal and vertical
%               will be separate xml files)
[surfacesfilename, surfacesfilepath] = uigetfile({'*.xml', 'XML files'}, 'Select surfaces file');
surfaces_data = xml2struct(surfacesfilename); 
%               Enter the degrees per pixel for the instrument
degreesperpixel = str2double(input('Please enter degrees per pixel > ','s'));
%               Enter the micrometres per pixel for the individual participant
micrometresperpixel = str2double(input('Please enter micrometres per pixel > ','s'));
%               Extract segmentation data from each individual layer
bscans_ILM = cell2mat(surfaces_data.surfaces.surface{1}.bscan(1));
bscans_RNFL = cell2mat(surfaces_data.surfaces.surface{2}.bscan(1));
bscans_GCL = cell2mat(surfaces_data.surfaces.surface{3}.bscan(1));
bscans_INL = cell2mat(surfaces_data.surfaces.surface{5}.bscan(1));
bscans_OPL = cell2mat(surfaces_data.surfaces.surface{6}.bscan(1));
bscans_BMEIS = cell2mat(surfaces_data.surfaces.surface{7}.bscan(1));
bscans_ISOSJ = cell2mat(surfaces_data.surfaces.surface{8}.bscan(1));
bscans_IB_OPR = cell2mat(surfaces_data.surfaces.surface{9}.bscan(1));
bscans_OB_OPR = cell2mat(surfaces_data.surfaces.surface{10}.bscan(1));
bscans_RPE = cell2mat(surfaces_data.surfaces.surface{11}.bscan(1));
bscans_OB_RPE = cell2mat(surfaces_data.surfaces.surface{12}.bscan(1));
%               Define locations/coordinates along 20 deg OCT scan for 
%               extraction of retinal layer thicknesses (in degrees)
OCT_coordinates = [3, 5, 7, 9, 11, 13, 15, 17];

%============================================================================
%               PROCESS SEGMENTATION DATA TO CALCULATE RETINAL THICKNESSES
%============================================================================
for i = 1:length(OCT_coordinates)
    pixel_coordinate = round(OCT_coordinates(i)/degreesperpixel);
    ILM_pixel = str2double(cell2mat(bscans_ILM.y(pixel_coordinate)).Text);
    RNFL_pixel = str2double(cell2mat(bscans_RNFL.y(pixel_coordinate)).Text);
    GCL_pixel = str2double(cell2mat(bscans_GCL.y(pixel_coordinate)).Text);
    INL_pixel = str2double(cell2mat(bscans_INL.y(pixel_coordinate)).Text);
    OPL_pixel = str2double(cell2mat(bscans_OPL.y(pixel_coordinate)).Text);
    BMEIS_pixel = str2double(cell2mat(bscans_BMEIS.y(pixel_coordinate)).Text);
    ISOSJ_pixel = str2double(cell2mat(bscans_ISOSJ.y(pixel_coordinate)).Text);
    IB_OPR_pixel = str2double(cell2mat(bscans_IB_OPR.y(pixel_coordinate)).Text);
    OB_OPR_pixel = str2double(cell2mat(bscans_OB_OPR.y(pixel_coordinate)).Text);
    RPE_pixel = str2double(cell2mat(bscans_RPE.y(pixel_coordinate)).Text);
    OB_RPE_pixel = str2double(cell2mat(bscans_OB_RPE.y(pixel_coordinate)).Text);
    
    RNFLthickness_pix(i) = RNFL_pixel - ILM_pixel;
    GCLthickness_pix(i) = GCL_pixel - RNFL_pixel;
    totalreceptorthickness_pix(i) = OB_RPE_pixel - INL_pixel;
end

%============================================================================
%               CONVERT TO MICROMETRES AND FINISH PROGRAM
%============================================================================
RNFLthickness_micrometres = RNFLthickness_pix*micrometresperpixel;
GCLthickness_micrometres = GCLthickness_pix*micrometresperpixel;
totalreceptorthickness_micrometres = totalreceptorthickness_pix*micrometresperpixel;

%               Print completion notification to the command window
disp('--------------------------------------------');
disp('Program complete');
