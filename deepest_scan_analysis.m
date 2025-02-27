[surfacesFileName, surfacesFilePath] = uigetfile({'*.xml', 'XML files'}, 'Select surfaces file');
surfaces_data = xml2struct(surfacesFileName); %read deepest xml file into script

pixel_cutting = true; % set this to true if the start at end of the image doesn't give us any useful information
max_pixel_cutoff = 900;
min_pixel_cutoff = 100;
horizontal = true;
eye = questdlg('Which eye?','Eye?','Right','Left','Right');

participantId = input('Please enter participant Id > ','s');

Bscans_ILM = cell2mat(surfaces_data.surfaces.surface{1}.bscan(1:end));
Bscans_RNFL = cell2mat(surfaces_data.surfaces.surface{2}.bscan(1:end));
Bscans_GCL = cell2mat(surfaces_data.surfaces.surface{3}.bscan(1:end));
Bscans_IPL = cell2mat(surfaces_data.surfaces.surface{4}.bscan(1:end));
Bscans_INL = cell2mat(surfaces_data.surfaces.surface{5}.bscan(1:end));
Bscans_OPL = cell2mat(surfaces_data.surfaces.surface{6}.bscan(1:end));
Bscans_BMEIS = cell2mat(surfaces_data.surfaces.surface{7}.bscan(1:end));
Bscans_ISOSJ = cell2mat(surfaces_data.surfaces.surface{8}.bscan(1:end));
Bscans_IB_OPR = cell2mat(surfaces_data.surfaces.surface{9}.bscan(1:end));
%Bscans_OB_OPR = cell2mat(surfaces_data.surfaces.surface{10}.bscan(1:end));
Bscans_RPE = cell2mat(surfaces_data.surfaces.surface{10}.bscan(1:end));
Bscans_OB_RPE = cell2mat(surfaces_data.surfaces.surface{11}.bscan(1:end));

pre_allocate_ScanProperties = 1:length(Bscans_ILM);
ScanProperties(pre_allocate_ScanProperties) = struct('Bscans1',[],'Bscans2',[],'Bscans3',[],'Bscans4',[],'Bscans5',[],'Bscans6',[],'Bscans7',[],'Bscans8',[],'Bscans9',[],'Bscans10',[],'Bscans11',[],'Bscans12',[],'ILM',[],'RNFL',[],'GCL',[],'IPL',[],'INL',[],'OPL',[],'BMEIS',[],'ISOSJ',[],'IB_OPR',[],'OB_OPR',[],'RPE',[],'OB_RPE',[]);

for j = 1:length(Bscans_ILM)
    ScanProperties(j).Bscans1 = cell2mat(Bscans_ILM(j).y);
    ScanProperties(j).ILM = nan(length(ScanProperties(j).Bscans1),1);
    for ii = 1:length(ScanProperties(j).Bscans1)
        ScanProperties(j).ILM(ii) = str2double(ScanProperties(j).Bscans1(ii).Text); % convert ILM pixel data into struct.
    if pixel_cutting
        ScanProperties(j).ILM_pixelcut = ScanProperties(j).ILM(min_pixel_cutoff:max_pixel_cutoff);
    end
    end
    
    ScanProperties(j).Bscans2 = cell2mat(Bscans_RNFL(j).y);
    ScanProperties(j).RNFL = nan(length(ScanProperties(j).Bscans2),1);
    for ii = 1:length(ScanProperties(j).Bscans2)
        ScanProperties(j).RNFL(ii) = str2double(ScanProperties(j).Bscans2(ii).Text); % convert RNFL xml pixel data into struct
    if pixel_cutting
        ScanProperties(j).RNFL_pixelcut = ScanProperties(j).RNFL(min_pixel_cutoff:max_pixel_cutoff);
    end
    end
    
    ScanProperties(j).Bscans3 = cell2mat(Bscans_GCL(j).y);
    ScanProperties(j).GCL = nan(length(ScanProperties(j).Bscans3),1);
    for ii = 1:length(ScanProperties(j).Bscans3)
        ScanProperties(j).GCL(ii) = str2double(ScanProperties(j).Bscans3(ii).Text); % convert GCL xml pixel data into struct
    if pixel_cutting
        ScanProperties(j).GCL_pixelcut = ScanProperties(j).GCL(min_pixel_cutoff:max_pixel_cutoff);
    end
    end
    
    ScanProperties(j).Bscans4 = cell2mat(Bscans_IPL(j).y);
    ScanProperties(j).IPL = nan(length(ScanProperties(j).Bscans4),1);
    for ii = 1:length(ScanProperties(j).Bscans4)
        ScanProperties(j).IPL(ii) = str2double(ScanProperties(j).Bscans4(ii).Text); % convert IPL xml pixel data into struct
    if pixel_cutting
        ScanProperties(j).IPL_pixelcut = ScanProperties(j).IPL(min_pixel_cutoff:max_pixel_cutoff);
    end
    end
    
    ScanProperties(j).Bscans5 = cell2mat(Bscans_INL(j).y);
    ScanProperties(j).INL = nan(length(ScanProperties(j).Bscans5),1);
    for ii = 1:length(ScanProperties(j).Bscans5)
        ScanProperties(j).INL(ii) = str2double(ScanProperties(j).Bscans5(ii).Text); % convert INL xml pixel data into struct
    if pixel_cutting
        ScanProperties(j).INL_pixelcut = ScanProperties(j).INL(min_pixel_cutoff:max_pixel_cutoff);
    end
    end
    
    ScanProperties(j).Bscans6 = cell2mat(Bscans_OPL(j).y);
    ScanProperties(j).OPL = nan(length(ScanProperties(j).Bscans6),1);
    for ii = 1:length(ScanProperties(j).Bscans6)
        ScanProperties(j).OPL(ii) = str2double(ScanProperties(j).Bscans6(ii).Text); % convert OPL xml pixel data into struct
    if pixel_cutting
        ScanProperties(j).OPL_pixelcut = ScanProperties(j).OPL(min_pixel_cutoff:max_pixel_cutoff);
    end
    end
    
    ScanProperties(j).Bscans7 = cell2mat(Bscans_BMEIS(j).y);
    ScanProperties(j).BMEIS = nan(length(ScanProperties(j).Bscans7),1);
    for ii = 1:length(ScanProperties(j).Bscans7)
        ScanProperties(j).BMEIS(ii) = str2double(ScanProperties(j).Bscans7(ii).Text); % convert BMEIS xml pixel data into struct
    if pixel_cutting
        ScanProperties(j).BMEIS_pixelcut = ScanProperties(j).BMEIS(min_pixel_cutoff:max_pixel_cutoff);
    end
    end
    
    ScanProperties(j).Bscans8 = cell2mat(Bscans_ISOSJ(j).y);
    ScanProperties(j).ISOSJ = nan(length(ScanProperties(j).Bscans8),1);
    for ii = 1:length(ScanProperties(j).Bscans8)
        ScanProperties(j).ISOSJ(ii) = str2double(ScanProperties(j).Bscans8(ii).Text); % convert ISOSJ xml pixel data into struct
    if pixel_cutting
        ScanProperties(j).ISOSJ_pixelcut = ScanProperties(j).ISOSJ(min_pixel_cutoff:max_pixel_cutoff);
    end
    end
    
    ScanProperties(j).Bscans9 = cell2mat(Bscans_IB_OPR(j).y);
    ScanProperties(j).IB_OPR = nan(length(ScanProperties(j).Bscans9),1);
    for ii = 1:length(ScanProperties(j).Bscans9)
        ScanProperties(j).IB_OPR(ii) = str2double(ScanProperties(j).Bscans9(ii).Text); % convert IB_OPR xml pixel data into struct
    if pixel_cutting
        ScanProperties(j).IB_OPR_pixelcut = ScanProperties(j).IB_OPR(min_pixel_cutoff:max_pixel_cutoff);
    end
    end
    
%     ScanProperties(j).Bscans10 = cell2mat(Bscans_OB_OPR(j).y);
%     ScanProperties(j).OB_OPR = nan(length(ScanProperties(j).Bscans10),1);
%     for ii = 1:length(ScanProperties(j).Bscans10)
%         ScanProperties(j).OB_OPR(ii) = str2double(ScanProperties(j).Bscans10(ii).Text); % convert OB_OPR xml pixel data into struct
%     if pixel_cutting
%         ScanProperties(j).OB_OPR_pixelcut = ScanProperties(j).OB_OPR(min_pixel_cutoff:max_pixel_cutoff);
%     end
%     end
    
    ScanProperties(j).Bscans11 = cell2mat(Bscans_RPE(j).y);
    ScanProperties(j).RPE = nan(length(ScanProperties(j).Bscans11),1);
    for ii = 1:length(ScanProperties(j).Bscans11)
        ScanProperties(j).RPE(ii) = str2double(ScanProperties(j).Bscans11(ii).Text); % convert RPE xml pixel data into struct
    if pixel_cutting
        ScanProperties(j).RPE_pixelcut = ScanProperties(j).RPE(min_pixel_cutoff:max_pixel_cutoff);
    end
    end
    
    ScanProperties(j).Bscans12 = cell2mat(Bscans_OB_RPE(j).y);
    ScanProperties(j).OB_RPE = nan(length(ScanProperties(j).Bscans12),1);
    for ii = 1:length(ScanProperties(j).Bscans12)
        ScanProperties(j).OB_RPE(ii) = str2double(ScanProperties(j).Bscans12(ii).Text); % convert OB_RPE xml pixel data into struct
    if pixel_cutting
        ScanProperties(j).OB_RPE_pixelcut = ScanProperties(j).OB_RPE(min_pixel_cutoff:max_pixel_cutoff);
    end
    end
end

if pixel_cutting
    x = (1:length(ScanProperties(j).ILM_pixelcut));
else
    x = (1:length(ScanProperties(j).ILM)); %#ok<UNRCH>
end
    % ILM rotation correction
    ILM_deg = 0; %**NEEDS MANUAL INPUT**% % This is the amount that tilted OCT scans need to be rotated by in order to achieve leveled scan. Usually will be between 145-155 deg for tilted scan. *TO DO:INPUT THIS AT THE BEGINNING OF SCRIPT*
    V_ILM = [x(:) ScanProperties(1).ILM_pixelcut zeros(size(ScanProperties(1).ILM_pixelcut(:)))];
    V_centre_ILM = mean(V_ILM,1); % centre of line
    Vc_ILM = V_ILM-ones(size(V_ILM,1),1)*V_centre_ILM; % centering coordinates
    
    a_rad_ILM = deg2rad(ILM_deg); % angle in radians
    E_ILM = [0  0 a_rad_ILM]; % euler angles for X,Y,Z-axis rotations
    Rz_ILM = [cos(E_ILM(3))  -sin(E_ILM(3)) 0;
        sin(E_ILM(3))  cos(E_ILM(3))  0;
        0        0        1]; % axis rotation
    Vrc_ILM = (Rz_ILM*Vc_ILM')'; % Centered coordinates for rotation
    Vr_ILM = Vrc_ILM + ones(size(V_ILM,1),1) * V_centre_ILM; % shifting back to original location
    figure
    mainplot = axes;
    plot(mainplot,Vr_ILM(:,1),Vr_ILM(:,2)); % check scan has been leveled in plot
    axis equal
    g = gcf;
    axesObjs_ILM = get(g,'Children'); % axis handles
    dataObjs_ILM = get(axesObjs_ILM,'Children'); % handles to low-level graphics objects in axes

    %BMEIS rotation correction
    BMEIS_deg = 0; %**NEEDS MANUAL INPUT**% % This is the amount that tilted OCT scans need to be rotated by in order to achieve leveled scan. Usually same as that for ILM. *TO DO:INPUT THIS AT THE BEGINNING OF SCRIPT*
    V_BMEIS = [x(:) ScanProperties(1).BMEIS_pixelcut zeros(size(ScanProperties(1).BMEIS_pixelcut(:)))];
    V_centre_BMEIS = mean(V_BMEIS,1); % centre of line
    Vc_BMEIS = V_BMEIS-ones(size(V_BMEIS,1),1)*V_centre_BMEIS; % centering coordinates
    a_rad_BMEIS = deg2rad(BMEIS_deg); % angle in radians
    E_BMEIS = [0  0 a_rad_BMEIS]; % euler angles for X,Y,Z-axis rotations
    Rz_BMEIS = [cos(E_BMEIS(3))  -sin(E_BMEIS(3)) 0;
        sin(E_BMEIS(3))  cos(E_BMEIS(3))  0;
        0        0        1]; % axis rotation
    Vrc_BMEIS = (Rz_BMEIS*Vc_BMEIS')'; % Centered coordinates for rotation
    Vr_BMEIS = Vrc_BMEIS + ones(size(V_BMEIS,1),1) * V_centre_BMEIS; % shifting back to original location
    hold on
    plot(mainplot,Vr_BMEIS(:,1),Vr_BMEIS(:,2),'k'); %#ok<NASGU> % check scan has been leveled in plot
    figure
    bm = plot(Vr_BMEIS(:,1),Vr_BMEIS(:,2),'k'); %#ok<NASGU> % check scan has been leveled in plot
    bm = gcf;
    axesObjs_BMEIS = get(bm,'Children'); % axis handles
    dataObjs_BMEIS = get(axesObjs_BMEIS,'Children'); % handles to low-level graphics objects in axes
    
    %RPE rotation correction
    RPE_deg = 0; %**NEEDS MANUAL INPUT**% % This is the amount that tilted OCT scans need to be rotated by in order to achieve leveled scan. Usually same as that for ILM. *TO DO:INPUT THIS AT THE BEGINNING OF SCRIPT*
    V_RPE = [x(:) ScanProperties(1).RPE_pixelcut zeros(size(ScanProperties(1).RPE_pixelcut(:)))];
    V_centre_RPE = mean(V_RPE,1); % centre of line
    Vc_RPE = V_RPE-ones(size(V_RPE,1),1)*V_centre_RPE; % centering coordinates
    a_rad_RPE = deg2rad(RPE_deg); % angle in radians
    E_RPE = [0  0 a_rad_RPE]; % euler angles for X,Y,Z-axis rotations
    Rz_RPE = [cos(E_RPE(3))  -sin(E_RPE(3)) 0;
        sin(E_RPE(3))  cos(E_RPE(3))  0;
        0        0        1]; % axis rotation
    Vrc_RPE = (Rz_RPE*Vc_RPE')'; % Centered coordinates for rotation
    Vr_RPE = Vrc_RPE + ones(size(V_RPE,1),1) * V_centre_RPE; % shifting back to original location
    hold on
    plot(mainplot, Vr_RPE(:,1),Vr_RPE(:,2)); %#ok<NASGU> % check scan has been leveled in plot
    figure
    rp = plot(Vr_RPE(:,1),Vr_RPE(:,2)); %#ok<NASGU> % check scan has been leveled in plot
    rp = gcf;
    axesObjs_RPE = get(rp,'Children'); % axis handles
    dataObjs_RPE = get(axesObjs_RPE,'Children'); % handles to low-level graphics objects in axes


%< GENERATE ANALYSIS>
[foveal_height_ILM, foveal_height_idx] = max(dataObjs_ILM.YData); % %**INPUT MANUALLY**%% calculate foveal_height and index within given region in the middle of the scan to avoid picking up artefacts
%foveal_height_idx = foveal_height_idx+300; %**MAKE SURE TO AMEND X INDEX FOR ABOVE SELECTION**%
%foveal_height_idx = 240;
%foveal_height_ILM = 130.8242;
foveal_height_x_index = dataObjs_ILM.XData(foveal_height_idx); % calculate x value at foveal height index for plotting
plot(mainplot,foveal_height_x_index,foveal_height_ILM,'r*')

corrected_XData_crest1 = dataObjs_ILM.XData(foveal_height_idx:(foveal_height_idx+200)); % go 100 pixels to left of the deepest foveal point for X (remember this is reversed in a rotated plot. The first crest point (1024 - 0) will be plotted to the left of the deepest foveal point)
corrected_YData_crest1 = dataObjs_ILM.YData(foveal_height_idx:(foveal_height_idx+200)); % go 100 pixels to left of the deepest foveal point for Y
[crest_height1, crest_height1_idx] = min(corrected_YData_crest1); % find the first max crest height
crest_height1_XvalILM = corrected_XData_crest1(crest_height1_idx); % find the x value at that corrected crest point for plotting and calculating foveal diameter.
hold on
plot(mainplot,crest_height1_XvalILM,crest_height1,'r*') % plot point of first max crest height for visualisation

corrected_XData_crest2 = dataObjs_ILM.XData((foveal_height_idx-200):foveal_height_idx); % go 100 pixels to right of the deepest foveal point for X (remember this is reversed in a rotated plot. The second crest point (0 - 1024) will be plotted to the right of the deepest foveal point)
corrected_YData_crest2 = dataObjs_ILM.YData((foveal_height_idx-200):foveal_height_idx); % go 100 pixels to right of the deepest foveal point for Y
[crest_height2, crest_height2_idx] = min(corrected_YData_crest2); % find the second max crest height
crest_height2_XvalILM = corrected_XData_crest2(crest_height2_idx); % find the x value at that crest point for plotting and calculating foveal diameter.
plot(mainplot,crest_height2_XvalILM,crest_height2,'r*') % plot point of second max crest height for visualisation

foveal_diameter = abs(diff(([crest_height1_XvalILM,crest_height2_XvalILM]))); % foveal diameter is difference between X crest values
plot(mainplot,[crest_height1_XvalILM,crest_height2_XvalILM],[crest_height1,crest_height2]); % plot foveal diameter distance for visualisation

figure
diam_line = plot([crest_height1_XvalILM,crest_height2_XvalILM],[crest_height1,crest_height2]);
diam_line = gcf;
axesObjs_diam_line = get(diam_line,'Children'); % axis handles
dataObjs_RPE_diam_line = get(axesObjs_diam_line,'Children'); % handles to low-level graphics objects in axes
interpolated_diam_line_Y = linspace(dataObjs_RPE_diam_line.YData(1), dataObjs_RPE_diam_line.YData(2), abs(crest_height1_XvalILM - crest_height2_XvalILM));
interpolated_diam_line_X = linspace(dataObjs_RPE_diam_line.XData(1), dataObjs_RPE_diam_line.XData(2), abs(crest_height1_XvalILM - crest_height2_XvalILM));

plot(mainplot,(crest_height1_XvalILM - crest_height1_idx),interpolated_diam_line_Y(crest_height1_idx),'r*')

foveal_depth = abs((interpolated_diam_line_Y(crest_height1_idx)) - foveal_height_ILM);
plot(mainplot, [(crest_height1_XvalILM - crest_height1_idx),foveal_height_idx], [interpolated_diam_line_Y(crest_height1_idx),foveal_height_ILM]);

foveal_height = abs(dataObjs_ILM.YData(foveal_height_idx) - dataObjs_BMEIS.YData(foveal_height_idx));
plot(mainplot, [foveal_height_idx,foveal_height_idx], [dataObjs_ILM.YData(foveal_height_idx),dataObjs_BMEIS.YData(foveal_height_idx)]);

% DELETED CODE (ORIGINAL)
% midpoint = [(crest_height1_XvalILM+crest_height2_XvalILM)/2,(crest_height1+crest_height2)/2]; % find midpoint of foveal diameter line
% plot(midpoint(1),midpoint(2),'r*') % plot midpoint of foveal diameter line for visualisation
% foveal_depth = abs(midpoint(2)-foveal_height_ILM); % calculate foveal depth
%foveal_height = abs(foveal_height_ILM - dataObjs_BMEIS(1).YData(foveal_height_idx));

if horizontal
    switch eye %#ok<UNRCH>
        case 'Left'
            actual_index_for_BMEIS = find(abs(dataObjs_ILM(1).YData - crest_height1) < 1e-100);
            parafoveal_retinalthickness_temporal = abs(crest_height1 - dataObjs_BMEIS(1).YData(actual_index_for_BMEIS(1))); % find temporal parafoveal retinal thickness for left eye (remember rotated image is inverted)
            plot(mainplot,[crest_height1_XvalILM,crest_height1_XvalILM],[crest_height1, dataObjs_BMEIS(1).YData(round(median(actual_index_for_BMEIS)))],'y');
        case 'Right'
            actual_index_for_BMEIS = find(abs(dataObjs_ILM(1).YData - crest_height1) < 1e-100);
            parafoveal_retinalthickness_nasal = abs(crest_height1 - dataObjs_BMEIS(1).YData(actual_index_for_BMEIS(1))); % find nasal parafoveal retinal thickness for right eye (remember rotated image is inverted)
            plot(mainplot,[crest_height1_XvalILM,crest_height1_XvalILM],[crest_height1, dataObjs_BMEIS(1).YData(round(median(actual_index_for_BMEIS)))],'y');
    end
else
    switch eye
        case 'Left'
            actual_index_for_BMEIS = find(abs(dataObjs_ILM(1).YData - crest_height1) < 1e-100);
            parafoveal_retinalthickness_inferior = abs(crest_height1 - dataObjs_BMEIS(1).YData(actual_index_for_BMEIS(1))); % find inferior parafoveal retinal thickness for left eye (remember rotated image is inverted)
            plot(mainplot,[crest_height1_XvalILM,crest_height1_XvalILM],[crest_height1, dataObjs_BMEIS(1).YData(round(median(actual_index_for_BMEIS)))],'y');
        case 'Right'
            actual_index_for_BMEIS = find(abs(dataObjs_ILM(1).YData - crest_height1) < 1e-100);
            parafoveal_retinalthickness_inferior = abs(crest_height1 - dataObjs_BMEIS(1).YData(actual_index_for_BMEIS(1))); % find inferior parafoveal retinal thickness for right eye (remember rotated image is inverted)
            plot(mainplot,[crest_height1_XvalILM,crest_height1_XvalILM],[crest_height1, dataObjs_BMEIS(1).YData(round(median(actual_index_for_BMEIS)))],'y');
    end
end

if horizontal
    switch eye %#ok<UNRCH>
        case 'Left'
            actual_index_for_BMEIS = find(abs(dataObjs_ILM(1).YData - crest_height2) < 1e-100);
            parafoveal_retinalthickness_nasal = abs(crest_height2 - dataObjs_RPE(1).YData(actual_index_for_BMEIS(1))); % find nasal parafoveal retinal thickness for left eye (remember rotated image is inverted)
            plot(mainplot,[crest_height2_XvalILM,crest_height2_XvalILM],[crest_height2, dataObjs_BMEIS(1).YData(round(median(actual_index_for_BMEIS)))],'y');
        case 'Right'
            actual_index_for_BMEIS = find(abs(dataObjs_ILM(1).YData - crest_height2) < 1e-100);
            parafoveal_retinalthickness_temporal = abs(crest_height2 - dataObjs_RPE(1).YData(actual_index_for_BMEIS(1))); % find temporal parafoveal retinal thickness for right eye (remember rotated image is inverted)
            plot(mainplot,[crest_height2_XvalILM,crest_height2_XvalILM],[crest_height2, dataObjs_BMEIS(1).YData(round(median(actual_index_for_BMEIS)))],'y');
    end
else
    switch eye
        case 'Left'
            actual_index_for_BMEIS = find(abs(dataObjs_ILM(1).YData - crest_height2) < 1e-100);
            parafoveal_retinalthickness_superior = abs(crest_height2 - dataObjs_RPE(1).YData(actual_index_for_BMEIS(1))); % find superior parafoveal retinal thickness for left eye (remember rotated image is inverted)
            plot(mainplot,[crest_height2_XvalILM,crest_height2_XvalILM],[crest_height2, dataObjs_BMEIS(1).YData(round(median(actual_index_for_BMEIS)))],'y');
        case 'Right'
            actual_index_for_BMEIS = find(abs(dataObjs_ILM(1).YData - crest_height2) < 1e-100);
            parafoveal_retinalthickness_superior = abs(crest_height2 - dataObjs_RPE(1).YData(actual_index_for_BMEIS(1))); % find superior parafoveal retinal thickness for left eye (remember rotated image is inverted)
            plot(mainplot,[crest_height2_XvalILM,crest_height2_XvalILM],[crest_height2, dataObjs_BMEIS(1).YData(round(median(actual_index_for_BMEIS)))],'y');
    end
end

% LayerThickness = struct('difference_RNFL',[],'difference_GCL',[],'difference_IPL',[],'difference_INL',[],'difference_OPL',[],'difference_IS',[],'difference_ISOSJ',[],'difference_OS',[],'difference_RPE',[]);
% 
% for i = 1:length(ScanProperties(1).RNFL)
%     LayerThickness.difference_RNFL(i) = ScanProperties(1).RNFL(i) - ScanProperties(1).ILM(i);
%     RNFL_thickness = mean([LayerThickness.difference_RNFL]);
%     
%     LayerThickness.difference_GCL(i) = ScanProperties(1).GCL(i) - ScanProperties(1).RNFL(i);
%     GCL_thickness = mean([LayerThickness.difference_GCL]);
%     
%     LayerThickness.difference_IPL(i) = ScanProperties(1).IPL(i) - ScanProperties(1).GCL(i);
%     IPL_thickness = mean([LayerThickness.difference_IPL]);
%     
%     LayerThickness.difference_INL(i) = ScanProperties(1).INL(i) - ScanProperties(1).IPL(i);
%     INL_thickness = mean([LayerThickness.difference_INL]);
%     
%     LayerThickness.difference_OPL(i) = ScanProperties(1).OPL(i) - ScanProperties(1).INL(i);
%     OPL_thickness = mean([LayerThickness.difference_OPL]);
%     
%     LayerThickness.difference_ONL(i) = ScanProperties(1).BMEIS(i) - ScanProperties(1).OPL(i);
%     ONL_thickness = mean([LayerThickness.difference_ONL]);
%     
%     LayerThickness.difference_IS_OS(i) = ScanProperties(1).ISOSJ(i) - ScanProperties(1).BMEIS(i);
%     IS_OSthickness = mean([LayerThickness.difference_ISOSJ]);
%     
%     LayerThickness.difference_ISOSJ(i) = ScanProperties(1).IB_OPR(i) - ScanProperties(1).ISOSJ(i);
%     InnerSegmentOuterSegmentJunction_thickness = mean([LayerThickness.difference_ISOSJ]);
%     
%     LayerThickness.difference_RPE(i) = ScanProperties(1).OB_RPE(i) - ScanProperties(1).RPE(i);
%     RPE_thickness = mean([LayerThickness.difference_RPE]);
% end

if horizontal
    report(:,1) = foveal_height; %#ok<UNRCH>
    report(:,2) = parafoveal_retinalthickness_temporal;
    report(:,3) = parafoveal_retinalthickness_nasal;
    report(:,4) = foveal_depth;
    report(:,5) = foveal_diameter;
%     report(:,6) = RNFL_thickness;
%     report(:,7) = GCL_thickness;
%     report(:,8) = IPL_thickness;
%     report(:,9) = INL_thickness;
%     report(:,10) = OPL_thickness;
%     report(:,11) = ONL_thickness;
%     report(:,12) = IS_OSthickness;
%     report(:,13) = InnerSegmentOuterSegmentJunction_thickness;
%     report(:,14) = RPE_thickness;
    
    fileName = [participantId,'Foveal_parameters.csv'];
    header = 'foveal_height,parafoveal_retinalthickness_temporal,parafoveal_retinalthickness_nasal,foveal_depth,RNFL_thickness,GCL_thickness,IPL_thickness,INL_thickness,OPL_thickness,IS_thickness,InnerSegmentOuterSegmentJunction_thickness,OS_thickness,RPE_thickness,foveal_diameter';
    generatecsv(report,fileName,header); % export report
    
else
    report(:,1) = foveal_height; %#ok<UNRCH>
    report(:,2) = parafoveal_retinalthickness_inferior;
    report(:,3) = parafoveal_retinalthickness_superior;
    report(:,4) = foveal_depth;
    report(:,5) = foveal_diameter;
%     report(:,6) = RNFL_thickness;
%     report(:,7) = GCL_thickness;
%     report(:,8) = IPL_thickness;
%     report(:,9) = INL_thickness;
%     report(:,10) = OPL_thickness;
%     report(:,11) = ONL_thickness;
%     report(:,12) = IS_OSthickness;
%     report(:,13) = InnerSegmentOuterSegmentJunction_thickness;
%     report(:,14) = RPE_thickness;
    
    fileName = [participantId,'Foveal_parameters.csv'];
    header = 'foveal_height,parafoveal_retinalthickness_inferior,parafoveal_retinalthickness_superior,foveal_depth,RNFL_thickness,GCL_thickness,IPL_thickness,INL_thickness,OPL_thickness,ONL_thickness,IS/OS_thickness,InnerSegmentOuterSegmentJunction_thickness,RPE_thickness,foveal_diameter';
    generatecsv(report,fileName,header); % export report
end
