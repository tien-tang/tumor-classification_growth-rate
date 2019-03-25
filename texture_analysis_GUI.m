%%%% All this at the beginning is just pre-made code that comes with
%%%% setting up the GUI:

function varargout = texture_analysis_GUI(varargin)
% TEXTURE_ANALYSIS_GUI MATLAB code for texture_analysis_GUI.fig
%      TEXTURE_ANALYSIS_GUI, by itself, creates a new TEXTURE_ANALYSIS_GUI or raises the existing
%      singleton*.
%
%      H = TEXTURE_ANALYSIS_GUI returns the handle to a new TEXTURE_ANALYSIS_GUI or the handle to
%      the existing singleton*.
%
%      TEXTURE_ANALYSIS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEXTURE_ANALYSIS_GUI.M with the given input arguments.
%
%      TEXTURE_ANALYSIS_GUI('Property','Value',...) creates a new TEXTURE_ANALYSIS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before texture_analysis_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to texture_analysis_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help texture_analysis_GUI

% Last Modified by GUIDE v2.5 22-Sep-2017 13:41:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @texture_analysis_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @texture_analysis_GUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before texture_analysis_GUI is made visible.
function texture_analysis_GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to texture_analysis_GUI (see VARARGIN)

% Choose default command line output for texture_analysis_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles.single = struct('imagestore',0,'mask',0,'on_off',0);
guidata(hObject,handles)

% UIWAIT makes texture_analysis_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = texture_analysis_GUI_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%%% Start of code:


%%%%%   SINGLE IMAGE   %%%%%

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, ~, handles)
%CHOOSE SINGLE IMAGE FILE
% Press this button to select single image file and then make an roi on it

% Clear out any data currently saved in handles.single
handles.single = []; 

% Set up data structure
handles.single.imagestore = struct('filename',[],'image',[], ...
    'mask',[],'roi0',[]);

% User interface to choose file (limited to .jpg and .dcm format)
% Modify this section to look for other image file type
[handles.single.imagestore.filename,path] = uigetfile('*.jpg;*.dcm',...
    'Please select an image file. Must be JPEG or DICOM file format.');

% Use filename from uigetfile; read file
name = handles.single.imagestore.filename;
if min(name(length(name) - 3 : length(name)) == '.jpg')
    I = imread(strcat(path,handles.single.imagestore.filename)); % Read file
elseif min(name(length(name) - 3 : length(name)) == '.dcm')
    I = dicomread(strcat(path,handles.single.imagestore.filename)); % Read file 
end
handles.single.imagestore.image = I;

% User will draw roi on image in GUI; mask variable will store roi info

%%%% The following code will use the default value of "imadjust" to enhance
%%%% image contrast by mapping the intensity value to new value so that 1%
%%%% of the data is saturated at low and high intensities. The WILL CHANGE
%%%% the pixel value.

% I = imadjust(I);
% mask = roipoly(I);

%%%% The following code will display the image with the default range, use
%%%% this with the MRI

imshow(I,[])
mask = roipoly;

%%%% Code should only run with one or the other above lines of code commented out
I2 = I;
I2(~mask) = 0; % Make new image with zeros outside ROI

% The following series of for loops go through the rows and columns of the
% roi image and finds the first and last row and column that doesn't
% contain only zeros. This is used to "cut off" the extra space around the
% roi, allowing the image to zoom in and extraneous info to be excluded
% from calculations.

for k = 1:size(I2,1)
    if max(I2(k,:)) ~= 0
        r1 = k; % Find first row that isn't all zeros
        break
    end
end
for j = 1:size(I2,2)
    if max(I2(:,j)) ~= 0
        c1 = j; % Find first column that isn't all zeros
        break
    end
end
for k = 1:size(I2,1)
    kk = size(I2,1)+1-k;
    if max(I2(kk,:)) ~= 0
        r2 = kk; % Find last row that isn't all zeros
        break
    end
end
for j = 1:size(I2,2)
    jj = size(I2,2)+1-j;
    if max(I2(:,jj)) ~= 0
        c2 = jj; % Find last column that isn't all zeros
        break
    end
end

I2 = I2(r1:r2,c1:c2); % Cut off all the rows/columns that are all 0
handles.single.imagestore.mask = mask(r1:r2,c1:c2); % Store trimmed mask
handles.single.imagestore.roi0 = I2; % Store trimmed roi image
imshow(I2,[]) % Display trimmed roi image in GUI

handles.single.on_off = 1; % Turn "single switch" on (basically just keeps
% track of the fact that the user is doing a
% single image)
guidata(hObject,handles)

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, ~, handles)
%ANALYZE SINGLE ROI
% Call the texture_analysis function after selecting a single image; return
% results structure

if handles.single.on_off ~= 1 %check if "single switch" is on
    errordlg('Select single image to analyze.')
else
    if isempty(handles.single.imagestore.roi0) % Check that user has made ROI
        errordlg('Please select ROI.')
    else
        % This runs the texture_analysis function on the single image, and
        % saves the results in the handles structure:
        [handles.single.results64,handles.single.results128,handles.single.results256,handles.single.results512]...
            =texture_analysis(handles.single.imagestore.filename,...
            handles.single.imagestore.roi0,handles.single.imagestore.mask);
        fprintf('\nProcessing...\n')
    end
end
msgbox('Processing Completed');
guidata(hObject,handles)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(~, ~, handles)
%EXPORT DATA (SINGLE)
exportdata(handles.single.results64,handles.single.results128,handles.single.results256,handles.single.results512,handles)


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(~, ~, handles)
% VIEW RESULTS (SINGLE FILE)
% This allows the user to view results within Matlab using the GUI. When
% pressed, two figures pop up -- one with a table of values, and the other
% displaying either the histogram or the GLCM plot, depending on what is
% selected in the drop down box.

combined = [handles.single.results64,handles.single.results128,handles.single.results256,handles.single.results512];

featurenames = fieldnames(handles.single.results64); % Name of file

% Creat the data structure that sperates the 1st and 2nd order image
% features
data = struct2cell(combined);
data(1,:)=data{1,:};
data1 = data(1:12,:);
data2 = cell2mat(data(18:36,:));

offset = {' 0 deg - 256',' 45 deg',' 90 deg',' 135 deg', ' 0 deg - 512',...
    ' 45 deg',' 90 deg',' 135 deg',' 0 deg - 1024',' 45 deg',' 90 deg',...
    ' 135 deg',' 0 deg - 2048',' 45 deg',' 90 deg',' 135 deg'};

figure % Create figure with table of imaging feature values:
t1 = uitable('Data', data1, 'RowName', featurenames(1:12,:), 'ColumnName',...
    {'256','512','1024','2048'},'unit','normalized','Position', [0 0 1 1],'ColumnWidth',{150});

figure
t2 = uitable('Data', data2, 'RowName', featurenames(18:36,:),'ColumnName',offset,...
    'unit','normalized','Position', [0 0 1 1]);

figure % Create another figure that lets you select a graph to show:
results = handles.single.results64;
name = results.name{1};
name = name(1:(length(name)-4));

% Create histogram (initially, when it pops up, histogram is shown)
h = axes;
hist(results.roivec,32)
title(strcat(name,' Results: Histogram'),'Interpreter','none')
set(h,'Position',[0.13 0.11 0.775 0.7])

% This sets up the drop down menu:
uicontrol('Style', 'popup',...
    'String', 'Histogram|GLCM 16|GLCM 32|GLCM 64|GLCM 128',...
    'Position', [20 340 100 50],...
    'Callback', {@switchdisplay,handles,h});

function switchdisplay(hObj,~,handles,h)
% Called when user activates popup menu (allows for switching between
% histogram and GLCM plots)

val = get(hObj,'Value'); % Determine which option is selected
results = handles.single.results64;
name = results.name{1};
name = name(1:(length(name)-4));
glcmmean64 = mean(results.glcm64,3);
glcmmean128 = mean(results.glcm128,3);
glcmmean256 = mean(results.glcm256,3);
glcmmean512 = mean(results.glcm512,3);

if val == 1
    hist(results.roivec)
    title(strcat(name,' Results: Histogram'),'Interpreter','none')
    set(h,'Position',[0.13 0.11 0.775 0.7])
elseif val == 2
    surf(glcmmean64,'FaceAlpha',0.5)
    colorbar
    title(strcat(name,' Results: Mean GLCM - 256 numlevel'),'Interpreter','none')
elseif val == 3
    surf(glcmmean128,'FaceAlpha',0.5)
    colorbar
    title(strcat(name,' Results: Mean GLCM - 512 numlevel'),'Interpreter','none')
elseif val == 4
    surf(glcmmean256,'FaceAlpha',0.5)
    colorbar
    title(strcat(name,' Results: Mean GLCM - 1024 numlevel'),'Interpreter','none')
else
    surf(glcmmean512,'FaceAlpha',0.5)
    colorbar
    title(strcat(name,' Results: Mean GLCM - 2048 numlevel'),'Interpreter','none')
end


%%%%%   FOLDER OF IMAGES 

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, ~, handles)
%CHOOSE FOLDER OF IMAGE FILES
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.n = [];
% handles.imagestore = [];
% handles.cancel = 0;
% handles.i = [];
% handles.results = [];

handles.single.on_off = 0; % Turn off single switch

% Select folder containing image files to analyze
folder = uigetdir('','Select Folder Containing Image Files. Image files must be in JPEG or DICOM file format.');
listing = dir(folder); % Get list of files in folder
filelist = {}; % Initialize array

% Go through the files in selected folder and make list of JPEG files
for i = 3:length(listing) %(the first two lines of listing are always "." and ".."; not sure why)
    name = listing(i).name;
    if min(name(length(name) - 3 : length(name)) == '.jpg') || ...
            min(name(length(name) - 3 : length(name)) == '.dcm') % Check if this file is JPEG or DICOM
        filelist = vertcat(filelist,name); % If so, add to filelist
    end
end

n = length(filelist); % Number of files to be analyzed
handles.n = n;

% Initialize data structure; create fields and input list of file names
handles.imagestore = struct('filename',filelist,'image',[], ...
    'mask',[],'roi0',[]);

handles.cancel = 0; % Boolean; i.e. 0 means don't cancel

% handles.i will be used in various functions to iterate through the list
% of image files.
handles.i = 1; % Initialize iteration
guidata(hObject,handles)

% This while loop checks if the user has "cancelled" (i.e. wants to stop
% selecting ROIs/analyzing that set of images). It also checks the
% iteration (handles.i), and stops when it has iterated over all of the
% images.
while handles.cancel == 0 && handles.i <= handles.n
    i = handles.i;
    for j = 1:n
        if handles.cancel == 1
            break
        end
    end
    
    %%% Without some form of "cancel", if you exit in the middle of this,
    %%% it will keep popping up with the images in new figure windows until
    %%% it has gone through all the images.
    %%% In the future, maybe look up a way to make it cancel when you close
    %%% the GUI? Maybe use waitfor function?
    
    imagename = filelist{i}; % Take image file name
    handles.imagestore(i).name = imagename; % Save image file name in structure
    if min(imagename(length(imagename) - 3 : length(imagename)) == '.jpg')
        I = imread(strcat(folder,'\',imagename)); % Read file
    elseif min(imagename(length(imagename) - 3 : length(imagename)) == '.dcm')
        I = dicomread(strcat(folder,'\',imagename)); % Read file
    end
    
    handles.imagestore(i).image = I; % Save image to structure
    % User will draw roi on image in GUI; mask variable will store roi info
    
    %%%% The following code will use the default value of "imadjust" to enhance
    %%%% image contrast by mapping the intensity value to new value so that 1%
    %%%% of the data is saturated at low and high intensities. That WILL CHANGE
    %%%% the pixel value.
    
    % I = imadjust(I);
    % mask = roupoly(I);
    
    %%%% The following code will display the image with the default range, use
    %%%% this with the MRI scans
    
    imshow(I,[])
    mask = roipoly;
    
    % % %Use one fo the above display code for ROI extraction 
    
    I2 = I;
    I2(~mask) = 0; % Make part of image not in ROI = 0
    
    % The following series of for loops go through the rows and columns of the
    % roi image and finds the first and last row and column that doesn't
    % contain only zeros. This is used to "cut off" the extra space around the
    % roi, allowing the image to zoom in and extraneous info to be excluded
    % from calculations. (Same as for single image, but inside while loop)
    
    for k = 1:size(I2,1)
        if max(I2(k,:)) ~= 0
            r1 = k; % Find first row that isn't all zeros
            break
        end
    end
    for j = 1:size(I2,2)
        if max(I2(:,j)) ~= 0
            c1 = j; % Find first column that isn't all zeros
            break
        end
    end
    for k = 1:size(I2,1)
        kk = size(I2,1)+1-k;
        if max(I2(kk,:)) ~= 0
            r2 = kk; % Find last row that isn't all zeros
            break
        end
    end
    for j = 1:size(I2,2)
        jj = size(I2,2)+1-j;
        if max(I2(:,jj)) ~= 0
            c2 = jj; % Find last column that isn't all zeros
            break
        end
    end
    
    I2 = I2(r1:r2,c1:c2); % Cut off all the rows/columns that are all 0
    handles.imagestore(i).mask = mask(r1:r2,c1:c2);
    handles.imagestore(i).roi0 = I2;
    imshow(I2,[]) % Show new, cropped roi image
    
    guidata(hObject,handles)
    
    uiwait % Wait to go to next image until after "next" button is pressed
    handles.i = i + 1; % Increment for next iteration
    guidata(hObject,handles)
end
guidata(hObject,handles)

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(~, ~, handles)
%NEXT FILE
% This just makes it go to the next file for choosing ROIs

i = handles.i; % Get current iteration number
n = handles.n; % Get total number of images

% These if statementS ensure that an roi was selected before moving on to
% the next image, and informs the user if they've reached the last image.
if ~isempty(handles.imagestore(i).roi0) %check if roi was selected
    if i < n % Check if at last image
        uiresume
    else
        errordlg('No more image files in list.')
    end
else
    errordlg('Please select ROI')
    uiwait
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, ~, handles)
%CANCEL
% This is supposed to prevent Matlab from popping up with all the other
% images if you don't want to continue selecting ROIs. Might want to look
% for better way to do this eventually.
handles.cancel = 1;
guidata(hObject,handles)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, ~, handles)
%ANALYZE ALL ROIS
% Call the texture_analysis function after choosing ROIs for each of list
% of images; return results structur
if handles.single.on_off == 1 % Check if "single switch" is on
    errordlg('Select folder of image files to analyze.')
else
    if isempty(handles.imagestore(1).roi0) % Check if ROIs have been selected
        errordlg('Please select ROIs.')
    else
        fprintf('\nProcessing...\n')
        for i = 1:handles.i
            % This runs the texture_analysis function on each image, and
            % saves the results in the handles structure:
            [handles.results64(i),handles.results128(i),handles.results256(i),handles.results512(i)]=...
              texture_analysis(handles.imagestore(i).filename,...
                handles.imagestore(i).roi0,handles.imagestore(i).mask);
        end
    end
end
msgbox('Processing Completed');
guidata(hObject,handles)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(~, ~, handles)
%EXPORT DATA (LIST)
exportdata(handles.results64,handles.results128,handles.results256,handles.results512,handles)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(~, ~, handles)
% VIEW RESULTS (MULTIPLE FILES)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This creates a new interactive figure that displays the results of the
% analysis, allowing the user to select which data to view.

results = struct2cell(handles.results64);
firstorder = results(1:12,:,:);
a = permute(firstorder,[1 3 2]);
a(1,:)=[a{1,:,:}];
 
featurenames = fieldnames(handles.results64); 

figure 
t1 = uitable('Data', a, 'RowName', featurenames(1:12,:), 'ColumnName',...
    {'32','256'},'unit','normalized','Position', [0 0 1 1],'ColumnWidth',{150});


%%%%%   TEXTURE FEATURE EXTRACTION 


function [results64,results128,results256,results512] = texture_analysis(imagename, roi0, mask)
% Outputs a structure with all the texture analysis info for inputted image
%input:
%   imagename: whatever you want to call the image
%   roi0: image with 0 outside of roi, since its not possible to get NaN
%   with uint16 or other format to have NaN
%   mask: image with 0 outside of roi
%output:
%   results: structure containing 1st and 2nd order image features

% Create structure to contain results of analysis of 1st order features
results = struct('name',0,'numpixels',0,'mean',0,'std',0,'median',0,'max',0,'min',0,...
    'range',0,'variance',0,'skew',0,'kurtosis',0,'zernike',0);
results.name = {imagename}; %record name
imagedoubleNaN = im2double(roi0); %convert images to doubles for analysis
imagedoubleNaN(~mask) = NaN;

%%%% 1ST ORDER TEXTURE ANALYSIS %%%%

% % To accurately calculate the mean, all elements are put into a single
% % vector, and then the mean (ignoring NaNs) of that vector is calculated
ntot = size(imagedoubleNaN,1)*size(imagedoubleNaN,2);
roiNaNvec = reshape(imagedoubleNaN,ntot,1);
results.numpixels = sum(~isnan(roiNaNvec));

%%% 1st ORDER FEATURES WITH IMAGE REPRESENTED AS DOUBLE PRECISION NUMBERS
%%% WHERE THE IMAGE IS SCALED FROM 0-1.

% Mean pixel intensity converted by the double precision numbers
results.mean = nanmean(roiNaNvec);

% Standard Deviation
results.std = nanstd(roiNaNvec);

% Median
results.median = nanmedian(roiNaNvec);

% Range, need to calculate the min and max
results.max = nanmax(roiNaNvec);
results.min = nanmin(roiNaNvec);
results.range = range(roiNaNvec);

% Variance
results.variance = nanvar(roiNaNvec);

% Histogram Data
% In order for a histogram to be constructed with the Export Results button
% this function must create a vector without any NaNs. Here, it takes out
% all the NaNs from the vector and saves it in the results structure:
j = 1; %initiate counter
for i = 1:ntot
    if ~isnan(roiNaNvec(i))
        roivec(j) = roiNaNvec(i);
        j = j+1;
    end
end
results.roivec = roivec;

% Skewness
% From Matlab Help document:
%   "S = skewness(X) returns the sample skewness of the values in X.  For a
%   vector input, S is the third central moment of X, divided by the cube
%   of its standard deviation."
results.skew = skewness(roivec);

% Kurtosis
% From Matlab Help document:
%    "K = kurtosis(X) returns the sample kurtosis of the values in X.  For a
%    vector input, K is the fourth central moment of X, divided by fourth
%    power of its standard deviation."
results.kurtosis = kurtosis(roivec);


%%%% 2ND ORDER TEXTURE ANALYSIS %%%%

% %Turn off the following warning: Warning: GLCM does not count
% %pixel pairs if either of their values is NaN.
id = 'images:graycomatrix:scaledImageContainsNan';
warning('off', id)


% Gray Level Co-occurrence Matrix (GLCM)
% Construct/save glcm; call function to extract out all the features
% *Note: default for graycomatrix() is to use 8 gray levels, but can set
% 'NumLevels' parameter to different integer for more/fewer levels.
% glcm = graycomatrix(imagedoubleNaN,'NumLevels',256,'Offset',...
%     [0 1; -1 1; -1 0; -1 -1],'Symmetric',true);
glcm64 = graycomatrix(imagedoubleNaN,'NumLevels',8,'Offset',...
    [0 1; -1 1; -1 0; -1 -1],'Symmetric',true);
results.glcm64 = glcm64;
stats64 = GLCM_Features2(glcm64,0); % Call function to extract texture features from glcm

glcm128 = graycomatrix(imagedoubleNaN,'NumLevels',16,'Offset',...
    [0 1; -1 1; -1 0; -1 -1],'Symmetric',true);
results.glcm128=glcm128;
stats128 = GLCM_Features2(glcm128,0);

glcm256 = graycomatrix(imagedoubleNaN,'NumLevels',32,'Offset',...
    [0 1; -1 1; -1 0; -1 -1],'Symmetric',true);
results.glcm256=glcm256;
stats256 = GLCM_Features2(glcm256,0);

glcm512 = graycomatrix(imagedoubleNaN,'NumLevels',512,'Offset',...
    [0 1; -1 1; -1 0; -1 -1],'Symmetric',true);
results.glcm512=glcm512;
stats512 = GLCM_Features2(glcm512,0);



% % MORPHOLOGICAL ANALYSIS

% Calculate Zernike moment (4th order, repetition of 2):
n = 4; m = 2;
[~, results.zernike, ~] = Zernikmoment(mask,n,m);


% % OUTPUTS

% Add structure containing glcm stats to structure containing other results
fields = [fieldnames(results); fieldnames(stats64)];

results64 = cell2struct([struct2cell(results); struct2cell(stats64)], fields, 1);

results128 = cell2struct([struct2cell(results); struct2cell(stats128)], fields, 1);

results256 = cell2struct([struct2cell(results); struct2cell(stats256)], fields, 1);

results512 = cell2struct([struct2cell(results); struct2cell(stats512)], fields, 1);



function [out] = GLCM_Features2(glcmin,pairs)
%
% GLCM_Features2 helps to calculate the features from the different GLCMs
% that are input to the function. The GLCMs are stored in a i x j x n
% matrix, where n is the number of GLCMs calculated usually due to the
% different orientation and displacements used in the algorithm. Usually
% the values i and j are equal to 'NumLevels' parameter of the GLCM
% computing function graycomatrix(). 
%
% Features computed
% Autocorrelation: [2]                      (out.autocorrelation)
% Contrast: matlab/[1,2]                    (out.contrast)
% Correlation: matlab                       (out.correlationmatlab)
% Correlation: [1,2]                        (out.corr)
% Cluster Prominence: [2]                   (out.cprom)
% Cluster Shade: [2]                        (out.cshad)
% Dissimilarity: [2]                        (out.dissi)
% Energy: matlab / [1,2]                    (out.energy)
% Entropy: [2]                              (out.entropy)
% Homogeneity: matlab                       (out.hommatlab)
% Homogeneity: [2]                          (out.homop)
% Maximum probability: [2]                  (out.maxpr)
% Sum of squares: Variance [1]              (out.sosvh)
% Sum average [1]                           (out.savgh)
% Sum variance [1]                          (out.svarh)
% Sum entropy [1]                           (out.senth)
% Difference variance [1]                   (out.dvarh)
% Difference entropy [1]                    (out.denth)
% Information measure of correlation1 [1]   (out.inf1h)
% Informaiton measure of correlation2 [1]   (out.inf2h)
% Inverse difference (INV) is homom [3]     (out.hommatlab)
% Inverse difference normalized (INN) [3]   (out.indnc)
% Inverse difference moment normalized [3]  (out.idmnc)
%
% Formulae from MATLAB site (some look different from
% the paper by Haralick but are equivalent and give same results)
% Example formulae:
% Contrast = sum_i(sum_j(  (i-j)^2 * p(i,j) ) ) (same in matlab/paper)
% Correlation = sum_i( sum_j( (i - u_i)(j - u_j)p(i,j)/(s_i.s_j) ) ) (m)
% Correlation = sum_i( sum_j( ((ij)p(i,j) - u_x.u_y) / (s_x.s_y) ) ) (p[2])
% Energy = sum_i( sum_j( p(i,j)^2 ) )           (same in matlab/paper)
% Homogeneity = sum_i( sum_j( p(i,j) / (1 + |i-j|) ) ) (as in matlab)
% Homogeneity = sum_i( sum_j( p(i,j) / (1 + (i-j)^2) ) ) (as in paper)
%
% Where:
% u_i = u_x = sum_i( sum_j( i.p(i,j) ) ) (in paper [2])
% u_j = u_y = sum_i( sum_j( j.p(i,j) ) ) (in paper [2])
% s_i = s_x = sum_i( sum_j( (i - u_x)^2.p(i,j) ) ) (in paper [2])
% s_j = s_y = sum_i( sum_j( (j - u_y)^2.p(i,j) ) ) (in paper [2])
%
%
% Normalize the glcm:
% Compute the sum of all the values in each glcm in the array and divide
% each element by its sum
%
% Haralick uses 'Symmetric' = true in computing the glcm
% no Symmetric flag in the Matlab version 
% add the diagonally opposite pairs to obtain the Haralick glcm
% Here it is assumed that the diagonally opposite orientations are paired
% one after the other in the matrix
% If the above assumption is true with respect to the input glcm then
% setting the flag 'pairs' to 1 will compute the final glcms that would result
% by setting 'Symmetric' to true. If your glcm is computed using the
% Matlab version with 'Symmetric' flag you can set the flag 'pairs' to 0
%
% References:
% 1. R. M. Haralick, K. Shanmugam, and I. Dinstein, Textural Features of
% Image Classification, IEEE Transactions on Systems, Man and Cybernetics,
% vol. SMC-3, no. 6, Nov. 1973
% 2. L. Soh and C. Tsatsoulis, Texture Analysis of SAR Sea Ice Imagery
% Using Gray Level Co-Occurrence Matrices, IEEE Transactions on Geoscience
% and Remote Sensing, vol. 37, no. 2, March 1999.
% 3. D A. Clausi, An analysis of co-occurrence texture statistics as a
% function of grey level quantization, Can. J. Remote Sensing, vol. 28, no.
% 1, pp. 45-62, 2002
% 4. http://murphylab.web.cmu.edu/publications/boland/boland_node26.html

% Following code has been pieced together by stuff found online and through
% matlab

% If 'pairs' not entered: set pairs to 0
if ((nargin > 2) || (nargin == 0))
    error('Too many or too few input arguments. Enter GLCM and pairs.');
elseif ( (nargin == 2) )
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
        error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
        error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end
elseif (nargin == 1) % only GLCM is entered
    pairs = 0; % default is numbers and input 1 for percentage
    if ((size(glcmin,1) <= 1) || (size(glcmin,2) <= 1))
        error('The GLCM should be a 2-D or 3-D matrix.');
    elseif ( size(glcmin,1) ~= size(glcmin,2) )
        error('Each GLCM should be square with NumLevels rows and NumLevels cols');
    end
end


format long e
%%% I'm unclear on what this is supposed to do, but it doesn't seem to work
if (pairs == 1)
    newn = 1;
    for nglcm = 1:2:size(glcmin,3)
        glcm(:,:,newn)  = glcmin(:,:,nglcm) + glcmin(:,:,nglcm+1);
        newn = newn + 1;
    end
elseif (pairs == 0)
    glcm = glcmin;
end

size_glcm_1 = size(glcm,1);
size_glcm_2 = size(glcm,2);
size_glcm_3 = size(glcm,3);

% checked
%%% This part is simply setting up the structure fields with a zero vector
%%% with as many elements as there are GLCM's (i.e. different offsets)
out.autocorrelation = zeros(1,size_glcm_3); % Autocorrelation: [2]
out.contrast = zeros(1,size_glcm_3); % Contrast: matlab/[1,2]
out.correlationmatlab = zeros(1,size_glcm_3); % Correlation: matlab
out.corr = zeros(1,size_glcm_3); % Correlation: [1,2]
out.cprom = zeros(1,size_glcm_3); % Cluster Prominence: [2]
out.cshad = zeros(1,size_glcm_3); % Cluster Shade: [2]
out.dissi = zeros(1,size_glcm_3); % Dissimilarity: [2]
out.energy = zeros(1,size_glcm_3); % Energy: matlab / [1,2]
out.entropy = zeros(1,size_glcm_3); % Entropy: [2]
out.hommatlab = zeros(1,size_glcm_3); % Homogeneity: matlab
out.homop = zeros(1,size_glcm_3); % Homogeneity: [2]
out.maxpr = zeros(1,size_glcm_3); % Maximum probability: [2]
out.sosvh = zeros(1,size_glcm_3); % Sum of sqaures: Variance [1]
out.savgh = zeros(1,size_glcm_3); % Sum average [1]
out.svarh = zeros(1,size_glcm_3); % Sum variance [1]
out.senth = zeros(1,size_glcm_3); % Sum entropy [1]
out.dvarh = zeros(1,size_glcm_3); % Difference variance [4]
out.denth = zeros(1,size_glcm_3); % Difference entropy [1]
out.inf1h = zeros(1,size_glcm_3); % Information measure of correlation1 [1]
out.inf2h = zeros(1,size_glcm_3); % Informaiton measure of correlation2 [1]
out.indnc = zeros(1,size_glcm_3); % Inverse difference normalized (INN) [3]
out.idmnc = zeros(1,size_glcm_3); % Inverse difference moment normalized [3]

%%% These variables are useful for many of the feature calculations:
glcm_sum  = zeros(size_glcm_3,1);
glcm_mean = zeros(size_glcm_3,1);
glcm_var  = zeros(size_glcm_3,1);

u_x = zeros(size_glcm_3,1);
u_y = zeros(size_glcm_3,1);
s_x = zeros(size_glcm_3,1);
s_y = zeros(size_glcm_3,1);

% checked p_x p_y p_xplusy p_xminusy
p_x = zeros(size_glcm_1,size_glcm_3); % Ng x #glcms[1]
p_y = zeros(size_glcm_2,size_glcm_3); % Ng x #glcms[1]
p_xplusy = zeros((size_glcm_1*2 - 1),size_glcm_3); %[1]
p_xminusy = zeros((size_glcm_1),size_glcm_3); %[1]
% checked hxy hxy1 hxy2 hx hy
hxy  = zeros(size_glcm_3,1);
hxy1 = zeros(size_glcm_3,1);
hx   = zeros(size_glcm_3,1);
hy   = zeros(size_glcm_3,1);
hxy2 = zeros(size_glcm_3,1);

corm = zeros(size_glcm_3,1);
corp = zeros(size_glcm_3,1);

%%% This for loop iterates through each of the GLCM's (from each offset
%%% used) in the inputted 3-D matrix (if more than one offset is used). For
%%% each GLCM, it calculates useful variables:
for k = 1:size_glcm_3
    
    glcm_sum(k) = sum(sum(glcm(:,:,k)));
    glcm(:,:,k) = glcm(:,:,k)./glcm_sum(k); % Normalize each glcm
    glcm_mean(k) = mean2(glcm(:,:,k)); % compute mean after norm
    glcm_var(k)  = (std2(glcm(:,:,k)))^2;
    
    for i = 1:size_glcm_1 %%% Iterates through each gray level within the GLCM
        
        for j = 1:size_glcm_2
            %%% p_x(i) = sum over j of P(i,j), p_y(j) = sum over i of P(i,j)
            p_x(i,k) = p_x(i,k) + glcm(i,j,k);
            p_y(i,k) = p_y(i,k) + glcm(j,i,k); % taking i for j and j for i
            if (ismember((i + j),[2:2*size_glcm_1]))
                p_xplusy((i+j)-1,k) = p_xplusy((i+j)-1,k) + glcm(i,j,k);
            end
            
            %%% This if statement sums to get the P_x-y(z) term, which is equal
            %%% to the sum of all P(i,j) where abs(i - j) = z
            if (ismember(abs(i-j),[0:(size_glcm_1-1)]))
                p_xminusy((abs(i-j))+1,k) = p_xminusy((abs(i-j))+1,k) +...
                    glcm(i,j,k);
            end
        end
    end
    
end

% marginal probabilities are now available [1]
% p_xminusy has +1 in index for matlab (no 0 index)
% computing sum average, sum variance and sum entropy:


%Q    = zeros(size(glcm));

%%% B = repmat(A,M,N) creates a large matrix B consisting of an M-by-N
%%% tiling of copies of A. The size of B is [size(A,1)*M, size(A,2)*N].
%%% The statement repmat(A,N) creates an N-by-N tiling.
%%% i_matrix and j_matrix are matrices with the same size as the GLCM, but
%%% with only i or j values as elements.
i_matrix  = repmat([1:size_glcm_1]',1,size_glcm_2);
j_matrix  = repmat([1:size_glcm_2],size_glcm_1,1);
% i_index = [ 1 1 1 1 1 .... 2 2 2 2 2 ... ]
i_index   = j_matrix(:);
% j_index = [ 1 2 3 4 5 .... 1 2 3 4 5 ... ]
j_index   = i_matrix(:);
xplusy_index = [1:(2*(size_glcm_1)-1)]';
xminusy_index = [0:(size_glcm_1-1)]';
mul_contr = abs(i_matrix - j_matrix).^2;
mul_dissi = abs(i_matrix - j_matrix);
%div_homop = ( 1 + mul_contr); % used from the above two formulae
%div_homom = ( 1 + mul_dissi);

for k = 1:size_glcm_3 % number glcms
    
    out.contrast(k) = sum(sum(mul_contr.*glcm(:,:,k)));
    out.dissi(k) = sum(sum(mul_dissi.*glcm(:,:,k)));
    out.energy(k) = sum(sum(glcm(:,:,k).^2));
    out.entropy(k) = - sum(sum((glcm(:,:,k).*log(glcm(:,:,k) + eps))));
    out.hommatlab(k) = sum(sum((glcm(:,:,k)./( 1 + mul_dissi))));
    out.homop(k) = sum(sum((glcm(:,:,k)./( 1 + mul_contr))));
    % [1] explains sum of squares variance with a mean value;
    % the exact definition for mean has not been provided in
    % the reference: I use the mean of the entire normalized glcm
    out.sosvh(k) = sum(sum(glcm(:,:,k).*((i_matrix - glcm_mean(k)).^2)));
    out.indnc(k) = sum(sum(glcm(:,:,k)./( 1 + (mul_dissi./size_glcm_1) )));
    out.idmnc(k) = sum(sum(glcm(:,:,k)./( 1 + (mul_contr./(size_glcm_1^2)))));
    out.maxpr(k) = max(max(glcm(:,:,k)));
    
    u_x(k)       = sum(sum(i_matrix.*glcm(:,:,k)));
    u_y(k)       = sum(sum(j_matrix.*glcm(:,:,k)));
    % using http://www.fp.ucalgary.ca/mhallbey/glcm_variance.htm for s_x
    % s_y : This solves the difference in value of correlation and might be
    % the right value of standard deviations required
    % According to this website there is a typo in [2] which provides
    % values of variance instead of the standard deviation hence a square
    % root is required as done below:
    s_x(k)  = (sum(sum( ((i_matrix - u_x(k)).^2).*glcm(:,:,k) )))^0.5;
    s_y(k)  = (sum(sum( ((j_matrix - u_y(k)).^2).*glcm(:,:,k) )))^0.5;
    
    corp(k) = sum(sum((i_matrix.*j_matrix.*glcm(:,:,k))));
    corm(k) = sum(sum(((i_matrix - u_x(k)).*(j_matrix - u_y(k)).*glcm(:,:,k))));
    
    out.autocorrelation(k) = corp(k);
    out.corr(k) = (corp(k) - u_x(k)*u_y(k))/(s_x(k)*s_y(k));
    out.correlationmatlab(k) = corm(k) / (s_x(k)*s_y(k));
    
    out.cprom(k) = sum(sum(((i_matrix + j_matrix - u_x(k) - u_y(k)).^4).*...
        glcm(:,:,k)));
    out.cshad(k) = sum(sum(((i_matrix + j_matrix - u_x(k) - u_y(k)).^3).*...
        glcm(:,:,k)));
  
    out.savgh(k) = sum((xplusy_index + 1).*p_xplusy(:,k));
    % the summation for savgh is for i from 2 to 2*Ng hence (i+1)
    out.senth(k) =  - sum(p_xplusy(:,k).*...
        log(p_xplusy(:,k) + eps));
    
    % compute sum variance with the help of sum entropy
    out.svarh(k) = sum((((xplusy_index + 1) - out.senth(k)).^2).*...
        p_xplusy(:,k));
    % the summation for savgh is for i from 2 to 2*Ng hence (i+1)
    
    % compute difference variance, difference entropy,
    out.denth(k) = - sum((p_xminusy(:,k)).*...
        log(p_xminusy(:,k) + eps));
    out.dvarh(k) = sum((xminusy_index.^2).*p_xminusy(:,k));
    
    % compute information measure of correlation(1,2) [1]
    hxy(k) = out.entropy(k);
    glcmk  = glcm(:,:,k)';
    glcmkv = glcmk(:);
    
    hxy1(k) =  - sum(glcmkv.*log(p_x(i_index,k).*p_y(j_index,k) + eps));
    hxy2(k) =  - sum(p_x(i_index,k).*p_y(j_index,k).*...
        log(p_x(i_index,k).*p_y(j_index,k) + eps));
    hx(k) = - sum(p_x(:,k).*log(p_x(:,k) + eps));
    hy(k) = - sum(p_y(:,k).*log(p_y(:,k) + eps));
    
    out.inf1h(k) = ( hxy(k) - hxy1(k) ) / ( max([hx(k),hy(k)]) );
    out.inf2h(k) = ( 1 - exp( -2*( hxy2(k) - hxy(k) ) ) )^0.5;
    
   
end

function [Z, A, Phi] = Zernikmoment(p,n,m)
% -------------------------------------------------------------------------
% Copyright C 2014 Amir Tahmasbi
% Texas A&M University
% amir.tahmasbi@tamu.edu
% http://people.tamu.edu/~amir.tahmasbi/index.html
%
% License Agreement: To acknowledge the use of the code please cite the
%                    following papers:
%
% [1] A. Tahmasbi, F. Saki, S. B. Shokouhi,
%     Classification of Benign and Malignant Masses Based on Zernike Moments,
%     Comput. Biol. Med., vol. 41, no. 8, pp. 726-735, 2011.
%
% [2] F. Saki, A. Tahmasbi, H. Soltanian-Zadeh, S. B. Shokouhi,
%     Fast opposite weight learning rules with application in breast cancer
%     diagnosis, Comput. Biol. Med., vol. 43, no. 1, pp. 32-41, 2013.
%
% -------------------------------------------------------------------------
% Function to find the Zernike moments for an N x N binary ROI
%
% [Z, A, Phi] = Zernikmoment(p,n,m)
% where
%   p = input image matrix (N should be an even number)
%   n = The order of Zernike moment (scalar)
%   m = The repetition number of Zernike moment (scalar)
% and
%   Z = Complex Zernike moment
%   A = Amplitude of the moment
%   Phi = phase (angle) of the mement (in degrees)
%
% Example:
%   1- calculate the Zernike moment (n,m) for an oval shape,
%   2- rotate the oval shape around its centeroid,
%   3- calculate the Zernike moment (n,m) again,
%   4- the amplitude of the moment (A) should be the same for both images
%   5- the phase (Phi) should be equal to the angle of rotation

r = size(p,1); c = size(p,2);
if r ~= c % (Kathleen edit: add extra rows/columns of zeros if matrix not NxN)
    diff = max(r,c) - min(r,c);
    if r > c
        addc = zeros(r,diff);
        p = [p addc];
    elseif c > r
        addr = zeros(diff,c);
        p = [p; addr];
    end
end

N = size(p,1);
if mod(N,2) ~=0 % (Kathleen edit: add extra row and column of zeros if N is not even)
    p = [p zeros(N, 1); zeros(1, N + 1)];
end

x = 1:N; y = x;
[X,Y] = meshgrid(x,y);
R = sqrt((2.*X-N-1).^2+(2.*Y-N-1).^2)/N;
Theta = atan2((N-1-2.*Y+2),(2.*X-N+1-2));
R = (R<=1).*R;
Rad = radialpoly(R,n,m);    % get the radial polynomial

Product = p(x,y).*Rad.*exp(-1i*m*Theta);
Z = sum(Product(:));        % calculate the moments

cnt = nnz(R)+1;             % count the number of pixels inside the unit circle
Z = (n+1)*Z/cnt;            % normalize the amplitude of moments
A = abs(Z);                 % calculate the amplitude of the moment
Phi = angle(Z)*180/pi;      % calculate the phase of the mement (in degrees)


function rad = radialpoly(r,n,m)
% -------------------------------------------------------------------------
% Copyright C 2014 Amir Tahmasbi
% Texas A&M University
% amir.tahmasbi@tamu.edu
% http://people.tamu.edu/~amir.tahmasbi/index.html
%
% License Agreement: To acknowledge the use of the code please cite the
%                    following papers:
%
% [1] A. Tahmasbi, F. Saki, S. B. Shokouhi,
%     Classification of Benign and Malignant Masses Based on Zernike Moments,
%     Comput. Biol. Med., vol. 41, no. 8, pp. 726-735, 2011.
%
% [2] F. Saki, A. Tahmasbi, H. Soltanian-Zadeh, S. B. Shokouhi,
%     Fast opposite weight learning rules with application in breast cancer
%     diagnosis, Comput. Biol. Med., vol. 43, no. 1, pp. 32-41, 2013.
%
% -------------------------------------------------------------------------
% Function to compute Zernike Polynomials:
%
% f = radialpoly(r,n,m)
% where
%   r = radius
%   n = the order of Zernike polynomial
%   m = the repetition of Zernike moment

rad = zeros(size(r));                     % Initilization
for s = 0:(n-abs(m))/2
    c = (-1)^s*factorial(n-s)/(factorial(s)*factorial((n+abs(m))/2-s)*...
        factorial((n-abs(m))/2-s));
    rad = rad + c*r.^(n-2*s);
end


function exportdata(results64,results128,results256,results512,handles)
% Called by export data button functions

% This export creates a few files for each image, so it starts by creating
% a folder to hold all the files. The user is asked to name the folder:
foldername = inputdlg(strcat('A new folder will be created to hold ex', ...
    'ported data. What would you like to name this folder?'));
foldername = foldername{1}; % Extracts string from cell array
mkdir(foldername) % Creates a folder with inputted name
save(strcat(foldername,'\','handles'),'handles') % Save handles structure in folder


%%%% Code below export the histogram/GLCM/ROI 

% for j = 1:length(results) % Go through all analyzed images
%     name = results(j).name{1};
%     name = name(1:(length(name)-4)); % Cuts off '.jpg' from name
%
%     % Create and save histogram
%     h = figure;
%     hist(results(j).roivec)
%     title(name,'Interpreter','none')
%     saveas(h,strcat(foldername,'\',name,'Hist','.fig')) %save as matlab fig
%     saveas(h,strcat(foldername,'\',name,'Hist','.jpg')) %save as jpg image
%     close(h)
%
%     % Create and save plot of glcm
%     % Here, I put a zero in the (1,1) element of the glcm, because otherwise
%     % the (1,1) element is much larger than any other elements, which
%     % affects the scale of the graph and makes other features hard to see.
%     offset = {' 0 deg',' 45 deg',' 90 deg',' 135 deg'};
%     for i = 1:4
%         h = figure;
%         glcm0 = results(j).glcm(:,:,i);
%         glcm0(1,1) = 0;
%         surf(glcm0)
%         title(strcat(name,offset{i}),'Interpreter','none')
%         saveas(h,strcat(foldername,'\',name,offset{i},'GLCM','.fig')) % Save Matlab figure of hist
%         saveas(h,strcat(foldername,'\',name,offset{i},'GLCM','.jpg')) % Save picture of hist
%         close(h)
%     end
%     h = figure; % Create figure to show image so that it can be saved
%     if handles.single.on_off == 1 % Check if dealing with single file or folder
%         imshow(handles.single.imagestore.roi0,[])
%     else
%         imshow(handles.imagestore(j).roi0,[]) % Go through each roi for multiple files
%     end
%     title(strcat(name,' ROI'),'Interpreter','none')
%     saveas(h,strcat(foldername,'\',name,'ROI','.fig')) % Save roi image as Matlab fig
%     saveas(h,strcat(foldername,'\',name,'ROI','.jpg')) % Save roi image as picture
%     close(h)
% end

% The rest of this function exports the data into a '.csv' file:

%list = {'results16' 'results32' 'results64' 'results128'};
%results = [results16, results32, results64, results128];

resultstble(foldername,'results64',results64)
resultstble(foldername,'results128',results128)
resultstble(foldername,'results256',results256)
resultstble(foldername,'results512',results512)

%%%% Creates the csv file with the extracted texture features 
function resultstble(foldername,list,fuckthis)
    
    fid = fopen(strcat(foldername,'\',foldername,list,'results','.csv'),'wt');
    names = fieldnames(fuckthis)';
    
    for k = 1:length(names)
        fprintf(fid,'%s,',names{1,k}); % Print fieldnames as column labels
    end
    
    fprintf(fid, '\n');
    
    for j = 1:length(fuckthis) % Print each image's data in a different row
        name = fuckthis(j).name{1};
        name = name(1:(length(name)-4));
        fprintf(fid, '%s,',name);
        fprintf(fid, '%d,',fuckthis(j).numpixels);
        fprintf(fid, '%d,',fuckthis(j).mean);
        fprintf(fid, '%d,',fuckthis(j).std);
        fprintf(fid, '%d,',fuckthis(j).median);
        fprintf(fid, '%d,',fuckthis(j).max);
        fprintf(fid, '%d,',fuckthis(j).min);
        fprintf(fid, '%d,',fuckthis(j).range);
        fprintf(fid, '%d,',fuckthis(j).variance);
        fprintf(fid, '%d,',fuckthis(j).skew);
        fprintf(fid, '%d,',fuckthis(j).kurtosis);
        fprintf(fid, '%d,',fuckthis(j).zernike);
        fprintf(fid, '--,');
        fprintf(fid, '--,');
        fprintf(fid, '--,');
        fprintf(fid, '--,');
        fprintf(fid, '--,');
        fprintf(fid, '%d,',mean(fuckthis(j).autocorrelation));
        fprintf(fid, '%d,',mean(fuckthis(j).contrast));
        fprintf(fid, '%d,',mean(fuckthis(j).correlationmatlab));
        fprintf(fid, '%d,',mean(fuckthis(j).corr));
        fprintf(fid, '%d,',mean(fuckthis(j).cprom));
        fprintf(fid, '%d,',mean(fuckthis(j).cshad));
        fprintf(fid, '%d,',mean(fuckthis(j).dissi));
        fprintf(fid, '%d,',mean(fuckthis(j).energy));
        fprintf(fid, '%d,',mean(fuckthis(j).entropy));
        fprintf(fid, '%d,',mean(fuckthis(j).hommatlab));
        fprintf(fid, '%d,',mean(fuckthis(j).homop));
        fprintf(fid, '%d,',mean(fuckthis(j).maxpr));
        fprintf(fid, '%d,',mean(fuckthis(j).sosvh));
        fprintf(fid, '%d,',mean(fuckthis(j).savgh));
        fprintf(fid, '%d,',mean(fuckthis(j).svarh));
        fprintf(fid, '%d,',mean(fuckthis(j).senth));
        fprintf(fid, '%d,',mean(fuckthis(j).dvarh));
        fprintf(fid, '%d,',mean(fuckthis(j).denth));
        fprintf(fid, '%d,',mean(fuckthis(j).inf1h));
        fprintf(fid, '%d,',mean(fuckthis(j).inf2h));
        fprintf(fid, '%d,',mean(fuckthis(j).indnc));
        fprintf(fid, '%d\n',mean(fuckthis(j).idmnc));
        fprintf(fid, '%d\n');
        
    end
    fclose(fid);
    
    
    
% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(~, ~, ~)
% INSTRUCTIONS
% This creates a new figure window that displays instructions on how to use
% the Texture Analysis GUI

% Create the figure
inst = figure();

% Create array of strings containing instructions:
m = {};
s = {};
s{1} = 'Instructions for using the Texture Analysis GUI:';
s{2} = '';
s{3} = 'Analyzing a Single Image:';
s{4} = '';
s{5} = '1. Click "Select Single Image File" button, and select the image file you wish to analyze.';
s{6} = '2. The image will be displayed on the GUI, and the cursor will appear as a small cross. Click a point on ';
s{7} = 'the image to start drawing a polygon around the ROI you wish to analyze.';
s{8} = '3. After you draw the ROI boundaries, you can click and drag them to adjust the vertices.';
s{9} = '4. Once you are happy with the ROI, double click on the image.';
s{10} = '5. Now that you have successfully selected the ROI, you can click "Analyze Single ROI" to perform';
s{11} = 'texture analysis on the ROI. The results will be printed in the Command Window.';
s{12} = '6. You can now click "Export Results" to create a folder containing files of the exported analysis';
s{13} = 'results, and/or select "View Results" to view the results in a Matlab figure window.';
s{14} = '';
s{15} = 'Analyzing Multiple Images:';
s{16} = '';
s{17} = '1. Click "Choose Folder of Image Files" button, and select the folder containing all the image files you ';
s{18} = 'wish to analyze.';
s{19} = '2. Each of the images contained in the folder will be displayed (one at a time) on the GUI, and the cursor';
s{20} = 'will appear as a small cross. Click a point on the image to start drawing a polygon around the ROI you ';
s{21} = 'wish to analyze.';
s{22} = '3. After you draw the ROI boundaries, you can click and drag them to adjust the vertices.';
s{23} = '4. Once you are happy with the ROI, double click on the image and press "Next File" to proceed.';
s{24} = '5. When the next image from the folder appears, repeat steps 2-4 to draw ROIs for all the images.';
s{25} = '6. Once you have selected an ROI on all the images, if you click "Next File", a dialogue box will alert you';
s{26} = 'that you have reached the last image.';
s{27} = '7. Now that you have successfully selected the ROIs, you can click "Analyze All ROIs" to perform';
s{28} = 'texture analysis on the ROIs. The results will be printed in the Command Window.';
s{29} = '8. You can now click "Export Results" to create a folder containing files of the exported analysis';
s{30} = 'results, and/or select "View Results" to view the results in a Matlab figure window.';
s{31} = '';

% Iterate over array to properly place each string of text in figure:
for i = 1:length(s)
    m{i} = uicontrol('style','text'); % Create a uicontrol of type "text"
    set(m{i},'String',s{i}) % Choose string of text
    set(m{i},'Units','characters')
    set(m{i},'Position',[5 32-i 100 1]) % Set position of text
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(~, ~, ~)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
