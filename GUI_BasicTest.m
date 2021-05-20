%%  Background code for "GUI_BasicTest.fig"
%   
%   This code runs data logging, the currently available options are:
%   - Trigger - The logger will start recording information when an
%   analogue input signal passes an initial threshold.
%   - Oscilloscope - A plot will show that shows the requested signals
%   live.
%   - Continuous logging - The logger will run until the stop logging
%   button is pushed. This mode is intended for use with the Oscilloscope.
%   No output logfile is created.
%
%
%   DW - 15/05/21 - Modified with trigger
%   DW - 20/05/21 - Tested with PD
%%  Main
function varargout = GUI_BasicTest(varargin)
% GUI_BASICTEST MATLAB code for GUI_BasicTest.fig
%      GUI_BASICTEST, by itself, creates a new GUI_BASICTEST or raises the existing
%      singleton*.
%
%      H = GUI_BASICTEST returns the handle to a new GUI_BASICTEST or the handle to
%      the existing singleton*.
%
%      GUI_BASICTEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_BASICTEST.M with the given input arguments.
%
%      GUI_BASICTEST('Property','Value',...) creates a new GUI_BASICTEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_BasicTest_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_BasicTest_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_BasicTest

% Last Modified by GUIDE v2.5 19-May-2021 15:26:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EpiC_DATAQ_Trig_OpeningFcn, ...
                   'gui_OutputFcn',  @EpiC_DATAQ_Trig_OutputFcn, ...
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
end
% End initialization code - DO NOT EDIT



% --- Executes just before GUI_BasicTest is made visible.
function EpiC_DATAQ_Trig_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_BasicTest (see VARARGIN)

% Choose default command line output for GUI_BasicTest
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_BasicTest wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end



% --- Outputs from this function are returned to the command line.
function varargout = EpiC_DATAQ_Trig_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end



%%  Sample Duration
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SampleRate = str2double(get(hObject, 'String'));
if isnan(SampleRate)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



%%  Sample rate
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SampleRate = str2double(get(hObject, 'String'));
if isnan(SampleRate)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



%%  Output Filename
function edit3_Callback(hObject, eventdata, handles)

% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



%% Continuous recording
%  Stop running push button
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles, s)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Stop data acquisition switch 
global StpRecord
StpRecord = 1;
end

% Continous running check-box
% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
end


%% Trigger channel
function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


%%  Trigger active
% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
end



%%  Trigger value
function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
end

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



%%  Oscilloscope
% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6

end



%%  Run Datalogger
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(~, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%   Initiation
daqreset

%   Extract input strings (duration, rate, filename)
for i = 1:5
    DAQInput{i} = ExtractString(handles, i);
end

%   Extract active channel and calibration data
ChanCal = handles.uitable1.Data;
ActChn = cell2mat(ChanCal(:,1));
ActiveChannels = find(ActChn)-1;    %

CalFct = cell2mat(ChanCal(:,2));

optionPlot = handles.checkbox6.Value;
optionTrig = handles.checkbox4.Value;
optionCont = handles.checkbox3.Value;

%   Trigger options
if optionTrig == 1
    Trig.Chan = int8(str2double(DAQInput{4})+1); % Trigger channel (+1 to make it a Matlab index)
    Trig.Val = str2double(DAQInput{5});                                 % Trigger value     % TODO - make this a textbox and checkbox
    
end

%   Set Session Properties - not using default values.
sampleFreq = str2double(DAQInput{2});     %   Session frequency
captureLength = str2double(DAQInput{1});    %   Capture X seconds after trigger, must be multiple of myBufferRefillPeriod
bufferRefillPeriod = 0.1;                 %   in seconds
numPretriggerBuffers = 1;                 %   How many refill periods before trigger?
numChannels = int8(sum(ActChn));          %   Number of active chane;s

%   The samples are stored in a circular buffer, which is a
%   computationally efficient way of storing the samples
bufferSize = ceil(bufferRefillPeriod * sampleFreq);
numBuffers = numPretriggerBuffers + (captureLength/bufferRefillPeriod) + 1; % plus 1 to account for mid refill trigger
currentBuffer = 1;

%   Create log file
TimeStamp = clock;
logFName = sprintf('%s_Date%0.2i%0.2i%0.2i_Time%0.2i%0.2i%0.2i%s', DAQInput{3}(1:end-4), TimeStamp(1:5) , ...
                round(TimeStamp(end)), DAQInput{3}(end-3:end));
fid1 = CreateLog(logFName, DAQInput, TimeStamp, ActiveChannels);

%   Log data, format
FrmtPrp = cell(1,length(ActiveChannels));
FrmtPrp(:) = {'%f\t'};
FrmtPrp(end+1) = {'\n'}; % OLD '\r\n'
Frmt = strcat(FrmtPrp{:});



%%  Create Data Acquisition Session
%   Session
s = daq.createSession('ni');
s.Rate = sampleFreq;             
s.IsContinuous = true;

%   Add channels
ch = addAnalogInputChannel(s,'Dev1',ActiveChannels,'Voltage');

%   Change recording to single ended
for i = 1:length(find(ActChn))
    ch(i).TerminalConfig = 'SingleEnded';
end

if optionPlot == 1
    %   Create log-chart
    Axe1 = RecordFig(ActiveChannels);
end



%%  Logging options 
%   Define the set of buffers
circularBuffers = zeros(bufferSize, numChannels, numBuffers);
circularTimeBuffer = zeros(bufferSize, 1, numBuffers);

%   Logic for logging
dataBeingLogged = false;

%   Add listeners
if optionPlot == 1
    if exist('Trig','var')
        LstPlt = addlistener(s,'DataAvailable', @(src, event)plotData(src, event, Axe1, Trig.Val)); % Live plotting
    else
        LstPlt = addlistener(s,'DataAvailable', @(src, event)plotData(src, event, Axe1)); % Live plotting
        
    end
end

if optionTrig == 1
    LstRec = addlistener(s,'DataAvailable',@(src, event)logData(src, event)); %   Trigger / record data
end

if optionCont == 1
    LstRec = addlistener(s,'DataAvailable',@(src, event)continuousData(src, event)); %   Trigger / record data    
end

%   Initiate logging variables
trigInd = 0;
global exitWait
exitWait = 0; % Logic for exiting the waiting loop

%   Variable for stop button
global StpRecord
StpRecord = 0;



%%  Check inputs 
if exist('Trig','var')     
    %   Check trigger channel entered correctly
    if any([isnan(Trig.Chan), gt(Trig.Chan,7), lt(Trig.Chan,0)])
        error('Trigger channel must be an integer between 0 and 7')
    end 
end

if and(optionTrig, optionCont)
    error('Cannot have trigger and continuous logging active')
end



%%  Logging
%   Start logging
startBackground(s);

%   Pause matlab execution while recording (until exitWait updates)
T = timer('TimerFcn',@(~,~)1+1,'StartDelay',0.01);
% The TimerFcn is a dummy function that does nothing

while exitWait == 0
    start(T)
    wait(T)  
end



%%  Post-process

if not(StpRecord == 1) % If exist recording don't make plot
    
    %   Frequency plot
    if not(optionCont == 1)
        %   Plot fft of logged data
        RecData = ExtractData(logFName);
        RecData.ProcData = RecData.Data-mean(RecData.Data);

        figure
        AxeFft = axes;
        AxeFft.Title.String = ['Recording Channels -', sprintf('%i,', ActiveChannels(:))];
        AxeFft.XLabel.String = 'Frequency [Hz]';
        AxeFft.YLabel.String = 'Power Spectral Density';
        hold on

        for i = 1:length(ActiveChannels)
            [Pxy,W] = cpsd(RecData.ProcData(:,i),RecData.ProcData(:,i),[],0,...
                RecData.Freq/0.0001, RecData.Freq);
            PltFFT(i) = plot(W,Pxy);
        end
        % AxeFft.Title.String = 'Fft of data-logger recording';
        AxeFft.XLim = [0 20];

        %   Update FFT labels to match time-series plot
        ColorBank = {[0 0 0],[1 0.1 0.1],[0.2 0.2 0.2],[1 0.3 0.3],[0.4 0.4 0.4],...
            [1 0.5 0.5],[0.6 0.6 0.6],[1 0.7 0.7]};
        LineBank = {'-', '--', '-.', ':', '-', '--', '-.', ':'};
        %   Cycle through lines 
        for i = 1:length(PltFFT)
            PltFFT(i).Color = ColorBank{i};
            PltFFT(i).LineStyle = LineBank{i};
        end

        for i = 1:length(PltFFT)    
            LgdText{i} = sprintf('Channel %i',ActiveChannels(i));
        end
        FFTLogLgd = legend(AxeFft, LgdText);
        FFTLogLgd.Location = 'eastoutside';
    end 

end

%   Tidy workspae
clear
fclose('all');
daqreset



%%   Nested function
%  Log data
function logData(src, event) % Log data in file

% Check stopped recording
if StpRecord == 1
    %   Stop logging session
    s.stop();
    s.release();
    
    %   Signal to exit waiting loop & return
    exitWait = 1;
    return
end
    
% Read new data
newData = event.Data;
newTime = event.TimeStamps;

% Refill the buffer and throw out data beyond pretrigger time (FIFO)
circularBuffers = cat(3, circularBuffers(:,:,2:end), newData);
circularTimeBuffer = cat(3, circularTimeBuffer(:,:,2:end), newTime);


%   Trigger
trig = newData(:,Trig.Chan) > Trig.Val; % triggers when digital signal is high


% Find the first data point that fits the trigger condition
if any(trig) && (dataBeingLogged == false)
    trigInd = find(trig, 1, 'first');
    dataBeingLogged = true;
end


if dataBeingLogged == true
    if currentBuffer == (captureLength/bufferRefillPeriod) + 1 % +1 to account for the padded buffer (line 27)
        % Reorganize data once trigger condition met
        circularBuffers = permute(circularBuffers, [1 3 2]); % rearrange 3D matrix for proper reshape
        circularBuffers = reshape(circularBuffers, numBuffers*bufferSize, numChannels);
        circularTimeBuffer = reshape(circularTimeBuffer, numBuffers*bufferSize, 1);
        
        % Clean up HW resources
        s.stop();
        s.release();
        delete(LstRec);     
        if exist('LstPlt')  %   Delete plot listener if it exists
            delete(LstPlt)
        end
        
        % Extract the pretrigger length and the capture length from
        % the buffer.
        numBufferSamples = numPretriggerBuffers*bufferSize; 
        fullCapTrigInd = numBufferSamples + trigInd;
        trigTime = circularTimeBuffer(fullCapTrigInd);
        numCapSamples = captureLength * sampleFreq;
        reqCapture = circularBuffers((fullCapTrigInd-numBufferSamples):(fullCapTrigInd+numCapSamples),:);
        reqCaptureTime = circularTimeBuffer((fullCapTrigInd-numBufferSamples):(fullCapTrigInd+numCapSamples),:) - trigTime;

        % 	Output    
        RecordData(reqCapture, fid1, Frmt)  %   Write data to logfile
        outputData = reqCapture; %   Save as .mat
        save([logFName(1:end-3) 'mat'], 'outputData', 'sampleFreq')  
        
        % Plot
        
        %   Line options
        ColorBank = {[0 0 0],[0.6 0.6 0.6],[1 0.1 0.1],[0.2 0.2 0.2],[1 0.3 0.3],[0.4 0.4 0.4],...
            [1 0.5 0.5],[1 0.7 0.7]};
        LineBank = {'-', '--', '-.', ':', '-', '--', '-.', ':'};
        for indFig = 1:size(reqCapture,2) 
            LgdTextEnd{indFig} = sprintf('Channel %i',ActiveChannels(indFig));
        end
        
        figure
        axeEnd = axes;
        hold on

        for indFig = 1:size(reqCapture,2)
            pltDatEnd{indFig} = plot(reqCaptureTime, reqCapture(:,indFig));
            pltDatEnd{indFig}.Color = ColorBank{indFig};
            pltDatEnd{indFig}.LineStyle = LineBank{indFig};
            pltDatEnd{indFig}.DisplayName = LgdTextEnd{indFig};
            
        end
        plotbrowser('on')

        %   Signal to exit waiting loop
        exitWait = 1;
    end
    currentBuffer = currentBuffer + 1;
    
end  


function continuousData(src, event)
%   Listener function for continuous data logging
%
if StpRecord == 1
    %   Stop logging session
    s.stop();
    s.release();
    
    %   Delete listeners 
    delete(LstCont);     
    if exist('LstPlt')  %   Delete plot listener if it exists
        delete(LstPlt)
    end       

    %   Signal to exit waiting loop
    exitWait = 1;
end

end
    


end
end


%%  Local functions
function CustString = ExtractString(handles, i)
%   Function to extract data from strings in the input textboxes
%
FName = char(join(['edit' string(i)], ''));
a = getfield(handles, FName);
CustString = a.String;

end


function Axe1 = RecordFig(ActiveChannels)
%   Function to setup the oscilloscope figure
%
figure;
Axe1 = axes;
Axe1.Title.String = ['Recording Channels -', sprintf('%i,', ActiveChannels(:))];
Axe1.XLabel.String = 'Time [s]';
Axe1.YLabel.String = 'Amplitude';
hold on

end


function [fid1] = CreateLog(FName, DAQInput, TimeStamp, ActiveChannels)
%   Function to setup the logging output text file
%
fid1 = fopen(FName,'w');

% Write log-file header
WriteHeader(fid1, FName, DAQInput, TimeStamp, ActiveChannels)

end


function WriteHeader(fid, FName, DAQInput, TimeStamp, ActiveChannels)
%   Function to write a header for the logging text file
%
% Write header to data log file
A{1} = '---------- DAQ NI USB-6003 Output File ----------';
A{2} = sprintf('%s - Filename', FName);
A{3} = sprintf('%0.2i/%0.2i/%0.2i - File opened date', TimeStamp(1:3));
A{4} = sprintf('%0.2i:%0.2i:%0.4g - File opened time', TimeStamp(4:5), TimeStamp(6));
A{5} = sprintf('%f - Sampling frequency [Hz]', str2double(DAQInput{2}));
A{6} = [sprintf('%i, ', ActiveChannels(:)), ' - DAQ active channels' ];
A{7} = '-------------------------------------------------';

% Write header to the open log-file
fprintf(fid, '%s \n', A{:});

end


function plotData(src, event, AxeNme, TrigVal)
%   Listener function for oscilloscope plotting
%

% Move axes with recording - only show reecent 5 seconds
AxeNme.XLim = [event.TimeStamps(end)-5 event.TimeStamps(end)];

%   Plot active channels
for i = 1:size(event.Data,2)
    PltDat{i} = plot(AxeNme, event.TimeStamps, event.Data(:,i));
    PltDat{i}.HandleVisibility = 'off';    % To prevent 
end

%   Line options
ColorBank = {[0 0 0],[0.6 0.6 0.6],[1 0.1 0.1],[0.2 0.2 0.2],[1 0.3 0.3],[0.4 0.4 0.4],...
    [1 0.5 0.5],[1 0.7 0.7]};
LineBank = {'-', '--', '-.', ':', '-', '--', '-.', ':'};

% Cycle through lines 
for i = 1:length(PltDat)
    PltDat{i}.Color = ColorBank{i};
    PltDat{i}.LineStyle = LineBank{i};
end

%   If 1st loop plot trigger and add legend 
if event.TimeStamps(1) == 0
    if exist('TrigVal', 'var')
        PltTrig = plot([0 600], [TrigVal TrigVal]);
        %   Assume logging will be less than 10 minutes
        PltTrig.Color = [1 0 0];
    end

    %   Turn line handles on for adding legend
    for i = 1:length(PltDat)
        PltDat{i}.HandleVisibility = 'on';
    end
    
    %   Add title
    TitleString = strsplit(AxeNme.Title.String, {'-',','});
    TitleActChn = TitleString(2:end-1);
    ActChn = cellfun(@(x)str2num(x), TitleActChn);

    %   Add legend elements
    for i = 1:length(PltDat)
        LgdText{i} = sprintf('Channel %i',ActChn(i));
    end

    if exist('TrigVal', 'var')
        LgdText{i+1} = 'Trigger Value';
    end

    %   Plot legend
    LogLgd = legend(AxeNme, LgdText);
    LogLgd.Location = 'eastoutside';
end


end
    

function RecordData(data, fid1, formatSpec) % Log data in file
% Log data in format specified by formatSpec in open fid1 with event.Data

fprintf(fid1,formatSpec,data');

end    




function RecData = ExtractData(FileName)
%%  Extract data from DAQ USB-6003 data logging session
%   Assume using "GUI_BasicTest" matlab scripts
%
%   Input:
%   - FileName - string with logging session filename, include .txt
%   extentsion
%
%   Output:
%   - RecData - Struct containing data extracted form data logger session
%
%   DW - 04/04/19 - Created
%%  Record data
%   Check input variable
% if ~isstring(FileName) % Input might be character vector
%     error('Filename must be entered as a string.')
% end

if ~exist(FileName,'file')
   error('Filename: %s, does not exist on Matlab path.', FileName)
end

%   Open data logger file
FileID = fopen(FileName,'r');

% Process header
fseek(FileID, 51, 'bof'); % Skip the header first line
FileName = textscan(fgetl(FileID), '%s %s'); % Read filename
textscan(fgetl(FileID), '%s %s'); % Read date
textscan(fgetl(FileID), '%s %s'); % Read time
SampleFreq = textscan(fgetl(FileID), '%f %s'); % Read sampling frequency 

% Evaluate active channels
ChnExt = textscan(fgetl(FileID), '%s', 'Delimiter',{','});
ChnExt{1}(end) = [];
Chn = cellfun(@(x)str2double(x), ChnExt{1});
fseek(FileID, 51, 'cof'); % Skip header end line

%  Read logged data
DataStrFmt = '%f\t';
for i = 1:length(Chn-1)
    DataStrFmt = [DataStrFmt '%f\t'];
end
A = fscanf(FileID, [DataStrFmt '\n']); 

for i = 1:length(Chn)
   Data(:,i) = A(i:length(Chn):end); 
end
A = [];
fclose('all');

% Format output data struct
RecData.FileName = FileName{1}{1};
RecData.Freq = SampleFreq{1};
RecData.Chn = Chn;
RecData.Data = Data;



end
