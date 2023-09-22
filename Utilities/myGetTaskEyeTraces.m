% getTaskEyeTraces.m
%
%        $Id:$ 
%      usage: e = getTaskEyeTraces(stimfileName,<taskNum=1>,<phaseNum=1>,<dispFig=0>,<dataPad=3>,<removeBlink=1>)
%         by: justin gardner
%       date: 06/23/10
%    purpose: get the eye traces for a task. Pass in a stimfile name and the
%             taskNum and phaseNum you want. The eye traces will be returned
%             in e.eye with units of degrees of visual angle.
%             Data will be exracted by up to dataPad(=3) seconds
%             beyond the end of the trial so that you
%             can extract data around the last event
%             in the trial if it's close to the trial end.
%
%             If you want to nan out blink data, set removeBlink to 1
%             e = getTaskEyeTraces(('101203_stim03','removeBlink=1');
%             This will remove blink data += 5 ms. To set a larger range
%             like say 100ms around the blinks to remove:
%             e = getTaskEyeTraces(('101203_stim03','removeBlink=0.1');
%             To remove say 5 samples from either side of the blink
%             e = getTaskEyeTraces(('101203_stim03','removeBlink=5');
%
%             % if you want to display a figure with the eye positon
%             e = getTaskEyeTraces(('101203_stim03','dispFig=1');
%
%             To get the eye traces relative to a segment other than 1
%             set the alternative field segNum equal to that segment
%             e = getTaskEyeTraces('100616_stim01','segNum=2');
%
%    e = getTaskEyeTraces('100616_stim01','taskNum=1','phaseNum=1','dataPad=3');
%
function e = myGetTaskEyeTraces(stimfileName,varargin)

getBothEyes=0;

% it would be nice if this was fully compatible with getTaskParameters, which
% would require it to be able to take the myscreen & task structs... and also
% simply provide the 'eye' struct when there is data and not without. this seems
% too complicated and not in line with the simplicity of the base task funs

e = [];

% check arguments
if nargin == 0
  help getTaskEyeTraces
  return
end

taskNum=[];phaseNum=[];dispFig=[];dataPad=[];removeBlink=[];
if exist('getArgs') == 2
  getArgs(varargin,{'taskNum=1','phaseNum=1','dispFig=0','dataPad=3','removeBlink=0','segNum=1'});
else
  disp(sprintf('(getTaskEyeTraces) To run this program you need functions from the mrTools distribution. \nSee here: http://gru.brain.riken.jp/doku.php/mgl/gettingStarted#initial_setup'));
  return
end

% if we are passed in the name of a file then assume it is a stimfile and load it through
% getTaskParameters
if isstr(stimfileName)
  [e stimfile] = getTaskParameters(stimfileName);
  if isempty(e)
    return
  end
  % make sure e/task are the correct form of cell array
  e = cellArray(e);
  stimfile.task = cellArray(stimfile.task,2);
else
  disp(sprintf('(getTaskEyeTraces) Must pass in stimfile name'));
  return
end

% get the correct task and phase
if taskNum > length(e)
  disp(sprintf('(getTaskEyeTraces) taskNum=%i out of range for this stimfile (1:%i)',taskNum,length(e)));
  return
end
if phaseNum > length(e{taskNum})
  disp(sprintf('(getTaskEyeTraces) phaseNum=%i out of range for this task (1:%i)',phaseNum,length(e{taskNum})));
  return
end

e = e{taskNum}(phaseNum);

% keep stimfile
e.stimfile = stimfile;
e.stimfile.taskNum = taskNum;
e.stimfile.phaseNum = phaseNum;

% check for taskID
if ~isfield(stimfile.task{taskNum}{phaseNum},'taskID')
  if ~isfield(stimfile.task{taskNum}{phaseNum},'collectEyeData')
    disp(sprintf('(getTaskEyeTraces) **** No taskID field found in task. This stimfile was probably generated with an older version of mgl/task. You need to update your mgl code. ****'));
    return
  else
    % for mglEyelink V1 messages only one task could collect data
    if (stimfile.task{taskNum}{phaseNum}.collectEyeData == 1)
      taskID = 0;
      phaseNum = 1; % start with phase 1 in the recorded mgl messages 
    else
      taskID = NaN;
    end
  end
else
  taskID = stimfile.task{taskNum}{phaseNum}.taskID;
end

% check eye tracker type
eyeTrackerType = stimfile.myscreen.eyeTrackerType;

% if isequal (eyeTrackerType,'None')
%     eyeTrackerType = 'Eyelink';
% end

if ~isequal(eyeTrackerType,'Eyelink')
  disp(sprintf('(getTaskEyeTraces) Loading of eye tracker type %s not implemented',eyeTrackerType));
  return
else
  eyeTrackerFileExt = 'edf';
end

% get the filename
if isfield(stimfile.myscreen.eyetracker,'datafilename')
  eyeTrackerFilename = stimfile.myscreen.eyetracker.datafilename;
else
  disp(sprintf('(getTaskEyeTraces) No associated eyetracker datafilename found in myscreen'));
  return
end

% check for file, should be in myscreen directory
origFilename = eyeTrackerFilename;
eyeTrackerFilename = fullfile(stimfile.stimfilePath,sprintf('%s.%s',eyeTrackerFilename,eyeTrackerFileExt));
if ~isfile(eyeTrackerFilename)
    disp(sprintf('(getTaskEyeTraces) Could not find eye tracker file %s',eyeTrackerFilename));
    % if not there, it may be in the same folder as the stimfile
    [filepath,name,ext] = fileparts(stimfileName);
    eyeTrackerFilename = fullfile(filepath,sprintf('%s.%s',origFilename,eyeTrackerFileExt));
    if ~isfile(eyeTrackerFilename)
        disp(sprintf('(getTaskEyeTraces) Could not find eye tracker file %s',eyeTrackerFilename));
        e = []; % June
        return
    end
end


% replace tilde
if exist('mlrReplaceTilde') == 2
  eyeTrackerFilename = mlrReplaceTilde(eyeTrackerFilename);
else
  if ~isempty(findstr('~',eyeTrackerFilename))
    disp(sprintf('(getTaskEyeTraces) The ~ in filename %s may not be parsed correctly',eyeTrackerFilename));
  end
end

% load the file
disppercent(-inf,sprintf('(getTaskEyeTraces) Opening edf file %s',eyeTrackerFilename));
sprintf('\n');
edf = myMglEyelinkEDFRead(eyeTrackerFilename);

% return empty e for invalid files -June, 02.16
if isempty(edf)
    e = [];
    return
end


% if ~getBothEyes
%     edf = rmfield(edf,'gazeBoth');
% end
edf.samplerate = edf.samplerate(1);
disppercent(inf);
% if isempty(edf),return,end

% get all the messages that match our taskID
% get the number of trials
edf.nTrials = max(edf.mgl.trialNum(find(edf.mgl.taskID == taskID)));

if removeBlink
    % blink window, extra padding (50ms)
    if removeBlink==1 % == true because we can't pass logical and test with getArgs, use 1-eps for 1s
        blinkWindow = ceil(0.005*edf.samplerate); % default
    elseif removeBlink > 1 % assume in samples
        blinkWindow = ceil(removeBlink);
    elseif removeBlink < 1 % assume seconds
        blinkWindow = ceil(removeBlink*edf.samplerate);        
    end
    for Bn = 1:numel(edf.blinks.startTime)
        blinks = (edf.gaze.time >= edf.blinks.startTime(Bn)-blinkWindow & ...
                         edf.gaze.time <= edf.blinks.endTime(Bn)+blinkWindow);
        edf.gaze.x(blinks) = NaN;
        edf.gaze.y(blinks) = NaN;
        edf.gaze.pupil(blinks) = NaN;
        
        if getBothEyes
            for eyeNum=1:2 %left then right
                blinks = (edf.gazeBoth{eyeNum}.time >= edf.blinks.startTime(Bn)-blinkWindow & ...
                    edf.gazeBoth{eyeNum}.time <= edf.blinks.endTime(Bn)+blinkWindow);
                edf.gazeBoth{eyeNum}.x(blinks) = NaN;
                edf.gazeBoth{eyeNum}.y(blinks) = NaN;
                edf.gazeBoth{eyeNum}.pupil(blinks) = NaN;
            end
        end
        
    end
end



% now process each trial
disppercent(-inf,sprintf('(getTaskEyeTraces) Extracting trial by trial data for %i trials',edf.nTrials));
% get the start times
for i = 1:max(edf.mgl.trialNum(edf.mgl.taskID == taskID)) 
  % get start time
  % find the segment 0 message
  %%% I think this should be segment 1--at least in my code segment 1 == seg1
  %%% and that seems to be what updateTask writes out. The seg==0 often includes
  %%% deadtime related to waiting for backtics, user start, etc
  segmentOneTime = edf.mgl.time((edf.mgl.taskID == taskID) &  ...
                                   (edf.mgl.trialNum==i) &  ...
                                   (edf.mgl.segmentNum==1));
  segNumTime = edf.mgl.time((edf.mgl.taskID == taskID) &  ...
			    (edf.mgl.trialNum==i) &  ...
			    (edf.mgl.segmentNum==segNum));
            
  if isempty(segNumTime) && isempty(segmentOneTime)%June
      e = [];%June
      return %June
  end%June
  
  % call this the startTime %what to do for segNumTime = []?
  if ~isempty(segNumTime)
    startTime(i) = segNumTime;
  else
    startTime(i) = nan;
  end
  % end time is the start of the next trial
  if i > 1
    if ~isempty(segmentOneTime)
      endTime(i-1) = segmentOneTime;
    else
      endTime(i-1) = nan;
    end
  else%June
      endTime = [];%June
  end
end

% make the end time of the last trial, at most as long as the longest trial so far
 
if ~isempty(endTime)
  maxTrialLen = max(endTime-startTime(1:end-1))+1;
  endTime(end+1) = min(max(edf.gaze.time),startTime(end)+maxTrialLen-1);
% endTime(end+1) = max(edf.gaze.time); %changed this because endTime was too short. ZR
else
  % if we have only one trial, then use till the end of the data (what if
  % the startTime is NaN? It returns NaN for maxTrialLen
  endTime(end+1) = max(edf.gaze.time);
  maxTrialLen = endTime-startTime+1;
end
maxTrialLen = maxTrialLen+dataPad;

% now get time between samples (isn't the sample rate available?)
% timeBetweenSamples = median(diff(edf.gaze.time));
% get the time between samples in milliseconds
timeBetweenSamples = (1/edf.samplerate)*1000;

% figure out how large to make data array
e.eye.xPos = nan(edf.nTrials,ceil(maxTrialLen/timeBetweenSamples));
e.eye.yPos = nan(edf.nTrials,ceil(maxTrialLen/timeBetweenSamples));
e.eye.pupil = nan(edf.nTrials,ceil(maxTrialLen/timeBetweenSamples));

% put in time in seconds
e.eye.time = (0:(size(e.eye.xPos,2)-1))/edf.samplerate;



% go through each trial and populate traces
warning('off','MATLAB:interp1:NaNinY');

for iTrial = 1:edf.nTrials
  disppercent(iTrial/edf.nTrials);
  % the times for this trial
  thisTrialTimes = startTime(iTrial):timeBetweenSamples:endTime(iTrial);
  % get data for xPos, yPos and pupil traces form edf data 
  e.eye.xPos(iTrial,1:length(thisTrialTimes)) = interp1(edf.gaze.time,edf.gaze.x,thisTrialTimes,'linear',nan);
  e.eye.yPos(iTrial,1:length(thisTrialTimes)) = interp1(edf.gaze.time,edf.gaze.y,thisTrialTimes,'linear',nan);
  e.eye.pupil(iTrial,1:length(thisTrialTimes)) = interp1(edf.gaze.time,edf.gaze.pupil,thisTrialTimes,'linear',nan);
end

warning('on','MATLAB:interp1:NaNinY');
disppercent(inf);

% convert to device coordinates
w = stimfile.myscreen.screenWidth;
h = stimfile.myscreen.screenHeight;
xPix2Deg = stimfile.myscreen.imageWidth/w;
yPix2Deg = stimfile.myscreen.imageHeight/h;

hDir = 1;vDir = 1;
if isfield(stimfile.myscreen,'flipHV') && (length(stimfile.myscreen.flipHV) >=2 )
  if stimfile.myscreen.flipHV(1) hDir = -1;end
  if stimfile.myscreen.flipHV(2) vDir = -1;end
end
e.eye.xPos = hDir * ((e.eye.xPos-(w/2))*xPix2Deg);
e.eye.yPos = vDir * (((h/2)-e.eye.yPos)*yPix2Deg);


if getBothEyes
    for eyeNum=1:2
        % figure out how large to make data array
        e.eyes{eyeNum}.xPos = nan(edf.nTrials,ceil(maxTrialLen/timeBetweenSamples));
        e.eyes{eyeNum}.yPos = nan(edf.nTrials,ceil(maxTrialLen/timeBetweenSamples));
        e.eyes{eyeNum}.pupil = nan(edf.nTrials,ceil(maxTrialLen/timeBetweenSamples));
        % put in time in seconds
        e.eyes{eyeNum}.time = (0:(size(e.eyes{eyeNum}.xPos,2)-1))/edf.samplerate;
        
        % go through each trial and populate traces
        warning('off','MATLAB:interp1:NaNinY');
        for iTrial = 1:edf.nTrials
            disppercent(iTrial/edf.nTrials);
            % the times for this trial
            thisTrialTimes = startTime(iTrial):timeBetweenSamples:endTime(iTrial);
            % get data for xPos, yPos and pupil traces form edf data
            e.eyes{eyeNum}.xPos(iTrial,1:length(thisTrialTimes)) = interp1(edf.gazeBoth{eyeNum}.time,edf.gazeBoth{eyeNum}.x,thisTrialTimes,'linear',nan);
            e.eyes{eyeNum}.yPos(iTrial,1:length(thisTrialTimes)) = interp1(edf.gazeBoth{eyeNum}.time,edf.gazeBoth{eyeNum}.y,thisTrialTimes,'linear',nan);
            e.eyes{eyeNum}.pupil(iTrial,1:length(thisTrialTimes)) = interp1(edf.gazeBoth{eyeNum}.time,edf.gazeBoth{eyeNum}.pupil,thisTrialTimes,'linear',nan);
        end
        warning('on','MATLAB:interp1:NaNinY');
        disppercent(inf);
        e.eyes{eyeNum}.xPos = hDir * ((e.eyes{eyeNum}.xPos-(w/2))*xPix2Deg);
        e.eyes{eyeNum}.yPos = vDir * (((h/2)-e.eyes{eyeNum}.yPos)*yPix2Deg);
    end
end

%%% save original complete timecourses (to be used for microsaccade analysis, edited)

% startTime has start times for all trials
% end times are the start times for the next trial
% we then use those 

%e.complete.startTime = ceil((startTime - startTime(1))/timeBetweenSamples + 1);%in samples
%e.complete.endTime = ceil((endTime - startTime(1))/timeBetweenSamples + 1);

% e.complete.xPos = edf.gaze.x(1:e.complete.endTime(end));%only until the end of the last trial
% e.complete.yPos = edf.gaze.y(1:e.complete.endTime(end));

for iTrial = 1:edf.nTrials
    eyeStartIdx = find(edf.gaze.time == startTime(iTrial)); 
    eyeEndIdx = find(edf.gaze.time == endTime(iTrial)); 
    if ~isempty(eyeEndIdx) 
        EndSample(iTrial) = eyeEndIdx;
    else
        endTimeTemp = interp1(edf.gaze.time,edf.gaze.time,endTime(iTrial),'nearest');
        EndSample(iTrial) = find(edf.gaze.time == endTimeTemp);
    end
    if ~isempty(eyeStartIdx)
        StartSample(iTrial) = eyeStartIdx;
    else
        startTimeTemp = interp1(edf.gaze.time,edf.gaze.time,startTime(iTrial),'nearest');
        StartSample(iTrial) = find(edf.gaze.time == startTimeTemp);
    end
end
e.complete.startTime = StartSample - StartSample(1) + 1;
e.complete.endTime = EndSample - StartSample(1) + 1;
e.complete.pupil = edf.gaze.pupil(StartSample(1):EndSample(end));
e.complete.xPos = edf.gaze.x(StartSample(1):EndSample(end));
e.complete.yPos = edf.gaze.y(StartSample(1):EndSample(end));


% convert to device coordinates (deg)
w = stimfile.myscreen.screenWidth;
h = stimfile.myscreen.screenHeight;
xPix2Deg = stimfile.myscreen.imageWidth/w;
yPix2Deg = stimfile.myscreen.imageHeight/h;

hDir = 1;vDir = 1;
if isfield(stimfile.myscreen,'flipHV') && (length(stimfile.myscreen.flipHV) >=2 )
  if stimfile.myscreen.flipHV(1) hDir = -1;end
  if stimfile.myscreen.flipHV(2) vDir = -1;end
end

e.complete.xPos = hDir * ((e.complete.xPos-(w/2))*xPix2Deg);
e.complete.yPos = vDir * (((h/2)-e.complete.yPos)*yPix2Deg);
%%%%%%

%e.complete.pupil = edf.gaze.pupil(1:e.complete.endTime(end));
%e.complete.time = edf.gaze.time;
%e.complete.time = ceil((e.complete.time - startTime(1))/timeBetweenSamples + 1);%in samples

e.fixations = edf.fixations;
e.saccades = edf.saccades;
e.mgl = edf.mgl;


% display figure
if dispFig
    displayEyeTraces(e);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    displayEyeTraces    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayEyeTraces(e)

% get all combination of parameters
[stimvol names trialNums] = getStimvolFromVarname('_every_',e.stimfile.myscreen,e.stimfile.task,e.stimfile.taskNum,e.stimfile.phaseNum);

% remove empty trialNums
for iTrialType = 1:length(trialNums)
  emptyTypes(iTrialType) = isempty(trialNums{iTrialType});
end
stimvol = {stimvol{find(~emptyTypes)}};
names = {names{find(~emptyTypes)}};
trialNums = {trialNums{find(~emptyTypes)}};
figure;

% display trial by trial.
for iTrialType = 1:length(trialNums)
  % get color
  c = getSmoothColor(iTrialType,length(trialNums),'hsv');
  % display horizontal eye trace
  subplot(2,3,1:2);
  plot(e.eye.time,e.eye.xPos(trialNums{iTrialType},:)','Color',c);
  hold on
  % display vertical eye trace
  subplot(2,3,4:5);
  plot(e.eye.time,e.eye.yPos(trialNums{iTrialType},:)','Color',c);
  hold on
  % display as an x/y plot the median eye trace for this condition
  subplot(2,3,[3 6]);
  xPos = nanmedian(e.eye.xPos(trialNums{iTrialType},:));
  yPos = nanmedian(e.eye.yPos(trialNums{iTrialType},:));
  plot(xPos,yPos,'.','Color',c);
  hold on
end

% figure trimmings
hMin = -15;hMax = 15;
vMin = -15;vMax = 15;
subplot(2,3,1:2);
% yaxis(hMin,hMax);
ylim([hMin,hMax]);
% xaxis(0,max(e.eye.time));
xlim([0,max(e.eye.time)]);
xlabel('Time (sec)');
ylabel('H. eyepos (deg)');
subplot(2,3,4:5);
% yaxis(vMin,vMax);
ylim([vMin,vMax]);
% xaxis(0,max(e.eye.time));
xlim([0,max(e.eye.time)]);
xlabel('Time (sec)');
ylabel('V. eyepos (deg)');
subplot(2,3,[3 6]);
% xaxis(hMin,hMax);
xlim([hMin,hMax]);
% yaxis(vMin,vMax);
ylim([vMin,vMax]);
xlabel('H. eyepos (deg)');
ylabel('V. eyepos (deg)');
title('Median eye position by trial type');
axis square

%%%% isfile.m%%%%
%
%      usage: isfile(filename)
%         by: justin gardner
%       date: 08/20/03
%       e.g.: isfile('filename')
%    purpose: function to check whether file exists
%
function retval = isfile(filename)

if (nargin ~= 1)
  help isfile;
  return
end

% open file
fid = fopen(filename,'r');

% check to see if there was an error
if (fid ~= -1)
  fclose(fid);
  retval = 1;
else
  retval = 0;
end
%%%%%%%%%%
