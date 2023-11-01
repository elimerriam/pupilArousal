%% rwdRapid.m

%        $Id: rwdRapid.m,v 1.0 2019/01/18 zvi
%      usage: rwdRapid('rewardType=''H''','runNum=1','useStaircase=1','currBal=0','numTrials=15','displayName=''rm315''', 'threshStair1=0.15','threshStair2=0.4');
%         by: Zvi Roth & June Hee Kim
%       date: 01/18/2019 (last modified: 12/15/2022)
%    purpose: Measures arousal effect on stimulus responses.
%    Experiment is different for first run and for all subsequent runs.
%    First run uses 2 staircases, has feedback, and reward is 0.
%    rwdRapid('runNum=1','useStaircase=1','currBal=0','numTrials=15','displayName=rm315', 'threshStair1=0.15','threshStair2=0.4');

%    Subsequent runs have a fixed threshold, no feedback, and pay a reward. 
%    rwdRapid('rewardType=''H''','runNum=2','useStaircase=0','currBal=0','numTrials=17','displayName=rm315');


function [] = rwdRapid(varargin)
clear fixStimulus
pauseDuration = 8;
minThreshold = 0.15;
maxThreshold = 0.5;
nullTrials = [0 0 0 1];

% % evaluate the input arguments
getArgs(varargin, [], 'verbose=0');


% set default parameters

if ieNotDefined('displayName'), displayName = 'laptop'; end %%
if ieNotDefined('waitForBacktick')
    if strcmp(displayName, 'rm315') || strcmp(displayName, 'laptop')
        waitForBacktick = 0;
    else
        waitForBacktick = 1;
    end
end
if ieNotDefined('useStaircase'), useStaircase = 1; end
if ieNotDefined('threshStair1'), threshStair1 = 0.1; end
if ieNotDefined('threshStair2'), threshStair2 = maxThreshold; end
if useStaircase==0
    threshStair1 = (threshStair1+threshStair2)/2;
    threshStair1 = max(threshStair1, minThreshold);%threshold must be >0
    threshStair2 = threshStair1;
end


interTime = 2;
cueTime = 0.5; %cue display duration
responseTime=2;
stimLen = interTime;%should be a multiple of frameLen

if ieNotDefined('trialLen'),trialLen = 15;end %in seconds
if ieNotDefined('frameLen'),frameLen = stimLen/4; end
if ieNotDefined('innerEdge'),innerEdge = 1.3; end
if ieNotDefined('outerEdge'),outerEdge = 50; end
if ieNotDefined('rewardType'), rewardType = 'L'; end
if ieNotDefined('runNum'), runNum = 1; end
if ieNotDefined('probRwd'), probRwd = 0; end
if ieNotDefined('currBal'), currBal = 00; end
% if ieNotDefined('fixThresh'), fixThresh = 0.2; end
if ieNotDefined('numTRs'), numTRs = 170; end
TR=1.5;
if ieNotDefined('numTrials'), numTrials = ceil(TR*numTRs/trialLen); end


incrRwdL = -0.01;%reward decreases on every low reward run
incrRwdH = 0.04;%reward increases on every high reward run
initRwdL = 0.06;%reward for first low reward run
initRwdH = 0.5;%reward for first high reward run

if runNum==1
    rewardValue = 0;
elseif rewardType == 'H'
    rewardValue = initRwdH + incrRwdH * (runNum-1)/2;
elseif rewardType == 'L'
    rewardValue = initRwdL + incrRwdL * (runNum-1)/2;
    rewardValue = max(rewardValue,0.01);%don't want to reach zero
end

global stimulus;

% store parameters in stimulus variable
% update stimulus every 250 ms
stimulus.frameLen = frameLen;
% stimulus is on for 1.5 seconds
stimulus.stimLen = stimLen;

% trial is 15 seconds
stimulus.trialLen = trialLen;
% inner and outer edges
stimulus.inner = innerEdge;
stimulus.outer = outerEdge;


%reward parameters
stimulus.currBal = currBal;
stimulus.rewardType = rewardType;
stimulus.rewardValue = rewardValue;
stimulus.probRwd = probRwd;
stimulus.runNum = runNum;

% initalize the screen
myscreen.background = 'gray';
myscreen.autoCloseScreen = 0;
myscreen.allowpause = 1;
myscreen.saveData = 1;
myscreen.displayName = displayName;

myscreen = initScreen(myscreen);

global fixStimulus

fixStimulus.trialLength = trialLen;
fixStimulus.cueTime = cueTime;
fixStimulus.interTime = interTime;
fixStimulus.responseTime = responseTime;
fixStimulus.useStaircase = useStaircase;
fixStimulus.diskSize = 0;
fixStimulus.stairStepSize = 0.05;
fixStimulus.threshStair1 = threshStair1;
fixStimulus.threshStair2 = threshStair2;

fixStimulus.waitForBacktick = waitForBacktick;

fixStimulus.targetDigitCount = 2;%number of targets to respond 1 to
fixStimulus.targetDigit = 0; %target digit
fixStimulus.digitColor = [1 1 1]; %digit color
fixStimulus.digitSizeX = 0.25;
fixStimulus.digitSizeY = 0.2;
fixStimulus.digitWidth = 2;
fixStimulus.cueSize =0.55; %circle size in degrees
fixStimulus.cueWidth =0.07; %circle width
fixStimulus.cueColor = [0 1 1];
fixStimulus.responseColor = [1 1 0];
fixStimulus.cueContrast = 0.5;

if stimulus.runNum>1
   fixStimulus.incorrectColor = myscreen.background*[1 1 1];
   fixStimulus.correctColor = myscreen.background*[1 1 1];
else
   fixStimulus.incorrectColor = [0.8 0 0];
   fixStimulus.correctColor = [0 0.8 0];
end

[task{1} myscreen] = myFixStairInitTask(myscreen);
task{1}{1}.numTrials = numTrials;
% task{1}{1}.nTrials = numTrials;
task{1}{1}.response = zeros(numTrials,1);
task{1}{1}.correctResponse = zeros(numTrials,1);
task{1}{1}.correctness = zeros(numTrials,1) - 1;
%both tasks should wait for backtick to begin
task{1}{1}.waitForBacktick = waitForBacktick;
task{2}{1}.waitForBacktick = waitForBacktick;

%each segment is a different phase. Entire trial is a single orientation, contrast,and spatial frequency.
seglen = [fixStimulus.cueTime stimulus.frameLen * ones(1,stimulus.stimLen/stimulus.frameLen) stimulus.trialLen-(stimulus.stimLen+fixStimulus.cueTime)];

task{2}{1}.seglen = seglen;
task{2}{1}.synchToVol = zeros(size(seglen));


orientations = linspace(0, 180, 7);
stimulus.orientations = orientations(1:end-1);

% stimulus properties, block randomized
stimulus.contrasts = logspace(-0.7,0,2);
task{2}{1}.randVars.block.contrast = 1:length(stimulus.contrasts);
% task{2}{1}.randVars.block.contrast = logspace(-0.7,0,5);
stimulus.freqs = logspace(-0.3,0.5,5);
task{2}{1}.randVars.block.spatFreq = 1:length(stimulus.freqs);
% task{2}{1}.randVars.block.spatFreq = logspace(-0.3,0.5,5);
task{2}{1}.randVars.len_ = 20;

% task{2}{1}.randVars.block.orientation = 1:length(orientations);
task{2}{1}.randVars.block.nullTrial = nullTrials; %determines proportion of null trials

task{2}{1}.random = 1;
task{2}{1}.numTrials = numTrials;
% task{2}{1}.nTrials = numTrials;
task{2}{1}.collectEyeData = true;

% initialize the task
for phaseNum = 1:length(task{2})
    [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end
% myscreen.imageWidth = myscreen.imageWidth+5;


% do our initialization which creates the gratings
stimulus = myInitStimulus(stimulus,myscreen,task);
% stimulus.orientations = orientations;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%myscreen = eyeCalibDisp(myscreen);

%initial screen
mglClearScreen;
totalRwd = task{1}{1}.numTrials * stimulus.rewardValue; %incr;
mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
%     if stimulus.rewardType == 'H';
%         text = sprintf('This is a high-reward run');
text = sprintf('You will gain or lose up to');
mglTextDraw(text,[0 2]);
text = sprintf('$%0.2f',totalRwd);
mglTextSet('Helvetica',80,[1 1 1],0,0,0,0,0,0,0);
mglTextDraw(text,[0 0]);
text = sprintf('in this run based on performance');
mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
mglTextDraw(text,[0 -2]);
mglFlush;

mglWaitSecs(pauseDuration);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mglClearScreen(); mglFlush; mglClearScreen();
% mglFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,[0 1 1], [0 0]) % default fixation cross
mglFlush
% mglFixationCross(fixStimulus.fixWidth,fixStimulus.fixLineWidth,[0 1 1], [0 0]) % default fixation cross
if ~waitForBacktick
    mglStrokeText( '1', 0, 0, fixStimulus.digitSizeX, fixStimulus.digitSizeY, fixStimulus.digitWidth, [1 1 1], 0 ); 
    mglFlush
    mglWaitSecs(3);
end

phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    [task{2} myscreen phaseNum] = updateTask(task{2},myscreen,phaseNum);
    % update the fixation task
    [task{1} myscreen] = updateTask(task{1},myscreen,1);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end


% print out command for next run
%calculate reward for this run
if probRwd
    randP = randperm(task{1}{1}.numTrials);
    randRun = randP(1)
    rwd = task{1}{1}.correctness(randRun) * task{1}{1}.numTrials * stimulus.rewardValue;
    rwd = task{1}{1}.correctness(randRun) * task{1}{1}.numTrials * stimulus.rewardValue;
else
    rwd = sum(task{1}{1}.correctness) * stimulus.rewardValue;
end

stimulus.currBal = stimulus.currBal + rwd;
disp(sprintf('\n% --------------------------------------------- %\n'));

%calculate reward for next run
newRunNum = stimulus.runNum+1;
if stimulus.rewardType == 'H'
    newRewardType = 'L';
    %     newRewardValue = initRwdL + incrRwdL * newRunNum/2;
elseif stimulus.rewardType == 'L'
    newRewardType = 'H';
    %     newRewardValue = initRwdH + incrRwdH * newRunNum/2;
end

disp(['rwdRapid(''rewardType=''''' newRewardType ...
    ''''''',''runNum=' num2str(newRunNum) ...
    ''',''useStaircase=0' ...
    ''',''currBal=' num2str(stimulus.currBal) ...
    ''',''numTrials=' num2str(task{1}{1}.numTrials) ...
    ''',''displayName=''''' myscreen.displayName ...
    ''''''', ''threshStair1=' num2str(fixStimulus.stair{1}.threshold) ...
    ''', ''threshStair2=' num2str(fixStimulus.stair{2}.threshold) ''');']);


disp(sprintf('\n% --------------------------------------------- %\n'));



mglClearScreen;
mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);

if rwd<0
    text = sprintf('You LOST');
else
    text = sprintf('You gained');
end
mglTextDraw(text,[0 3]);
text = sprintf('$%0.2f', abs(rwd));
mglTextSet('Helvetica',70,[1 1 1],0,0,0,0,0,0,0);
mglTextDraw(text,[0 1.5]);
text = sprintf('in this run');
mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
mglTextDraw(text,[0 0]);
text = sprintf('Current Balance: $%0.2f',stimulus.currBal);
mglTextDraw(text,[0 -2]);
mglFlush;
mglWaitSecs(pauseDuration);


% if we got here, we are at the end of the experiment
myscreen.stimulus = stimulus; %save stimulus
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = startSegmentCallback(task, myscreen)
global stimulus;

% if ~task.thistrial.nullTrial


if any(task.thistrial.thisseg == stimulus.stimulusSegments) && (~task.thistrial.nullTrial)
    icontrast = task.thistrial.contrast;
    ifreq = task.thistrial.spatFreq;
    newOri = stimulus.oriNum;
    while stimulus.oriNum == newOri
        newOri = ceil(rand(1)*length(stimulus.orientations));
    end
    stimulus.oriNum = newOri;
    stimulus.tex = stimulus.allGratings{icontrast,ifreq};
    stimulus.rotation = stimulus.orientations(newOri);
end

if task.numTrials == task.trialnum && task.thistrial.thisseg==length(task.seglen)
    disp('last trial');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = screenUpdateCallback(task, myscreen)

global stimulus;

% clear the screen
mglClearScreen;

if any(task.thistrial.thisseg == stimulus.stimulusSegments) && (~task.thistrial.nullTrial)
    % draw the texture
    mglBltTexture(stimulus.tex, [0 0 stimulus.height stimulus.height], 0, 0, stimulus.rotation);
    mglBltTexture(stimulus.innerMaskTex, [0 0 stimulus.height stimulus.height], 0, 0, 0);
    mglBltTexture(stimulus.outerMaskTex, [0 0 stimulus.height stimulus.height], 0, 0, 0);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the dot stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task)
% keep an array that lists which of the segments we
% are presenting the stimulus in.
% stimulus.stimulusSegments = 1:stimulus.stimLen/stimulus.frameLen;
stimulus.stimulusSegments = 2:(length(task{2}{1}.seglen) - 1);

stimulus.pixRes = min(myscreen.screenHeight/myscreen.imageHeight, myscreen.screenWidth/myscreen.imageWidth);

% scale factor
stimulus.sFac = 1;
% spatial frequency
stimulus.sf = 1.4;

% which phases we will have
stimulus.numPhases = 16;
stimulus.phases = 0:(360-0)/stimulus.numPhases:360;

% size of stimulus
% stimulus.height = 0.5*floor(myscreen.imageHeight/0.5)+2;
% stimulus.width = 0.5*floor(myscreen.imageWidth/0.5)+2;
stimulus.height = 0.5*floor(myscreen.imageHeight/0.5)+18;
stimulus.height = 0.5*floor(myscreen.imageHeight/0.5)+22;
stimulus.width = stimulus.height;

stimulus.width = 0.5*floor(myscreen.imageWidth/0.5)+7;
stimulus.height = stimulus.width;

% size of annulus
% stimulus.outer = stimulus.height;
stimulus.outTransition = 0;

% stimulus.inner = 0.75;
stimulus.inTransition = 0;

% chose a sin or square
stimulus.square = 0;

% initial phase number
stimulus.phaseNum = 1;
% initial orientation number
stimulus.oriNum = 1;

% make a grating just to get the size
tmpGrating = mglMakeGrating(stimulus.width, stimulus.height, stimulus.sf, 0, stimulus.phases(1), stimulus.pixRes, stimulus.pixRes);
sz = size(tmpGrating,2);


% create mask for fixation and edge
out = stimulus.outer/stimulus.width;
in = stimulus.inner/stimulus.width;
twOut = stimulus.outTransition/stimulus.width;
twIn = stimulus.inTransition/stimulus.width;
outerMask = mkDisc(sz,(out*sz)/2,[(sz+1)/2 (sz+1)/2],twOut*sz,[1 0]);
innerMask = mkDisc(sz,(in*sz)/2,[(sz+1)/2 (sz+1)/2],twIn*sz,[0 1]);

% rescale mask to max out at 1
outerMask = outerMask/max(outerMask(:));
innerMask = innerMask/max(innerMask(:));

outerMask(:,:,4) = (-1*(outerMask*255))+255;
innerMask(:,:,4) = (-1*(innerMask*255))+255;

outerMask(:,:,1:3) = 128;
innerMask(:,:,1:3) = 128;

innerMask = uint8(permute(innerMask, [3 1 2]));
outerMask = uint8(permute(outerMask, [3 1 2]));

stimulus.innerMaskTex = mglCreateTexture(innerMask);
stimulus.outerMaskTex = mglCreateTexture(outerMask);


% make a grating again, but now scale it
tmpGrating = mglMakeGrating(stimulus.width/stimulus.sFac, stimulus.height/stimulus.sFac, stimulus.sf, 0, stimulus.phases(1), stimulus.pixRes, stimulus.pixRes);
r = uint8(permute(repmat(tmpGrating, [1 1 4]), [3 1 2]));
stimulus.tex = mglCreateTexture(r,[],1);


for icontrast = 1:length(stimulus.contrasts)
    for ifreq=1:length(stimulus.freqs)
        newPhase = ceil(rand(1)*stimulus.numPhases);
        newOri = 1;
        grating = mglMakeGrating(stimulus.width/stimulus.sFac, stimulus.height/stimulus.sFac, stimulus.freqs(ifreq)*stimulus.sFac, ...
            stimulus.orientations(newOri), stimulus.phases(newPhase), stimulus.pixRes, stimulus.pixRes);
        grating = grating*stimulus.contrasts(icontrast);
        % scale to range of display
        grating = 255*(grating+1)/2;
        % make it rgba
        grating = uint8(permute(repmat(grating, [1 1 4]), [3 1 2]));
        grating(4,:,:) = 256;
        stimulus.allGratings{icontrast,ifreq} = mglCreateTexture(grating);
    end
end

% fixStairInitTask.m
%
%        $Id$
%      usage: [fixTask myscreen] = fixStairInitTask(myscreen)
%         by: justin gardner
%       date: 09/07/06
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: Implements a fixation task. In this task, a stream of digits appears
%             at the screen center. a blue circle flashes around fixation, signalling 
%             to the subject to begin counting zeros, and 
%             after a few seconds a yellow circle appears. 
%             If 2 zeros were counted subject presses 2. Otherwise 
%             subject presses 1.
%             With feedback, the circle will then turn green or red to indicate
%             correct or incorrect responses. The speed of the digits is
%             controlled by a 2 down 1 up staircase to keep task difficulty
%             the same.
%
%             See testExperiment.m for how this is used in a task. If you want
%             to change parameters, before you call fixStairInitTask, set
%             appropriate fields of the global variable fixStimulus. e.g.:
%
%             global fixStimulus
%             fixStimulus.interTime = 1;
%
%             See the code, for a list of all parameters that can be changed.
%
function [task, myscreen] = myFixStairInitTask(myscreen)

% check arguments
if ~any(nargin == [1])
    help fixDispStairInitTask
    return
end

% create the stimulus for the experiment, use defaults if they are
% not already set
global fixStimulus;
myscreen = initStimulus('fixStimulus',myscreen);

minThreshold=0.05;

% if ~isfield(fixStimulus,'precueTime'); fixStimulus.precueTime = 3; end%trial length in seconds
if ~isfield(fixStimulus,'trialLength'); fixStimulus.trialLength = 18; end%trial length in seconds
if ~isfield(fixStimulus,'cueTime'); fixStimulus.cueTime = 0.1; end%cue length in seconds
if ~isfield(fixStimulus,'interTime'); fixStimulus.interTime = 3; end%time between cues in seconds
if ~isfield(fixStimulus,'responseTime'); fixStimulus.responseTime = 2; end%time for response after 2nd cue
% if ~isfield(fixStimulus,'responseTime'); fixStimulus.responseTime = 2; end%time for response after 2nd cue
% fixStimulus.postResponseTime = fixStimulus.trialLength - (fixStimulus.precueTime + 2*fixStimulus.cueTime + fixStimulus.interTime + fixStimulus.responseTime);

fixStimulus.postResponseTime = fixStimulus.trialLength - (2*fixStimulus.cueTime + fixStimulus.interTime + fixStimulus.responseTime);

if ~isfield(fixStimulus,'targetDigitCount'); fixStimulus.targetDigitCount = 2; end%number of targets to respond 1 to
if ~isfield(fixStimulus,'targetDigit'); fixStimulus.targetDigit = 0; end%target digit

if ~isfield(fixStimulus,'digitColor'); fixStimulus.digitColor = [1 1 1]; end%digit color
if ~isfield(fixStimulus,'digitSize'); fixStimulus.digitSize = 30; end%digit size
if ~isfield(fixStimulus,'digitFont'); fixStimulus.digitFont = 'Helvetica'; end%digit font
if ~isfield(fixStimulus,'cueSize'); fixStimulus.cueSize =0.6; end%circle size in degrees
if ~isfield(fixStimulus,'cueWidth'); fixStimulus.cueWidth =0.07; end%circle width
if ~isfield(fixStimulus,'digitSizeX'); fixStimulus.digitSizeX =0.5; end%
if ~isfield(fixStimulus,'digitSizeY'); fixStimulus.digitSizeY =0.4; end%
if ~isfield(fixStimulus,'digitWidth'); fixStimulus.digitWidth =2; end%

if ~isfield(fixStimulus,'pedestal'); fixStimulus.pedestal = 0.4; end
if ~isfield(fixStimulus,'stairUp'); fixStimulus.stairUp = 1; end
if ~isfield(fixStimulus,'stairDown'); fixStimulus.stairDown = 2; end
if ~isfield(fixStimulus,'stairStepSize'); fixStimulus.stairStepSize = 0.05; end
if ~isfield(fixStimulus,'stairUseLevitt'); fixStimulus.stairUseLevitt = 0; end
% if ~isfield(fixStimulus,'stimColor'); fixStimulus.stimColor = [0 1 1]; end
if ~isfield(fixStimulus,'responseColor'); fixStimulus.responseColor = [1 1 0]; end
% if ~isfield(fixStimulus,'interColor'); fixStimulus.interColor = [0 1 1]; end
if ~isfield(fixStimulus,'correctColor'); fixStimulus.correctColor = [0 0.8 0]; end
if ~isfield(fixStimulus,'incorrectColor'); fixStimulus.incorrectColor = [0.8 0 0]; end
if ~isfield(fixStimulus,'cueColor'); fixStimulus.cueColor = [0 1 1]; end
if ~isfield(fixStimulus,'cueContrast'); fixStimulus.cueContrast = 0.5; end

if ~isfield(fixStimulus,'diskSize'); fixStimulus.diskSize = 1; end
if ~isfield(fixStimulus,'pos'); fixStimulus.pos = [0 0]; end
if ~isfield(fixStimulus,'verbose'); fixStimulus.verbose = 1;end



% set text properties
mglTextSet(fixStimulus.digitFont,fixStimulus.digitSize,fixStimulus.digitColor,0,0,0,0,0,0,0);

% mglTextSet('Helvetica',64,[1 1 1],0,0,0,0,0,0,0);

% create a fixation task
task{1}.seglen = [fixStimulus.cueTime fixStimulus.interTime ...
    fixStimulus.cueTime fixStimulus.responseTime fixStimulus.postResponseTime];

% if fixStimulus.waitForBacktick
%     task{1}.seglen(end) = task{1}.seglen(end) - 0.5;
% end
task{1}.getResponse = [0 0 1 1 0];
task{1}.synchToVol = zeros(size(task{1}.seglen));
% if fixStimulus.waitForBacktick
%     task{1}.synchToVol(end) = 1;
% end


[task{1}, myscreen] = addTraces(task{1}, myscreen, 'segment', 'phase', 'response');

% init a 2 down 1 up staircase
fixStimulus.stair{1} = upDownStaircase(1,2,fixStimulus.threshStair1,fixStimulus.stairStepSize,0);
fixStimulus.stair{1}.minThreshold = minThreshold;
fixStimulus.stair{2} = upDownStaircase(1,2,fixStimulus.threshStair2,fixStimulus.stairStepSize,0);
fixStimulus.stair{2}.minThreshold = minThreshold;

% init the task
[task{1}, myscreen] = initTask(task{1},myscreen,@fixStartSegmentCallback,@fixDrawStimulusCallback,@fixTrialResponseCallback,@fixTrialStartCallback);

[task{1}, myscreen] = addTraces(task{1}, myscreen, 'fixStair');

% keep the correct and incorrect counts
task{1}.correct = 0;
task{1}.incorrect = 0;

% set cue contrast
fixStimulus.cueColor = (fixStimulus.cueColor - myscreen.background*[0.5 0.5 0.5]) * fixStimulus.cueContrast + myscreen.background*[0.5 0.5 0.5];
fixStimulus.responseColor= (fixStimulus.responseColor - myscreen.background*[0.5 0.5 0.5]) * fixStimulus.cueContrast + myscreen.background*[0.5 0.5 0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = fixTrialStartCallback(task, myscreen)

global fixStimulus;

% which staircase will be used this trial
fixStimulus.whichStair = round(rand)+1;

% determine whether the number of target digit will be equal to the target number
task.thistrial.correctResponse = round(rand)+1;

%determine digit stream speed, according to staircase
task.thistrial.SOA = fixStimulus.stair{fixStimulus.whichStair}.threshold;

task.thistrial.ISI = task.thistrial.SOA/2;
task.thistrial.numDigits = ceil(fixStimulus.trialLength/task.thistrial.SOA);
task.thistrial.numDigitTypes = floor(fixStimulus.interTime/(task.thistrial.SOA*fixStimulus.targetDigitCount));
task.thistrial.numDigitTypes = min(task.thistrial.numDigitTypes,9);

numCueDigits = ceil((fixStimulus.cueTime*2+fixStimulus.interTime)/task.thistrial.SOA);%floor?
badStream=1;
if task.thistrial.correctResponse==1 %not target number
    while badStream
        tempDigitStream = floor((task.thistrial.numDigitTypes+1)*rand(task.thistrial.numDigits,1));
        tempCueTargets = sum(tempDigitStream(1:numCueDigits)==fixStimulus.targetDigit);
        if tempCueTargets ~= fixStimulus.targetDigitCount
            badStream=0;
        end
    end
else % target number
    while badStream
        tempDigitStream = floor((task.thistrial.numDigitTypes+1)*rand(task.thistrial.numDigits,1));
        tempCueTargets = sum(tempDigitStream(1:numCueDigits)==fixStimulus.targetDigit);
        if tempCueTargets == fixStimulus.targetDigitCount
            badStream=0;
        end
    end
end
task.thistrial.digitStream = tempDigitStream;
task.thistrial.numCueTargets = tempCueTargets;
task.numCueTargets(task.trialnum) = tempCueTargets;
% numCueTargets = sum(task.thistrial.digitStream(1:numCueDigits)==fixStimulus.targetDigit);
% task.thistrial.numCueTargets = numCueTargets;
if fixStimulus.verbose
    disp(['trial = ' num2str(task.trialnum) ' stair = ' num2str(fixStimulus.whichStair) ...
        ' threshold = ' num2str(fixStimulus.stair{fixStimulus.whichStair}.threshold,'%.2f')...
        ' #targets = ' num2str(task.thistrial.numCueTargets) ]);
    %   disp(sprintf('sigint = %i threshold = %0.2f',task.thistrial.sigInterval,fixStimulus.threshold));
end

% task.thistrial.targetCounter = 0;
task.thistrial.digitIndex = 0;
% task.thistrial.t0 = mglGetSecs;%start stopwatch for this trial


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = fixStartSegmentCallback(task, myscreen)

global fixStimulus;

% if this is the cue segment
if task.thistrial.thisseg == 1
    fixStimulus.thisColor = fixStimulus.cueColor;
elseif task.thistrial.thisseg == 3
    fixStimulus.thisColor = fixStimulus.responseColor;
else %beginning of another segment, no cue
    fixStimulus.thisColor = myscreen.background*[1 1 1];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called every frame udpate to draw the fixation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = fixDrawStimulusCallback(task, myscreen)

global fixStimulus;

timeIntoTrial = mglGetSecs(task.thistrial.trialstart);
trialDigitIndex = ceil(timeIntoTrial/task.thistrial.SOA);

%check whether we're moving on to the next digit
if trialDigitIndex > task.thistrial.digitIndex
    task.thistrial.digitIndex = task.thistrial.digitIndex+1;%should now equal trialDigitIndex
%     mglDeleteTexture(fixStimulus.displayText);
    task.thistrial.currentDigit = task.thistrial.digitStream(task.thistrial.digitIndex);
    task.thistrial.currentText = num2str(task.thistrial.currentDigit);

    %if this is a new digit we start a digit timer to know when to stop presenting it
    task.thistrial.digitTimer = mglGetSecs;
end
% the digit is presented until the ISI begins
if mglGetSecs(task.thistrial.digitTimer)<(task.thistrial.SOA - task.thistrial.ISI)
    mglStrokeText( task.thistrial.currentText, 0, 0, fixStimulus.digitSizeX, fixStimulus.digitSizeY, fixStimulus.digitWidth, [1 1 1], 0 ); 
end

if task.thistrial.thisseg == 1 || task.thistrial.thisseg == 3 || task.thistrial.thisseg == 4
    nslices = 32;
    mglGluAnnulus(0,0,fixStimulus.cueSize-fixStimulus.cueWidth/2, fixStimulus.cueSize+fixStimulus.cueWidth/2, fixStimulus.thisColor',nslices);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called when subject responds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = fixTrialResponseCallback(task, myscreen)

global fixStimulus;

% % % check whether targetCounter is equal to numCueTargets
% disp(task.thistrial.targetCounter == task.thistrial.numCueTargets);
% disp(num2str(task.thistrial.targetCounter));
% disp(num2str(task.thistrial.numCueTargets));


%determine what the correct response is
%1 if not the target number, 2 if it is the target number
correctResponse= 1 + (task.numCueTargets(task.trialnum)==fixStimulus.targetDigitCount);
task.correctResponse(task.trialnum)  = correctResponse;

%save response 
task.response(task.trialnum) = find(task.thistrial.buttonState);

% get correct or incorrect
response = task.response(task.trialnum) == correctResponse;
response = response(1);

% task.correctness(task.trialnum) = 2*response(1)-1;%1 or -1
task.correctness(task.trialnum) = 2*response-1;%1 or -1

if response
    % set to correct fixation color
    fixStimulus.thisColor = fixStimulus.correctColor;
    % set trace to 2 to indicate correct response
    myscreen = writeTrace(2,task.fixStairTrace,myscreen);
    % and update correct count
    task.correct = task.correct+1;
    if fixStimulus.verbose
        disp('YES');
    end
else
    % set to incorrect fixation color
    fixStimulus.thisColor = fixStimulus.incorrectColor;
    % set trace to -2 to indicate incorrect response
    myscreen = writeTrace(-2,task.fixStairTrace,myscreen);
    % and update incorrect count
    task.incorrect = task.incorrect+1;
    if fixStimulus.verbose
        disp('NO');
    end
end

% update staircase
if fixStimulus.useStaircase
    fixStimulus.stair{fixStimulus.whichStair} = upDownStaircase(fixStimulus.stair{fixStimulus.whichStair}, response);
end