%%  paramorg.m

%   usage: [T, weird] = paramorg(PupilData,N,R)
%   by: June Hee Kim
%   date: 11/01/2022 (last modified: 01/05/2023)
%   purpose: Indicates subjects with unwanted task conditions and produces
%   a table to check for task parameters (can modify which task parameters
%   to look for as needed).

%   INPUT: PupilData - extracted pupil data (cell)
%          R - number of runs for each subjects (double)
%          input variables could be extracted from extractPupilData.m 

%          minSFCount - minimum number of spatial frequency levels 
%          minContrastCount: 

%   OUTPUT: T - table consisting of subject number, number of contrast &
%           spatial frequency levels, whether recorded contrast levels matches with
%           that of task stimfile, nummber of runs, and sample rate. 

%           weird - list subjects with task variables not in interest or
%           pupil responses mostly consisting of NaNs.


%% 

function [T, weird] = paramorg(PupilData, R, varargin)

% initialize output variables 
weird = [];
contrastMismatch = [];

numSubj = size(R,2);

% check if PupilData and R is provided 
if nargin < 2 || isempty(PupilData) || isempty(R)
    error('PupilData and R are required inputs.');
end

% check if criteria values are provided as input
if nargin == 3
    minSFCount = varargin{1};
    choice = input('Define minContrastCount: ');
    minContrastCount = choice;
elseif nargin == 4
    minSFCount = varargin{1};
    minContrastCount = varargin{2};
else
    % use default values if not provided
    minSFCount = 5;
    minContrastCount = 2;
end

    for ii = 1:numSubj   
        if  ~isempty(PupilData{ii,1})
            clear e s x rewardType
            e = PupilData{ii,1};
            s = PupilData{ii,2};

            % initialize spatial frequency and contrast levels count
            numSpatialFreq = 0;
            numContrast = 0;

            for r = 1:sum(~isnan(R(:,ii)))
                samplerate{ii}(r,:) = s{r}.sampleRate;
                if ~isfield(e{r}, 'stimulus') || ~isfield(e{r}.stimulus,'freqVals') || ~isfield(e{r}.stimulus,'contrastVals')
                    spatFreq{ii}(r,:) = nan(1,5);
                    contrast{ii}(r,:) = nan(1,5);
                    clevels{ii}(r,:) = nan(1,5);
                    match{ii}(r,:) = nan;
                else
                    spatFreq{ii}(r,:) = e{r}.stimulus.freqVals;
                    contrast{ii}(r,:) = e{r}.stimulus.contrastVals;
                    a = length(unique(contrast{ii}(r,:)));
                    if a < 5
                        clevels{ii}(r,:) = [unique(contrast{ii}(r,:)) nan(1,5-a)];
                    else
                        clevels{ii}(r,:) = unique(contrast{ii}(r,:));                    
                    end
                    % check if values stored in task's stim files match with those in recorded data
                    if length(e{r}.block.contrast) ~= a
                        match{ii}(r,:) = 0; 
                        contrastMismatch = [contrastMismatch; ii];
                    else
                        match{ii}(r,:) = 1;
                    end
                end
                % store run and subject number 
                run{ii,:}(r,:) = r;
                sub{ii,:}(r,:) = ii;
            end
            
            numSpatialFreq = sum(~isnan(unique(spatFreq{ii}))); % subjects without stimulus info would be numSpatialFreq = 0
            numContrast = sum(~isnan(unique(contrast{ii})));
            
            % record subjects with less than two runs or doesn't meet the
            % specified criteria
            if sum(~isnan(R(:,ii))) <= 2 || numSpatialFreq ~= minSFCount || numContrast ~= minContrastCount 
                weird = [weird; ii];
            end
        else
            % record subjects with empty cells
            weird = [weird; ii]; 
        end         
    end

    % table
    varNames = {'Subject','Contrast','SpatialFrequency','ContrastMatch','Run','SampleRate'};    
    T = table(cell2mat(sub),cell2mat(clevels'),cell2mat(spatFreq'),cell2mat(match'),cell2mat(run),cell2mat(samplerate'),...
    'VariableNames',varNames);

    % check if the list includes subjects with contrast mismatch
    valuesNotInWeird = setdiff(unique(contrastMismatch), weird);
    weird = [weird; valuesNotInWeird];
    
     % list of subjects to exclude 
    weird = sort(weird);
end