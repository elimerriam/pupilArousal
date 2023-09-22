% extractPupilData
%
%   usage: extractPupilData()
%   by: June Hee Kim
%   date: 12/01/2022 (last modified: 06/08/2022)
%   purpose: extracts and organizes eye tracking data from eye tracking experiment
%   INPUT:
%   D : working directory (make sure it contains stimfile .mat and eyelink data .edf)
%       Note: files for each subject must be named as the following format 
%       first three letters of first and last name, date collected ex. AleElc_2019-02-04
%   OUTPUT:
%   PupilData{N, 1}     e(1, ~isnan(R)) -- extracted eye(position, diameter, saccade, etc..) by viable run
%   PupilData{N, 2}     s(1, ~isnan(R)) -- stimulus, screen, and task info by viable run       
%   PupilData{N, 3}     x(1, ~isnan(R)) -- task variable information by viable run       
%   PupilData{N, 4}     rewardType(1, ~isnan(R)) -- reward info by viable run
%   (could be commented out depending on the experiment)
%   R - viable runs for each subject
%   N - subject(file) names

%[PupilData,R,N] = extractPD();

function [PupilData,R,N] = extractPupilData()
% retrieve directory information
D = pwd; %'/Users/kimjuh/Library/CloudStorage/Box-Box/backup/PCDM-main/June';
cd(D);
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'});

%  ollect all the runs from stim file
for ii = 1:numel(N) % subjects
    T = dir(fullfile(D,N{ii},'*.mat'));
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for r = 1:numel(C)
        R(r,ii)= str2double(C{r}(12:13));
        temp{ii}{r} = fullfile(D,N{ii},C{r});
    end
    R(R == 0) = NaN;
end

% extract e s x, check for readable runs, update R
for ii = 1:numel(N)
    clear e s x rewardType
    for r = 1:sum(~isnan(R(:,ii)))
        stimfile = temp{ii}{r};
        e{r} = myGetTaskEyeTraces(stimfile, 'removeBlink=3');
        % excluding empty or runs w/ less than 17 trials or 2/3 of runs all
        % consisting NaNs
        if isempty(e{r}) || ~isfield(e{r},'eye') || size(e{r}.eye.pupil,1)<17 || sum(sum(isnan(e{r}.eye.pupil))) >= (numel(e{r}.eye.pupil)*(2/3))
            temp{ii}{r} = NaN; %for record keeping 
            R(r,ii) = NaN;
            s{r} = [];
            x{r} = [];
            rewardType{r} = [];
            e{r} = [];
        else
            s{r} = load(stimfile);
            x{r} = getTaskParameters(stimfile);
            samplerate = s{r}.myscreen.eyetracker.params.sampleRate;
            if isfield(s{r}.myscreen,'stimulus')
                rewardType{r} = s{r}.myscreen.stimulus.rewardType;
            else
                if mod(r,2)==0
                    rewardType{r} = 'H';
                else
                    rewardType{r} = 'L';
                end
            end         
        end
    end
    PupilData{ii,1} = e(~cellfun('isempty',e));
    PupilData{ii,2} = s(~cellfun('isempty',s));
    PupilData{ii,3} = x(~cellfun('isempty',x));
    PupilData{ii,4} = rewardType(~cellfun('isempty',rewardType));
end
% uncomment below if you need to save the dataset
%save('PDVCSSMRI500blink3_2.mat','PupilData','N','R');
end


