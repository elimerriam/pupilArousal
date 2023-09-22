%% Bootstrapping 

%   usage: diffmean = bootstrapPA(PupilData,R,weird,numboots)
%   by: June Hee Kim
%   date: 12/01/2022 (last modified: 01/05/2023)
%   purpose: Validates the regression analysis used to calculate response amplitudes 
%            that it did not skew the results, diminishing any true effect of reward
%            on the stimulus-evoked response. 

%   Specifically does the following: Randomly splits null trial pupil responses into
%   pseudo null and stim. Computes pseudo null response template (mean response) and
%   subtracts from the pseudo stim trials. Computes peudo stim response template and response
%   amplitudes for appropriate trial types (same method in PupilArousal.m). Finds averaged difference between high
%   and low reward conditions within each trial types. 


%   INPUT: PupilData - extracted pupil data (cell)
%          R - number of viable runs for each subjects (double)
%          weird - list of subjects to exclude
%          numboots = number of bootstrap iteration. 10000 iterationns if
%          not specified
%          input variables could be extracted from extractPupilData.m 

%   OUTPUT: diffmean - averaged difference between high and low reward 
%           conditions within each trial types. 
%           diffmean(# of bootsrap iteration , trial type)

%%
function diffmean = bootstrapPA(PupilData,R,weird,numboots)

if ieNotDefined('numboots')
    numboots = 10000; % how many bootstrap iterations to perform per subject
end

numSubj = size(R,2);

for ii = 1:numSubj
    if  ~ismember(ii,weird)        
        clear e s x rewardType contrastLevels spatrFreqLavels
        e = PupilData{ii,1};
        s = PupilData{ii,2};
        x = PupilData{ii,3};
        rewardType = PupilData{ii,4};      
        
        %%% z-score every run %%%
        clear pupil trialLength pupilCat Zpupil
        for r = 1:sum(~isnan(R(:,ii)))
            pupil{r,:} = e{r}.eye.pupil(1:17,:);
        end
            trialLength = min(cellfun('size',pupil, 2));
        for r = 1:sum(~isnan(R(:,ii)))
            pupil{r,:} = pupil{r}(1:17,1:trialLength);
        end
        pupilCat = cell2mat(pupil);
        for r = 1:sum(~isnan(R(:,ii)))
            Zpupil{r,:} = (pupil{r,:} - mean(mean(pupilCat,'omitnan')))/std(pupilCat(:),0,'omitnan');
        end
        
        clear ZpupilResp stim null stimtonic nulltonic pupilResp
        for r = 1:sum(~isnan(R(:,ii)))
            pupilResp = Zpupil{r};
            trials(r,:) = x{r}.randVars.nullTrial(:,1:17);

            %%% Lowpass Butterworth Filter (3rd order) %%%
            ZpupilResp = preprocessPupil(pupilResp);
            
            for n = 0:1 % 0=stim, 1=null, extracting stim and null trials for pupil resp
                idx = trials(r,:) == n;
                h = ZpupilResp(idx,:);
                if n == 0
                    stim{r} = h;
                else
                    null{r} = h;
                end
            end    
        end 

        minTrial = min(cellfun('size',null, 2));
        clear nullmat rwdmat 
        for r = 1:sum(~isnan(R(:,ii)))
            % trim trial length 
            null{r} = null{r}(:,1:minTrial);
            col = size(null{r},1);
            % store reward variables 
            if rewardType{r} == 'H'
                rwd{r} = repelem('H',col,1);
            else
                rwd{r} = repelem('L',col,1);
            end
        end
        rwdmat = cell2mat(rwd');
        nullmat = cell2mat(null');

        for boot = 1:numboots
            fprintf('subject: %d bootnum: %d\n',ii, boot)
            clear k null1 null0 rwd1 rwd0 IRF_null1 IRF_null0 zIRF_null1 zIRF_null0 
            % randomly split null trials and respective reward conditions into half 
            k = randperm(size(nullmat,1));
            null1 = nullmat(k(1:round(length(k)/2)),:); % pesudo null trials
            rwd1 = rwdmat(k(1:round(length(k)/2)),:);
            null00 = nullmat(k((round(length(k)/2)+1):end),:); % pseudo stim trials
            rwd0 = rwdmat(k((round(length(k)/2)+1):end),:);
            
            % treat null1 as null trials and null00 as stim trials, compute IRF (z-scored)
            reward = ['H';'L'];
            for m = 1:2 
                clear IRF_null1 IRF_null0 idxnull0
                IRF_null1 = nanmean(null1(rwd1 == reward(m),:),1);
                zIRF_null1(m,:) = (IRF_null1 - nanmean(IRF_null1))/nanstd(IRF_null1);

                idxnull0 = rwd0 == reward(m);
                IRF_null0 = nanmean((null00(idxnull0,:) - IRF_null1),1);
                null0(idxnull0,:) = null00(idxnull0,:) - IRF_null1; % subtract pesudo null trial IRF
                zIRF_null0(m,:) = (IRF_null0 - nanmean(IRF_null0))/nanstd(IRF_null0);
            end

            % calculate response amplitude 
            clear HrwdBeta LrwdBeta
            HrwdBeta = zeros(size(null1,1),2);
            LrwdBeta = zeros(size(null1,1),2);
            for n = 0:1
                if n == 0
                    vars = null0; 
                    z_IRF = zIRF_null0;
                    rwdLevels = rwd0;
                else
                    vars = null1; 
                    z_IRF = zIRF_null1;
                    rwdLevels = rwd1;
                end
                clear beta
                for t = 1:size(vars,1)
                    clear Y ZIRF
                    if rwdLevels(t) == 'H'
                        ZIRF = z_IRF(1,:);
                    else
                        ZIRF = z_IRF(2,:);
                    end
                    Y = vars(t,:);
                    betaVal = calculateBeta(Y, ZIRF);
                    beta(t,:) = betaVal;
                end

                % differentiate response amplitude by reward condition 
                for t = 1:size(vars,1)
                    if rwdLevels(t) == 'H'
                        HrwdBeta(t,n+1) = beta(t,:); %col1 = stim, col2 =null;
                    else
                        LrwdBeta(t,n+1) = beta(t,:);
                    end
                end
                clear vars z_IRF rwdLevels 
            end 
            % difference between high and low reward (per each iteration for each subj)
             
            for n= 0:1
                % column 1~3 = stim betas column 4~6 = null betas
                diff{ii}(boot,(n*3)+1) = nanmean(nonzeros(LrwdBeta(:,n+1))); % low reward betas averaged over trials
                diff{ii}(boot,(n*3)+2) = nanmean(nonzeros(HrwdBeta(:,n+1))); % high reward betas averaged over trials
                diff{ii}(boot,(n*3)+3) = nanmean(nonzeros(HrwdBeta(:,n+1))) - nanmean(nonzeros(LrwdBeta(:,n+1))); % difference                  
            end
        end

    end        
end

% get rid of empty cells
diff = diff(~cellfun('isempty',diff));

% average differences in response amplitudes across repetitions for each subject 
clear result diffmean
for n = 0:1
    for r = 1:numboots
        for i = 1:length(diff)
            result{r}(i,:) = diff{i}(r,(n*3)+3);
        end
        diffmean(r,n+1) = nanmean(result{r});
    end
    clear result
end

end

%% Functions 

% function to preprocess pupil data (Z-score and run through Butterworth 3rd order filter)
function ZpupilResp = preprocessPupil(pupilResp)

    fc = 4; % cutoff frequency
    fs = 500; % data sampling frequency
    [b, a] = butter(3, fc / (fs / 2));

    ZpupilResp = NaN(size(pupilResp));
    
    for tl = 1:size(pupilResp, 1)
        dataOriginal = pupilResp(tl, :);
        dataIn = [NaN(1, 200), dataOriginal, NaN(1, 200)];
        dataIn(1, 1:250) = fillmissing(dataIn(1, 1:250), 'linear', 2, 'EndValues', 'nearest');
        dataIn(1, end - 250:end) = fillmissing(dataIn(1, end - 250:end), 'linear', 2, 'EndValues', 'nearest');
        
        if sum(~isnan(pupilResp(tl, :))) > 1 && sum(isnan(pupilResp(tl, :))) ~= size(pupilResp, 2)
            tf = isnan(dataIn);
            ix = 1:numel(dataIn);
            dataIn(tf) = interp1(ix(~tf), dataIn(~tf), ix(tf)); % Interpolate missing data
            dataOut = filter(b, a, dataIn); % Filter forward
            dataInRev = fliplr(dataOut);
            dataInRev = fillmissing(dataInRev, 'linear', 2, 'EndValues', 'nearest');
            dataOutRev = fliplr(filter(b, a, dataInRev)); % Filter in reverse
            ZpupilResp(tl, :) = dataOutRev(201:end - 200);
        elseif sum(~isnan(pupilResp(tl, :))) == 1
            ZpupilResp(tl, :) = dataOriginal;
        elseif sum(isnan(pupilResp(tl, :))) == size(pupilResp, 2)
            ZpupilResp(tl, :) = dataOriginal;
        end
    end
end

% function to calculate beta weights
function betaVal = calculateBeta(Y, ZIRF)
    X = [ZIRF; ones(1, size(ZIRF, 2))];
    nanIdx = isnan(Y) | isnan(ZIRF);
    B = (X(:, ~nanIdx))' \ (Y(~nanIdx))';   
    if B(1) ~= 0
        betaVal = B(1);
    else
        betaVal = NaN;
    end
end


