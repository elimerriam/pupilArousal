%%   PupilArousal.m

%   usage: PupilArousal(PupilData, R)
%   by: June Hee Kim & Zvi Roth
%   date: 06/03/2022 (last modified: 01/05/2023)
%   purpose: Filters eye data and calculates response amplitudes per trials 
%   separately for task-related and stimulus-evoked responses to compare
%   between reward conditions and stimulus properties.
%   Kim J, Yin C, Merriam EP, Roth ZN (2023) Pupil Size Is Sensitive to Low-Level Stimulus Features, Independent of Arousal-Related Modulation. ENeuro 


%   Specifically does the following: z-score pupil sizes within each run,
%   low pass Butterworth (3rd order) filtering of pupil data,
%   calculating stim & null response template, and estimating trial wise
%   response amplitude through linear regression. Our regression analysis for 
%   estimating response amplitude has been evaluated using bootstrapping
%   procedure.

%   INPUT: PupilData - extracted pupil data (cell)
%          R - number of viable runs for each subjects (double)
%          input variables could be extracted from extractPupilData.m 
%          boots - set boots = 0 to turn off bootstrapping option.
%          Otherwise the function runs bootstrapping for 10000 iterations. 

%   OUTPUT: out - structure containing 
%           1) parameters - table of task parameters for each subject 
%           (refer to paramorg.m for details)
%           2) betasAllSubjects - cells with calculated response amplitudes
%           for each trials (row,trials)
%           3) betasAverage - averaged betas under each condition per
%           subjects (row = subject, column 1 = stim low reward, column 2 =
%           stim high reward, column 3 = null low reward, column 4 = null
%           high reward)
%           4) stim - statistical results for stim trial 
%           (stimulus evoked pupil response) analysis; t-test and ANOVA each containing
%           p value, C.I (tbl for ANOVA), and stats
%           5) null - statistical results for null trial (task evoked pupil
%           response) analysis
%           6) actualDifference - averaged low reward beta subtracted from high
%           reward beta for each trial types. (column 1 = stim, column 2 =
%           null)
%           7) Respective graphs and ANOVA tables



function [out] = PupilArousal(PupilData,R,boots)
%% Response amplitude 
% calculating trialwise response amplitude (beta weight) and comparing 
% between reward and trial type conditions 

if nargin < 2 || isempty(PupilData) || isempty(R)
    error('PupilData and R are required inputs.');
end

numSubj = size(R,2);

[T, weird] = paramorg(PupilData,R);
out.parameters = T;

for ii = 1:numSubj
    %%% subjects to exclude based on parameterorganizer.m %%%
    if  ~ismember(ii,weird)
        %%% call in data for each subject from PupilData %%%
        clear e s x rewardType
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
        
        clear pupilResp idx h stim null minTrial 

        for r = 1:sum(~isnan(R(:,ii)))
            trials{ii}(r,:) = x{r}.randVars.nullTrial(:,1:17);
            spatFreq{ii}(r,:) = x{r}.randVars.spatFreq(:,1:17);
            contrast{ii}(r,:) = x{r}.randVars.contrast(:,1:17);
            react{ii}(r,:) = e{r}.reactionTime(:,1:17);
            if sum(mod(contrast{ii}(r,:),1)) ~= 0 || sum(mod(spatFreq{ii}(r,:),1)) ~= 0%
                [contrastVal,rankidx] = sort(unique(x{r}.randVars.contrast),'descend');
                [SFVal,rankidxsf] = sort(unique(x{r}.randVars.spatFreq),'descend');
                for replace = 1:length(contrastVal)
                    contrast{ii}(r,contrast{ii}(r,:) == contrastVal(replace)) = rankidx(replace);
                    spatFreq{ii}(r,spatFreq{ii}(r,:) == SFVal(replace)) = rankidxsf(replace);
                end
            end                        
            %correct{ii}(r,:) = s{r}.correctness(1:17,:)';
            pupilResp = Zpupil{r};

            %%% Lowpass Butterworth Filter (3rd order) %%%
            ZpupilResp = preprocessPupil(pupilResp);
            %%%%
            
            for n = 0:1 % n = 0 indicates stim trials, n = 1 indicates null trials 
                idx = trials{ii}(r,:) == n;
                h = ZpupilResp(idx,:);
                if n == 0
                    stim{r} = h;
                else
                    null{r} = h;
                end
            end
        end
        
        %%% trim trials to minimum trial length %%%
        minTrial = min(cellfun('size',null, 2));
        for r = 1:sum(~isnan(R(:,ii)))
            null{r} = null{r}(:,1:minTrial);
            stim{r} = stim{r}(:,1:minTrial);
        end
        
        clear nullmat stimmat IRF IRF_stim zIRF zIRF_stim
        reward = ['L';'H'];
        
        %%% calculate stim and null response templates for each reward level %%%
        for m = 1:2 % reward level; m = 1 indicates low reward and m = 2 high reward
            for n = 1:2 
                % null trial IRF (z-scored)
                if n == 1
                    clear nullRwd nullmat IRF
                    nullRwd = null(contains(rewardType,reward(m)));
                    nullmat = cell2mat(nullRwd');
                    IRF = nanmean(nullmat);
                    zIRF(m,:) = (IRF - nanmean(IRF))/nanstd(IRF);
                else
                    % stim trial IRF (z-scored)
                    for r = 1:sum(~isnan(R(:,ii)))
                        if rewardType{r} == reward(m)
                            stim{r} = stim{r} - IRF; % subtract null trial response template from stim trials 
                        else
                            stim{r} = stim{r}; 
                        end
                    end
                end
            end
            clear stimRwd stimmat IRF_stim
            stimRwd = stim(contains(rewardType,reward(m)));
            stimmat = cell2mat(stimRwd');
            IRF_stim = nanmean(stimmat);
            zIRF_stim(m,:) = (IRF_stim - nanmean(IRF_stim))/nanstd(IRF_stim);
        end             
        %%% calculate the beta weights separately for stim and null
        %%% trials for each reward condition via linear regression %%%
        for n = 0:1
            if n == 0
                vars = stim;
                z_IRF = zIRF_stim;
            else
                vars = null;
                z_IRF = zIRF;
            end
            for r = 1:sum(~isnan(R(:,ii)))
                clear idxtemp
                idxtemp = find(trials{ii}(r,:) == n );
                if rewardType{r} == 'L'
                    ZIRF = z_IRF(1,:);
                else 
                    ZIRF = z_IRF(2,:);
                end
                for t = 1:size(vars{r},1)
                    Y = vars{r}(t,:);
                    betaVal = calculateBeta(Y, ZIRF);
                    beta{ii}(r,idxtemp(t)) = betaVal;
                end
            end
        end
   else
       beta{ii}=[];
   end
end
out.betasAllSubjects = beta;

% combine and categorize all the variables for statistical analysis
subj = [];
betas = [];
TL = [];
RWD = [];
CT = [];
SF = [];

for ii = 1:length(beta)
    if ~isempty(beta{ii})
        subj = [subj;repelem(ii,size(reshape(trials{ii}',[],1),1))'];
        betas = [betas;reshape(beta{ii}',[],1)];
        TL = [TL;reshape(trials{ii}',[],1)];
        RWD = [RWD;repelem(PupilData{ii,4}',17)];
        CT = [CT;reshape(contrast{ii}',[],1)];
        SF = [SF;reshape(spatFreq{ii}',[],1)];
    end
end

% average beta for each subject based on trial type and reward condition 
clear betaSubjAvg
subjnum = unique(subj);
rwdLevel = {'L','H'}; 
for n = 0:1  
    for rwd = 1:2 
        for ii = 1:numel(unique(subj))
            betaSubjAvg(ii,(n*2)+rwd) = nanmean(betas(subj == subjnum(ii) & TL == n & cellfun(@(x) ismember(rwdLevel(rwd),x),RWD)));
        end
    end
end
out.betasAverage = betaSubjAvg;

% t-test
for n = 0:1
    clear x y
    x = betaSubjAvg(:,(n*2)+1);
    y = betaSubjAvg(:,(n*2)+2);
    [h,p,ci,stats] = ttest(x,y);
    if n == 0
        out.stim.beta_ttest = {p;ci;stats};
    else
        out.null.beta_ttest = {p;ci;stats};
    end
end

% difference in beta between reward conditions (high reward beta - low reward beta) 
meanBeta = betaSubjAvg;
for n = 0:1
    actualdiff(n+1) = nanmean(meanBeta(:,(n*2)+2)-meanBeta(:,(n*2)+1));
end
out.actualDifference = actualdiff;

%%% graph %%%
figure
clf 
for n = 0:1
    trialTitle = {'Stim' 'Null'};
    if n == 0
        subplot1 = subplot(1,2,2);
    else
        subplot1 = subplot(1,2,1);
    end
    for r = 1:size(meanBeta,1)
        if meanBeta(r,(n*2)+1) - meanBeta(r,(n*2)+2) < 0
            p = plot(meanBeta(r,(n*2)+1:(n*2)+2),'-','LineWidth',0.5,'Marker','o','MarkerSize',3,'MarkerFaceColor',[84/256,125/256,207/256],'Color',[84/256,125/256,207/256]);
            ax=gca;
            ax.FontSize = 28;

            hold on
        else
            p = plot(meanBeta(r,(n*2)+1:(n*2)+2),'-','LineWidth',0.5,'Marker','o','MarkerSize',3,'MarkerFaceColor',[204/256,78/256,102/256],'Color',[204/256,78/256,102/256]);
            ax=gca;
            ax.FontSize = 28;

            hold on
        end
        meanAll = nanmean(meanBeta);
    end
    b = bar(1:2,meanAll(1,(n*2)+1:(n*2)+2),0.5);
    b.FaceColor = 'flat';
    b.CData(2,:) = [204/256,78/256,102/256];
    b.CData(1,:) = [84/256,125/256,207/256];
    
    for rwd = 1:2
        clear err
        err = (std(meanBeta(:,(n*2)+rwd),[],1)/sqrt(size(meanBeta(:,(n*2)+rwd),1)));
        er = errorbar(rwd,meanAll(1,(n*2)+rwd),err);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';
    end
    xticks([1 2])
    xlim([0.5 2.5])
    ylim([0 1.2])
    xticklabels({'Low','High'})
    title(trialTitle{n+1},'FontSize',28)
    xlabel('Reward Level','FontSize',28)
    ylabel('Response Amplitude','FontSize',28)
    uistack(b,'bottom')
end


%% Accessing the effect of stimulus properties and reward levels on beta

% two-way (and three-way) ANOVA to access the effect of stimulus properties
% and reward on beta 
clear SFstim CTstim betas_stim RWDstim subjstim
SFstim = SF(TL == 0 & ~isnan(betas));
CTstim = CT(TL == 0 & ~isnan(betas));
betas_stim = betas(TL == 0 & ~isnan(betas)); %
RWDstim = RWD(TL == 0 & ~isnan(betas));
subjstim = subj(TL == 0 & ~isnan(betas));

[pvals,tbl,stats] = anovan(betas_stim,{SFstim CTstim RWDstim subjstim},'model',3,'random',4,'varnames',{'spatialFrequency','contrast','reward','subject'});
out.stim.beta_SFxCTxRWD = {pvals;tbl;stats};
[pvals,tbl,stats] = anovan(betas_stim,{CTstim RWDstim subjstim},'model',2,'random',3,'varnames',{'contrast','reward','subject'});
out.stim.beta_CTxRWD = {pvals;tbl;stats};
[pvals,tbl,stats] = anovan(betas_stim,{SFstim RWDstim subjstim},'model',2,'random',3,'varnames',{'spatialFrequency','reward','subject'});
out.stim.beta_SFxCT = {pvals;tbl;stats};

%%% graphs (two levels contrast) %%%
% response amplitude as a function of contrast for different reward
% conditions. 
figure
clf 
subplot(1,3,1)
linecolor(1,:) = [84/256, 125/256, 206/256];
linecolor(2,:) = [203/256, 78/256, 102/256];
for rwd = 1:2
    clear stdError graphline
    for c = 1:2     
        for ii = 1:numel(unique(subj))
            contrastRP(ii,(rwd-1)*2+c) = nanmean(betas(subj == subjnum(ii) & ...
                CT == c & TL == 0 & cellfun(@(x) ismember(rwdLevel(rwd),x),RWD)));    
        end   
    end
    avgCRP = mean(contrastRP);
    graphline = avgCRP((rwd-1)*2+1:(rwd-1)*2+2);
    p2 = plot([1.15 1.55],graphline,'color',linecolor(rwd,:),'Marker','.','MarkerSize',30,'LineWidth',1.5);
    hold on
    stdError = std(contrastRP,0,1,'omitnan')/sqrt(size(contrastRP,1));
    dsErrorsurface([1.15 1.55], graphline, stdError(:,(rwd-1)*2+1:(rwd-1)*2+2), p2.Color, 0.15)  
end
xticks([1.15 1.55])
xlim([1.1 1.6])
ylim([0 1])
xticklabels({'Low','High'})
xlabel('Contrast Level','FontSize',20)
ylabel('Response Amplitude','FontSize',20)
legend('low reward','','high reward')

% response amplitude as a function of spatial frequency for different reward
% conditions. 
subplot(1,3,2)
 linecolor(1,:) = [84/256, 125/256, 206/256];
 linecolor(2,:) = [203/256, 78/256, 102/256];
for rwd = 1:2
    clear stdError graphline
    for sf = 1:5    
        for ii = 1:numel(unique(subj))
            sfRP(ii,(rwd-1)*5+sf) = nanmean(betas(subj == subjnum(ii) & ...
                TL == 0 & SF == sf & cellfun(@(x) ismember(rwdLevel(rwd),x),RWD)));    
        end   
    end
    avgsfRP = mean(sfRP);
    graphline = avgsfRP((rwd-1)*5+1:(rwd-1)*5+5);
    p2 = plot(graphline,'color',linecolor(rwd,:),'Marker','.','MarkerSize',30,'LineWidth',1.5);
    hold on
    stdError = std(sfRP,0,1,'omitnan')/sqrt(size(sfRP,1));
    dsErrorsurface(1:length(graphline), graphline, stdError(:,(rwd-1)*5+1:(rwd-1)*5+5), p2.Color, 0.15)  
end
xticks([1 2 3 4 5])
xlim([0.5 5.5])
ylim([0 1])
xticklabels({0.5,0.8,1.3,2.0,3.2})
xlabel('Spatial Frequency (cpd)','FontSize',20)
ylabel('Response Amplitude','FontSize',20)
legend('low reward','','high reward')

% response amplitude as a function of spatial frequency for different contrast
% conditions.
subplot(1,3,3)
linecolor(1,:) = [236/256, 177/256, 33/256];
linecolor(2,:) = [127/256, 47/256, 141/256];
clear sfRP avgsfRP
for c = 1:2
    clear stdError graphline
    for sf = 1:5    
        for ii = 1:numel(unique(subj))
            sfRP(ii,(c-1)*5+sf) = nanmean(betas(subj == subjnum(ii) & ...
                CT == c & TL == 0 & SF == sf));    
        end   
    end
    avgsfRP = mean(sfRP);
    graphline = avgsfRP((c-1)*5+1:(c-1)*5+5);
    p2 = plot(graphline,'color',linecolor(c,:),'Marker','.','MarkerSize',30,'LineWidth',1.5);
    hold on
    stdError = std(sfRP,0,1,'omitnan')/sqrt(size(sfRP,1));
    dsErrorsurface(1:length(graphline), graphline, stdError(:,(c-1)*5+1:(c-1)*5+5), p2.Color, 0.15)  
end
xticks([1 2 3 4 5])
xlim([0.5 5.5])
ylim([0 1])
xticklabels({0.5,0.8,1.3,2.0,3.2})
xlabel('Spatial Frequency (cpd)','FontSize',20)
ylabel('Response Amplitude','FontSize',20)
legend('low contrast','','high contrast')

%% Bootstrapping 
% validates the regression analysis used to calculate response amplitudes 
% that it did not skew the results, diminishing any true effect of reward
% on the stimulus-evoked response.

if nargin == 2
    choice = input('Do you want to run bootstrap analysis? Y/N: ', 's');
    if strcmp(choice, 'Y')
        boots = 1;
    else
        boots = 0;
    end
end

if boots == 0
     warning('Bootstrap option not chosen. Ending analysis')
elseif boots == 1
    numboots = 2;
    diffmean = bootstrapPA(PupilData,R,weird,numboots);
    
    % histogram for null distributions with actual difference overlapped       
    figure
    clf
    clear actualdiff 
    meanBeta = betaSubjAvg;
    trialTitle = {'Pseudo Stim' 'Pseudo Null'};

    for n = 0:1
        subplot(1,2,n+1)
        hist(diffmean(:,n+1),44)
        hold on
        actualdiff = nanmean(meanBeta(:,(n*2)+2)-meanBeta(:,(n*2)+1)); %
        xline(actualdiff,'LineWidth',1,'Color','r');
        aboveActual =(sum(diffmean(:,n+1)>=actualdiff)/numel(diffmean(:,n+1)))*100;
        txt = [sprintf('%.2f',aboveActual),'%'];
        text(actualdiff + 0.01,620,txt,'FontSize',15)
        belowActual =100 - aboveActual;
        txt2 = [sprintf('%.2f',belowActual),'%'];
        text(actualdiff - 0.03,620,txt2,'FontSize',15)
        %xlim([-0.15 0.15])
        xlabel('High reward - low reward (z-scored)')
        %ylim([0 800])
        ylabel('Number of data')
        title(trialTitle{n+1})
    end
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


