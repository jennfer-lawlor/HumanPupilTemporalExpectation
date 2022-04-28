function Human_Block_Pupil()
Rootaddress = 'home/jennifer/Dropbox/Pilote12/'%'C:/Users/jlawlor3/Dropbox/Pilote12/';%
load([Rootaddress 'ElectrophyStruct/20_07_03_DataBehaviorGLMlatestPsychophysicsPupilHuman.mat'])
load([Rootaddress 'ElectrophyStruct/20_07_03_DataGLMlatestPsychophysicsPupilHuman.mat'])

% load('/home/jennifer/Dropbox/Pilote12/ElectrophyStruct/20_07_03_DataGLMlatestPsychophysicsPupilHuman.mat')
% load '/home/jennifer/Dropbox/Pilote12/ElectrophyStruct/20_07_03_DataBehaviorGLMlatestPsychophysicsPupilHuman.mat'
% load('C:\Users\jlawlor3\Desktop\PsychophysicsJen\Pilote12\ElectrophyStruct\DataBehaviorGLMlatestPsychophysicsPupilHuman_190209.mat')
% load 'C:\Users\jlawlor3\Desktop\PsychophysicsJen\Pilote12\ElectrophyStruct\DataGLMlatestPsychophysicsPupilHuman_190209.mat'

% general parameters
BinSize = data(1).BinSize;
ListSubject = [2:9 12:18]; % reointroduced subject 1
BlockValue =[0,1.5,3];
RemoveMotor = 0; % was 32
BaselineLength = round(0.79/BinSize); % 800ms of baseline
baselinewindow = 10;
xwindowslope = [33:66]+baselinewindow;%43:76;% 50 bins is s +10 bins for baseline until 3s That's 1.2s of signal
% general tools
coutcome = {[0.8,0.1,0.4],[0.3,0.8,0.4],[0.1,0.4,0.8]};
colorblockid = cellfun(@(x) x/255,{[102,178,255],[0,128,255],[0,76,153]},'UniformOutput',0);
yvalues =[-0.2,0.8];
xvalues = [-0.3+BinSize:BinSize:(250*BinSize)-0.3]; % Up to 5.2s
colordiff = cellfun(@(x) x/255,{[255,153,51],[255,108,0],[255,70,0]},'UniformOutput',0);

% process pupil data
for SubjectNum = ListSubject
    clear StartTrial OnsetTrial weirdpointsline weirdpointscol MatRawPupiltmp MatRawPupilChange MatRawPupilFa RealBaselineLength BaselineSignal weirdpoints TrialLength MatRawPupil vectorstart weirdpointsline weirdpointscol allrawpupil stdthreshold
    % behavioral and stim events
    StartTrial = find(data(SubjectNum).t == 0);
    OnsetTrial = find(data(SubjectNum).s(1,:) == 1);
    RealBaselineLength = unique(OnsetTrial-StartTrial);
    ListTrials{SubjectNum} = unique(trial(SubjectNum).TrialId);
    DiffLvl{SubjectNum} = trial(SubjectNum).Diff(OnsetTrial);
    Outcome{SubjectNum} = trial(SubjectNum).Outcome(OnsetTrial);
    BlockId{SubjectNum} = trial(SubjectNum).ChangeTimeBlock(ListTrials{SubjectNum});
    % intialise ButtonTime and change time
    ButtonTime{SubjectNum} = nan(1,length(ListTrials{SubjectNum}));
    ChangeTime{SubjectNum} = nan(1,length(ListTrials{SubjectNum}));
    ButtonTime{SubjectNum}(find(Outcome{SubjectNum} == 2)) = find(data(SubjectNum).s(5,ismember(trial(SubjectNum).TrialId,ListTrials{SubjectNum}(find(Outcome{SubjectNum} == 2))))==1)-find(data(SubjectNum).s(1,ismember(trial(SubjectNum).TrialId,ListTrials{SubjectNum}(find(Outcome{SubjectNum} == 2))))==1);
    ButtonTime{SubjectNum}(find(Outcome{SubjectNum}==1)) = find(data(SubjectNum).s(5,ismember(trial(SubjectNum).TrialId,ListTrials{SubjectNum}(find(Outcome{SubjectNum} == 1))))==1)-find(data(SubjectNum).s(1,ismember(trial(SubjectNum).TrialId,ListTrials{SubjectNum}(find(Outcome{SubjectNum} == 1))))==1);
    ChangeTime{SubjectNum}(ismember(ListTrials{SubjectNum},trial(SubjectNum).TrialId(find(data(SubjectNum).s(2,:)==1)))) = find(data(SubjectNum).s(2,:)==1) - OnsetTrial(1,ismember(ListTrials{SubjectNum},trial(SubjectNum).TrialId(find(data(SubjectNum).s(2,:)==1))));
    ChangeSize{SubjectNum} = nan(1,length(ListTrials{SubjectNum}));
    ChangeSize{SubjectNum} = trial(SubjectNum).Diff(OnsetTrial);
    
    
    % correct pupil size for megative value 
%     if min(data(SubjectNum).r)<0
%     data(SubjectNum).r = data(SubjectNum).r-min(data(SubjectNum).r);
%     end
    
    BaselineStartIndex = zeros(1,length(data(SubjectNum).r));
    BaselineStartIndex(OnsetTrial-BaselineLength) = 1; % should be the same as data(SubjectNum).t
    BaselineSignal = GetSnippets(data(SubjectNum).r,BaselineStartIndex,1,BaselineLength);
    
    
    % Remove trash trials from trial list
    BaselineTrial{SubjectNum} = nanmean(BaselineSignal,2)';
    StdBaselineTrial{SubjectNum} = nanstd(BaselineSignal');
    
    TrialTrash{SubjectNum} = unique([find(sum(~isnan(BaselineSignal)')<2)]);%, find(StdBaselineTrial{SubjectNum} == 0)]);%, find(StdBaselineTrial{SubjectNum}<1)]);%prctile(StdBaselineTrial{SubjectNum},2.5)),find(StdBaselineTrial{SubjectNum}>prctile(StdBaselineTrial{SubjectNum},97.5))]); % at least 5 points for the baseline, 150ms
    
    
    CorrectTrials{SubjectNum} = 1:length(StartTrial);%ListTrials{SubjectNum};
    CorrectTrials{SubjectNum}(TrialTrash{SubjectNum}) = [];
    % Build pupil evoked response matrix from sound onset -300ms to button
    % press or end of trial +1s
    % First build raw pupil response matrix
    TrialLength = find(data(SubjectNum).s(5,:)==1)-(StartTrial+(RealBaselineLength-baselinewindow));
    vectorstart =  zeros(1,length(data(SubjectNum).r));
    vectorstart(OnsetTrial-baselinewindow) = 1;
    %Intialisprint(gcf, '-dsvg', ['/home/jennifer/Desktop/PupilPaper/FigAttempt/Fig1/HitRate_ChangeTime_BlockId']);
    % e matevoked matrice to all trials are the same length
    MatRawPupiltmp = nan(length(ListTrials{SubjectNum}),max(TrialLength));
    
    % put all data(Subject).r as positive becaus ethe negative value don't
    % mean anything as an output, so shoft distrubtion by the minimum value
    MatRawPupiltmp = GetSnippets(data(SubjectNum).r,vectorstart,1,max(TrialLength));
    
    MatRawPupil = MatRawPupiltmp;
    
    %     for tt = 1:length(StartTrial)
    %         clear F TF
    %         [F,TF] = fillmissing(MatRawPupiltmp(tt,:),'spline','SamplePoints',1:size(MatRawPupiltmp(tt,:),2));
    %         MatRawPupil(tt,:) = F;
    %     end
    % Put Nan when trial ends
    for tt = 1:length(StartTrial)
        MatRawPupil(tt,TrialLength(tt):end) = nan;
        if ~ismember(tt,CorrectTrials{SubjectNum})
            MatRawPupil(tt,1:end) = nan;
        end
    end
    matrawsubject{SubjectNum} = MatRawPupil;
        
    % first remove baseline per trial
    MatEvokedPupilNorm{SubjectNum} = (MatRawPupil-BaselineTrial{SubjectNum}');%./StdBaselineTrial{SubjectNum}';
    % plot distribution to check
%     figure(1); hold on;
%     subplot(max(ListSubject),1,SubjectNum); hold on; hist(reshape(MatEvokedPupilNorm{SubjectNum},1,size(MatEvokedPupilNorm{SubjectNum},1)*size(MatEvokedPupilNorm{SubjectNum},2)),100)
    % then zscore per subject (one value for the mean one value for the std)
    MatEvokedPupilNorm{SubjectNum} = MatEvokedPupilNorm{SubjectNum}./(nanstd(reshape(MatEvokedPupilNorm{SubjectNum},1,size(MatEvokedPupilNorm{SubjectNum},1)*size(MatEvokedPupilNorm{SubjectNum},2))));%  %-min(reshape(MatEvokedPupilNorm{SubjectNum},1,size(MatEvokedPupilNorm{SubjectNum},1)*size(MatEvokedPupilNorm{SubjectNum},2))))./(max(reshape(MatEvokedPupilNorm{SubjectNum},1,size(MatEvokedPupilNorm{SubjectNum},1)*size(MatEvokedPupilNorm{SubjectNum},2))))%./nanstd(reshape(MatEvokedPupilNorm{SubjectNum},1,size(MatEvokedPupilNorm{SubjectNum},1)*size(MatEvokedPupilNorm{SubjectNum},2)));
%     figure(2); hold on;
%     subplot(max(ListSubject),1,SubjectNum); hold on; hist(reshape(MatEvokedPupilNorm{SubjectNum},1,size(MatEvokedPupilNorm{SubjectNum},1)*size(MatEvokedPupilNorm{SubjectNum},2)),100)
    
    %    -nanmean(reshape(MatEvokedPupilNorm{SubjectNum},1,size(MatEvokedPupilNorm{SubjectNum},1)*size(MatEvokedPupilNorm{SubjectNum},2)))

%MatEvokedPupilNorm{SubjectNum} = (MatEvokedPupilNorm{SubjectNum}-nanmean(MatEvokedPupilNorm{SubjectNum}(:,1:10),2))
    
    
    MatEvokedWithMotor{SubjectNum} = MatEvokedPupilNorm{SubjectNum}; % Change for hit amd ;iss and bp for Fa
   
   
    
    % Remove 1.2s before button press
    
    for iii = find(Outcome{SubjectNum} == 1)
        if (ButtonTime{SubjectNum}(iii)>RemoveMotor)
            MatEvokedPupilNorm{SubjectNum}(iii,ButtonTime{SubjectNum}(iii)-RemoveMotor:end) = nan;
        else
            MatEvokedPupilNorm{SubjectNum}(iii,1:end) = nan;
        end
    end
    for iii = find(ismember(ListTrials{SubjectNum},trial(SubjectNum).TrialId(find(data(SubjectNum).s(2,:)==1))))
        if ChangeTime{SubjectNum}(iii)==0
            MatEvokedPupilNorm{SubjectNum}(iii,1:end) = nan;
        else
            MatEvokedPupilNorm{SubjectNum}(iii,ChangeTime{SubjectNum}(iii):end) = nan;
        end
    end
    
    
    % save baseline norm for baseline analysis
    BaselineTrialNorm{SubjectNum} = (BaselineTrial{SubjectNum}-nanmean(reshape(MatRawPupil,1,size(MatRawPupil,1)*size(MatRawPupil,2))))./(nanstd(reshape(MatRawPupil,1,size(MatRawPupil,1)*size(MatRawPupil,2))));
    
end

% Plot Outcome mean
MeanEvokedResponseHit = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z) x(intersect(find(y == 2),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),Outcome(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
MeanEvokedResponseMiss = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z) x(intersect(find(y == 3),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),Outcome(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
MeanEvokedResponseFa = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z) x(intersect(find(y == 1),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),Outcome(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');

h = figure(); hold on;
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 14);
A = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseHit),(nanstd(MeanEvokedResponseHit)/sqrt(length(ListSubject))),{'Color',coutcome{2},'LineWidth',1.5},1); ...
    B = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseMiss),(nanstd(MeanEvokedResponseMiss)/sqrt(length(ListSubject))),{'Color',coutcome{3},'LineWidth',1.5},1); ...
    C = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseFa),(nanstd(MeanEvokedResponseFa)/sqrt(length(ListSubject))),{'Color',coutcome{1},'LineWidth',1.5},1); ...
    xlim([xvalues(1),xvalues(end)]); ylim([yvalues(1),yvalues(2)]);
plot([0,0],[yvalues(1),yvalues(2)],'--k');

% title('Pupil evoked response per ouctome all trials');
legend([A.mainLine,B.mainLine,C.mainLine],'Hit','Miss','False Alarm'); legend boxoff;
xlabel('Time within trial prior to change [s]');ylabel('Pupil size [corr.]');
% saveas(h,[Rootaddress 'PupilPaper/EvokedOutcomeAllTrials.pdf'])
% saveas(h,[Rootaddress 'PupilPaper/EvokedOutcomeAllTrials.fig'])


% Plot evoked response per block
MeanEvokedResponseEarlyhit = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z,w) x(intersect(find(y == 0),intersect(z,find(w == 2))),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
MeanEvokedResponseInterhit = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z,w) x(intersect(find(y == 1.5),intersect(z,find(w == 2))),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
MeanEvokedResponseLatehit = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z,w) x(intersect(find(y == 3),intersect(z,find(w == 2))),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject),'UniformOutput',false)),'UniformOutput',false)');

MeanEvokedResponseEarly = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z) x(intersect(find(y == 0),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
MeanEvokedResponseInter = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z) x(intersect(find(y == 1.5),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
MeanEvokedResponseLate = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z) x(intersect(find(y == 3),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');

MeanEvokedResponseEarly14 = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z,w) x(intersect(intersect(find(y == 0),[find(y == 0,1):find(y == 0,1)+29]),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
MeanEvokedResponseInter14 = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z,w) x(intersect(intersect(find(y == 1.5),[find(y == 1.5,1):find(y == 1.5,1)+29]),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
MeanEvokedResponseLate14 = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z,w) x(intersect(intersect(find(y == 3),[find(y == 3,1):find(y == 3,1)+29]),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');

VarEvokedResponseEarly14 = cell2mat(cellfun(@nanvar,(cellfun(@(x,y,z) x(intersect(intersect(find(y == 0),[find(y == 0,1):find(y == 0,1)+29]),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
VarEvokedResponseInter14 = cell2mat(cellfun(@nanvar,(cellfun(@(x,y,z) x(intersect(intersect(find(y == 1.5),[find(y == 1.5,1):find(y == 1.5,1)+29]),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
VarEvokedResponseLate14 = cell2mat(cellfun(@nanvar,(cellfun(@(x,y,z) x(intersect(intersect(find(y == 3),[find(y == 3,1):find(y == 3,1)+29]),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');

MeanEvokedResponseEarly34 = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z) x(intersect(intersect(find(y == 0),[find(y == 0,1)+30:find(y == 0,1)+119]),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
MeanEvokedResponseInter34 = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z) x(intersect(intersect(find(y == 1.5),[find(y == 1.5,1)+30:find(y == 1.5,1)+119]),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
MeanEvokedResponseLate34 = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z) x(intersect(intersect(find(y == 3),[find(y == 3,1)+30:find(y == 3,1)+119]),z),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),'UniformOutput',false)),'UniformOutput',false)');


VarEvokedResponseEarly34 = cell2mat(cellfun(@(w,x,y,z) nanvar(w(intersect(intersect(find(x==0),find(y>100 & y<266)),z),1:length(xvalues))),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),ChangeTime(ListSubject),CorrectTrials(ListSubject),'UniformOutput',0)');
VarEvokedResponseInter34 = cell2mat(cellfun(@(w,x,y,z) nanvar(w(intersect(intersect(find(x==1.5),find(y>100 & y<266)),z),1:length(xvalues))),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),ChangeTime(ListSubject),CorrectTrials(ListSubject),'UniformOutput',0)');
VarEvokedResponseLate34 =  cell2mat(cellfun(@(w,x,y,z) nanvar(w(intersect(intersect(find(x==3),find(y>100 & y<266)),z),1:length(xvalues))),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),ChangeTime(ListSubject),CorrectTrials(ListSubject),'UniformOutput',0)');

% this pannel A figure 2/3
% plot pupil as a function of time for hits per block
h = figure(); hold on;
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);
A = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseEarlyhit),(nanstd(MeanEvokedResponseEarlyhit)/sqrt(length(ListSubject))),{'Color',colorblockid{1},'LineWidth',1.5},1); ...
    B = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseInterhit),(nanstd(MeanEvokedResponseInterhit)/sqrt(length(ListSubject))),{'Color',colorblockid{2},'LineWidth',1.5},1); ...
    C = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseLatehit),(nanstd(MeanEvokedResponseLatehit)/sqrt(length(ListSubject))),{'Color',colorblockid{3},'LineWidth',1.5},1); ...
    xlim([xvalues(1),xvalues(length(xvalues))]); ylim([yvalues(1),yvalues(2)]);
title('Pupil evoked response per block for hits');
legend([A.mainLine,B.mainLine,C.mainLine],'Early','Interm.','Late'); legend boxoff;
xlabel('Time within trial prior to change [s]');ylabel('Pupil size [corr.]');
plot([0,0],[yvalues(1),yvalues(2)],'--k');
legend([A.mainLine,B.mainLine,C.mainLine],'Early','Interm.','Late'); legend boxoff;
% saveas(h,[Rootaddress 'PupilPaper/EvokedBlockHits.pdf'])
% saveas(h,[Rootaddress 'PupilPaper/EvokedBlockHits.fig'])


% plot pupil as a function of time for all trials per block
h = figure(); hold on;
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);
A = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseEarly),(nanstd(MeanEvokedResponseEarly)/sqrt(length(ListSubject))),{'Color',colorblockid{1},'LineWidth',1.5},1); ...
    B = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseInter),(nanstd(MeanEvokedResponseInter)/sqrt(length(ListSubject))),{'Color',colorblockid{2},'LineWidth',1.5},1); ...
    C = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseLate),(nanstd(MeanEvokedResponseLate)/sqrt(length(ListSubject))),{'Color',colorblockid{3},'LineWidth',1.5},1); ...
    xlim([xvalues(1),xvalues(length(xvalues))]); ylim([yvalues(1),yvalues(2)]);
title('Pupil evoked response per block for all trials');
xlabel('Time within trial prior to change [s]');ylabel('Pupil size [corr.]');
plot([0,0],[yvalues(1),yvalues(2)],'--k');
legend([A.mainLine,B.mainLine,C.mainLine],'Early','Interm.','Late'); legend boxoff;
% saveas(h,[Rootaddress 'PupilPaper/EvokedBlockAllTrials.pdf'])
% saveas(h,[Rootaddress 'PupilPaper/EvokedBlockAllTrials.fig'])
% 
% variable2analysis = 'EvokedBlockAllTrials'; VAR(:,:,1) = MeanEvokedResponseEarly(:,1:200);
% VAR(:,:,2) = MeanEvokedResponseInter(:,1:200); VAR(:,:,3) = MeanEvokedResponseLate;
% 
% t = table(squeeze(VAR(1,1,:)),squeeze(VAR(1,2,:)),squeeze(VAR(1,3,:)),...
%     squeeze(VAR(2,1,:)),squeeze(VAR(2,2,:)),squeeze(VAR(2,3,:)),...
%     squeeze(VAR(3,1,:)),squeeze(VAR(3,2,:)),squeeze(VAR(3,3,:)),...
%     'VariableNames',{'Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8','Y9'}); % Create a table storing the responses
% within = table(categorical({'1';'2';'3';'1';'2';'3';'1';'2';'3'}),...
%     categorical({'1';'1';'1',;'2';'2';'2';'3';'3';'3'}),'VariableNames',{'block','time'});
% rm = fitrm(t,'Y1-Y9~1','WithinDesign',within);
% ranovaTable = ranova(rm,'WithinModel','block+time');
% pVblock = ranovaTable.pValue(3);
% pVtime = ranovaTable.pValue(5);
% disp(variable2analysis);
% fprintf('p_{block}=%3.3f / p_{time}=%3.3f\n',pVblock,pVtime);
% 
% multcompare(rm,'block')
% % friedman
% meanhit=squeeze(nanmean(VAR,1))';
% [fpblock,b,cblock]=friedman(meanhit);


% Plot first quarter and the rest quarter
h = figure(); hold on;
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);
subplot(1,2,1); hold on;
A = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseEarly14),(nanstd(MeanEvokedResponseEarly14)/sqrt(length(ListSubject))),{'Color',colorblockid{1},'LineWidth',1.5},1); ...
    B = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseInter14),(nanstd(MeanEvokedResponseInter14)/sqrt(length(ListSubject))),{'Color',colorblockid{2},'LineWidth',1.5},1); ...
    C = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseLate14),(nanstd(MeanEvokedResponseLate14)/sqrt(length(ListSubject))),{'Color',colorblockid{3},'LineWidth',1.5},1); ...
    xlim([xvalues(1),xvalues(length(xvalues))]); ylim([yvalues(1),yvalues(2)]);
title('1/4 all trials');
xlabel('Time within trial prior to change [s]');ylabel('Pupil size [corr.]');
plot([0,0],[yvalues(1),yvalues(2)],'--k');
subplot(1,2,2); hold on;
A = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseEarly34),(nanstd(MeanEvokedResponseEarly34)/sqrt(length(ListSubject))),{'Color',colorblockid{1},'LineWidth',1.5},1); ...
    B = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseInter34),(nanstd(MeanEvokedResponseInter34)/sqrt(length(ListSubject))),{'Color',colorblockid{2},'LineWidth',1.5},1); ...
    C = shadedErrorBar(xvalues,nanmean(MeanEvokedResponseLate34),(nanstd(MeanEvokedResponseLate34)/sqrt(length(ListSubject))),{'Color',colorblockid{3},'LineWidth',1.5},1); ...
    xlim([xvalues(1),xvalues(length(xvalues))]); ylim([yvalues(1),yvalues(2)]);
title('3/4 all trials');
xlabel('Time within trial prior to change [s]');ylabel('Pupil size [corr.]');
plot([0,0],[yvalues(1),yvalues(2)],'--k');
legend([A.mainLine,B.mainLine,C.mainLine],'Early','Interm.','Late'); legend boxoff;
% saveas(h,[Rootaddress 'PupilPaper/EvokedBlockAllTrialsQuarter.pdf'])


% Plot first quarter and last quarter variance per block
h = figure(); hold on;
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);
subplot(1,2,1); hold on;
A = shadedErrorBar(xvalues,nanmean(VarEvokedResponseEarly14),2*(nanstd(VarEvokedResponseEarly14)/sqrt(length(ListSubject))),{'Color',colorblockid{1},'LineWidth',1.5},1); ...
    B = shadedErrorBar(xvalues,nanmean(VarEvokedResponseInter14),2*(nanstd(VarEvokedResponseInter14)/sqrt(length(ListSubject))),{'Color',colorblockid{2},'LineWidth',1.5},1); ...
    C = shadedErrorBar(xvalues,nanmean(VarEvokedResponseLate14),2*(nanstd(VarEvokedResponseLate14)/sqrt(length(ListSubject))),{'Color',colorblockid{3},'LineWidth',1.5},1); ...
    xlim([xvalues(1),xvalues(length(xvalues))]);
title('1/4 all trials variance');
xlabel('Time within trial prior to change [s]');ylabel('Pupil size [corr.]');
ylim([min(min(VarEvokedResponseEarly14)),max(max(VarEvokedResponseEarly34/2))]);
plot([0,0],[min(min(VarEvokedResponseEarly14)),max(max(VarEvokedResponseEarly34/2))],'--k');
subplot(1,2,2); hold on;
A = shadedErrorBar(xvalues,nanmean(VarEvokedResponseEarly34),2*(nanstd(VarEvokedResponseEarly34)/sqrt(length(ListSubject))),{'Color',colorblockid{1},'LineWidth',1.5},1); ...
    B = shadedErrorBar(xvalues,nanmean(VarEvokedResponseInter34),2*(nanstd(VarEvokedResponseInter34)/sqrt(length(ListSubject))),{'Color',colorblockid{2},'LineWidth',1.5},1); ...
    C = shadedErrorBar(xvalues,nanmean(VarEvokedResponseLate34),2*(nanstd(VarEvokedResponseLate34)/sqrt(length(ListSubject))),{'Color',colorblockid{3},'LineWidth',1.5},1); ...
    xlim([xvalues(1),xvalues(length(xvalues))]);
title('3/4 all trials variance');
xlabel('Time within trial prior to change [s]');ylabel('Pupil size [corr.]');
plot([0,0],[min(min(VarEvokedResponseEarly34)),max(max(VarEvokedResponseEarly34/2))],'--k');
ylim([min(min(VarEvokedResponseEarly34)),max(max(VarEvokedResponseEarly34/2))]);
legend([A.mainLine,B.mainLine,C.mainLine],'Early','Interm.','Late'); legend boxoff;
% only plot 3/4
h = figure(); hold on;
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);
A = shadedErrorBar(xvalues,nanmean(VarEvokedResponseEarly34),2*(nanstd(VarEvokedResponseEarly34)/sqrt(length(ListSubject))),{'Color',colorblockid{1},'LineWidth',1.5},1); ...
    B = shadedErrorBar(xvalues,nanmean(VarEvokedResponseInter34),2*(nanstd(VarEvokedResponseInter34)/sqrt(length(ListSubject))),{'Color',colorblockid{2},'LineWidth',1.5},1); ...
    C = shadedErrorBar(xvalues,nanmean(VarEvokedResponseLate34),2*(nanstd(VarEvokedResponseInter34)/sqrt(length(ListSubject))),{'Color',colorblockid{3},'LineWidth',1.5},1); ...
    xlim([xvalues(1),xvalues(length(xvalues))]);
title('3/4 all trials variance');
xlabel('Time within trial prior to change [s]');ylabel('Var. pupil size [corr.]');
plot([0,0],[min(min(VarEvokedResponseEarly34)),max(max(VarEvokedResponseEarly34))/2],'--k');
ylim([min(min(VarEvokedResponseEarly34)),max(max(VarEvokedResponseEarly34))/2]);
legend([A.mainLine,B.mainLine,C.mainLine],'Early','Interm.','Late'); legend boxoff;
% saveas(h,[Rootaddress '/PupilPaper/FigAttempt/VarEPRBlock34Quarter.pdf'])

% plot variance per time bin per block only for hits

for ct = 1:6
    for b = 1:3
        varianceblockct(ct,b,:) = cellfun(@(w,x,y,z) nanvar(nanmean(z(find(y(intersect(find(x == BlockValue(b)),find(w == 2)))>ct),10:floor(ct/0.03)+10),2)),Outcome(ListSubject),BlockId(ListSubject),ChangeTime(ListSubject),matrawsubject(ListSubject),'UniformOutput',true);
    end
end


for b = 1:3
    varianceblockbaseline(b,:) = cellfun(@(w,x,y,z) nanvar(nanmean(z(find(y(intersect(find(x == BlockValue(b)),find(w == 2)))>ct),1:10),2)),Outcome(ListSubject),BlockId(ListSubject),ChangeTime(ListSubject),matrawsubject(ListSubject),'UniformOutput',true);
end
% figure(); hold on; plot(vbuparianceblockbaseline(1,:)-varianceblockbaseline(3,:))
h = figure(); hold on; shadedErrorBar(1:6,squeeze(nanmean(squeeze(varianceblockct(:,1,:)./varianceblockct(:,3,:)),2))',(nanstd(squeeze(varianceblockct(:,1,:)./varianceblockct(:,3,:))')));
% saveas(h,[Rootaddress '/PupilPaper/FigAttempt/Diffvarhit.pdf'])



% plot response per block per subject
h = figure(); hold on;
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);
for snb = 1:length(ListSubject)
    subplot(5,4,snb); hold on;
    A = plot(xvalues,MeanEvokedResponseEarly(snb,:),'Color',colorblockid{1},'LineWidth',1.5); ...
        B = plot(xvalues,MeanEvokedResponseInter(snb,:),'Color',colorblockid{2},'LineWidth',1.5); ...
        C = plot(xvalues,MeanEvokedResponseLate(snb,:),'Color',colorblockid{3},'LineWidth',1.5); ...
        plot([0,0],[yvalues(1),yvalues(2)],'--k');
    xlim([xvalues(1),xvalues(length(xvalues))]); ylim([yvalues(1)-0.15,yvalues(2)+0.15]);
    xlabel('Time relative to sound onset [s]');
    ylabel('Pupil size norm.');
end


figure(); hold on; 

% linear fit per subject per trial
for snb = ListSubject
    clear removetrialswithnans interptrials transitiontrials evokedearly pcamat evokedinterm evokedlate propertrials aearly bearly ainterm binterm alate blate
    % Careful remove Fa and change time in the window
    %     removetrialswithnans = find(isnan(mean(MatEvokedPupilNorm{(snb)}(:,xwindowslope),2)))';
    % interp trials with less than 10 nans
    %     for tt = interptrials
    %         bla = MatEvokedPupilNorm{(snb)}(tt,xwindowslope);
    %     end
    % remove trials early in the block
    propertrials = 1:size(MatEvokedPupilNorm{(snb)},1);
    propertrials(unique([TrialTrash{(snb)}, unique([find(ChangeTime{(snb)}<=xwindowslope(end)),find(ButtonTime{(snb)}<=xwindowslope(end))])])) = [];
    slopemat{snb} = MatEvokedPupilNorm{(snb)}(propertrials,xwindowslope);
    
    %     slope per trial % polyfit not functional for nans
    
    slopetrial{snb} = nan(size(MatEvokedPupilNorm{(snb)},1),2);
    areatrial{snb} = nan(size(MatEvokedPupilNorm{(snb)},1),1);
    
    
    for tt = 1:length(propertrials)
        slopetrial{snb}(propertrials(tt),:)= polyfit(xwindowslope,slopemat{snb}(tt,:),1);
        areatrial{snb}(propertrials(tt),:)= (trapz(slopemat{snb}(tt,:)));
        %         figure(snb); hold on; plot(MatPupilRawSave{ListSubject(snb)}(propertrials(tt),xwindow)); plot(slopetrial{snb}(propertrials(tt),1)*xwindow+slopetrial{snb}(propertrials(tt),2));
    end
    
    % slope per trial category block and outcome
    for b = 1:3
        slopeblock(snb,b,:) = polyfit(1:length(xwindowslope),nanmean(MatEvokedPupilNorm{snb}(intersect(find(BlockId{snb} == BlockValue(b)),propertrials),xwindowslope)),1);
        areablock(snb,b)= (trapz(nanmean(MatEvokedPupilNorm{snb}(intersect(find(BlockId{snb} == BlockValue(b)),propertrials),xwindowslope))));
        areaoutcome(snb,b)= (trapz(nanmean(MatEvokedPupilNorm{snb}(intersect(find(Outcome{snb} == b),propertrials),xwindowslope))));

        for o = 1:3
            if isempty((intersect(intersect(find(BlockId{snb} == BlockValue(b)),find(Outcome{snb} == o)),propertrials))) || length((intersect(intersect(find(BlockId{snb} == BlockValue(b)),find(Outcome{snb} == o)),propertrials)))==1
                slopeblockoutcome(snb,b,o,:) = [nan,nan];
                areablockoutcome(snb,b,o)= nan;
                nbtrialsslope(snb,b,o)= nan;
            else
                slopeblockoutcome(snb,b,o,:) = polyfit(1:length(xwindowslope),nanmean(MatEvokedPupilNorm{snb}(intersect(intersect(find(BlockId{snb} == BlockValue(b)),find(Outcome{snb} == o)),propertrials),xwindowslope)),1);
                areablockoutcome(snb,b,o)= (trapz(nanmean(MatEvokedPupilNorm{snb}(intersect(intersect(find(BlockId{snb} == BlockValue(b)),find(Outcome{snb} == o)),propertrials),xwindowslope))));
                nbtrialsslope(snb,b,o)= numel(intersect(intersect(find(BlockId{snb} == BlockValue(b)),find(Outcome{snb} == o)),propertrials));
            end
        end
        subplot(6,3,snb);hold on; plot(1:size(xwindowslope,2),nanmean(MatEvokedPupilNorm{snb}(intersect(find(BlockId{snb} == BlockValue(b)),propertrials),xwindowslope)),'Color',colorblockid{b}); 
        plot(([1:size(xwindowslope,2)]*slopeblock(snb,b,1))+slopeblock(snb,b,2),'--k');
    end
end
xoutcome = [2 5 8];

xblock = [2 5 8];
subjitt = linspace(-0.1,0.1,length(ListSubject));

h =figure(); hold on;
for b = 1:3
        for snb= 1:length(ListSubject)
            plot(xblock(b)+subjitt(snb),areablock(ListSubject(snb),b),'.','Color',colorblockid{b},'Markersize',10);
        end
        plot(xblock(b),nanmean(areablock(ListSubject,b)),'.','Color',colorblockid{b},'Markersize',30);
end
xticks(xoutcome)
xticklabels({'E','I','L'})
ylabel('AUC [1-2s]')
% saveas(h,[Rootaddress 'PupilPaper/BlockAUC.pdf'])
% saveas(h,[Rootaddress 'PupilPaper/BlockAUC.fig'])% plot distribution per block all trials 


% plot distribution slope trial per block per subject per outcome
for b = 1:3
    for o = 1:3
        TrialnbOutcome (b,o,:) = (cellfun(@(y,z,w) numel(intersect(find(y == BlockValue(b)),intersect(z,find(w == o)))),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject)));
        TrialnbOutcomecorrected (b,o,:) = (cellfun(@(x,y,z,w) numel(intersect(find(~isnan(x)),intersect(find(y == BlockValue(b)),intersect(z,find(w == o))))),areatrial(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject)));
        MeanSlopeBlockOutcome(b,o,:) = (cellfun(@(x,y,z,w) nanmean(x(intersect(find(y == BlockValue(b)),intersect(z,find(w == o))),1)),slopetrial(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject)));
        MeanAreaBlockOutcome(b,o,:) = (cellfun(@(x,y,z,w) nanmean(x(intersect(find(y == BlockValue(b)),intersect(z,find(w == o))),1)),areatrial(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject)));
        TrialAreaBlockOutcome(b,o,:) = (cellfun(@(x,y,z,w) (x(intersect(find(y == BlockValue(b)),intersect(z,find(w == o))))),areatrial(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject),'UniformOutput',0));     
        TrialAreaBlockOutcomeID(b,o,:) = (cellfun(@(y,z,w) intersect(find(y == BlockValue(b)),intersect(z,find(w == o))),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject),'UniformOutput',0));     
    end
end
xblock = [-0.5 0 0.5];
xoutcome = [2 5 8];
subjitt = linspace(-0.1,0.1,length(ListSubject));

h =figure(); hold on;
for o = 1:3
    for b = 1:3
        for snb= 1:length(ListSubject)
            plot((xblock(b)+xoutcome(o))+subjitt(snb),MeanAreaBlockOutcome(b,o,snb),'.','Color',colorblockid{b},'Markersize',10);
        end
        plot((xblock(b)+xoutcome(o)),nanmean(MeanAreaBlockOutcome(b,o,:)),'.','Color',colorblockid{b},'Markersize',30);
    end
end
xticks(xoutcome)
xticklabels({'FA','Hit','Miss'})
ylabel('AUC [1-2s]')
saveas(h,[Rootaddress 'PupilPaper/BlockOutcomeAUC.pdf'])
saveas(h,[Rootaddress 'PupilPaper/BlockOutcomeAUC.fig'])


normalityareablock = swtest(reshape(areablock(ListSubject,:),1,3*15),0.05);
% test anova 1 (repeated) without outcome

variable2analysis = 'BlockAUC'; VAR = (areablock(ListSubject,:));
t = table(squeeze(VAR(:,1)),squeeze(VAR(:,2)),squeeze(VAR(:,3)),...
    'VariableNames',{'Y1','Y2','Y3'}); % Create a table storing the responses
within = table(categorical({'1';'2';'3'}),'VariableNames',{'block'});
rm = fitrm(t,'Y1-Y3~1','WithinDesign',within);
ranovaTable = ranova(rm,'WithinModel','block');
pVblock = ranovaTable.pValue(3);
disp(variable2analysis);
fprintf('p_{block}=%3.3f\n',pVblock);
multcompare(rm,'block')

% test anova 2 (repeated)

normalitytestArea = swtest(reshape(MeanAreaBlockOutcome,1,9*15),0.05);

variable2analysis = 'BlockOutcomeAUC'; VAR = (MeanAreaBlockOutcome);
t = table(squeeze(VAR(1,1,:)),squeeze(VAR(1,2,:)),squeeze(VAR(1,3,:)),...
    squeeze(VAR(2,1,:)),squeeze(VAR(2,2,:)),squeeze(VAR(2,3,:)),...
    squeeze(VAR(3,1,:)),squeeze(VAR(3,2,:)),squeeze(VAR(3,3,:)),...
    'VariableNames',{'Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8','Y9'}); % Create a table storing the responses
within = table(categorical({'1';'2';'3';'1';'2';'3';'1';'2';'3'}),...
    categorical({'1';'1';'1',;'2';'2';'2';'3';'3';'3'}),'VariableNames',{'outcome','block'});
rm = fitrm(t,'Y1-Y9~1','WithinDesign',within);
ranovaTable = ranova(rm,'WithinModel','block+outcome');
pVoutcome = ranovaTable.pValue(3);
pVblock = ranovaTable.pValue(5);
disp(variable2analysis);
fprintf('p_{diff}=%3.3f / p_{block}=%3.3f\n',pVoutcome,pVblock);

multcompare(rm,'block')


normalitytestSlope = swtest(reshape(MeanSlopeBlockOutcome,1,9*15),0.05);

h =figure(); hold on;
for o = 1:3
    for b = 1:3
        for snb= 1:length(ListSubject)
            plot((xblock(b)+xoutcome(o))+subjitt(snb),MeanSlopeBlockOutcome(b,o,snb),'.','Color',colorblockid{b},'Markersize',10);
        end
        plot((xblock(b)+xoutcome(o)),nanmean(MeanSlopeBlockOutcome(b,o,:)),'.','Color',colorblockid{b},'Markersize',30);
    end
end
xticks(xoutcome)
xticklabels({'FA','Hit','Miss'})
ylabel('Slope (linear fit 1-2s)')
saveas(h,[Rootaddress 'PupilPaper/BlockOutcomeSlope.pdf'])
variable2analysis = 'BlockOutcomeAUC'; VAR = (MeanSlopeBlockOutcome);
t = table(squeeze(VAR(1,1,:)),squeeze(VAR(1,2,:)),squeeze(VAR(1,3,:)),...
    squeeze(VAR(2,1,:)),squeeze(VAR(2,2,:)),squeeze(VAR(2,3,:)),...
    squeeze(VAR(3,1,:)),squeeze(VAR(3,2,:)),squeeze(VAR(3,3,:)),...
    'VariableNames',{'Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8','Y9'}); % Create a table storing the responses
within = table(categorical({'1';'2';'3';'1';'2';'3';'1';'2';'3'}),...
    categorical({'1';'1';'1',;'2';'2';'2';'3';'3';'3'}),'VariableNames',{'outcome','block'});
rm = fitrm(t,'Y1-Y9~1','WithinDesign',within);
ranovaTable = ranova(rm,'WithinModel','block+outcome');
pVoutcome = ranovaTable.pValue(3);
pVblock = ranovaTable.pValue(5);
disp(variable2analysis);
fprintf('p_{diff}=%3.3f / p_{block}=%3.3f\n',pVoutcome,pVblock);

% plot correlation farate and slope or AUC
load([Rootaddress '/PupilPaper/PupilBehaviorTestStruct_122820.mat'])
h = figure(); hold on;
for snb= 1:length(ListSubject)
    subplot(1,2,1); hold on; plot([mean(PupilBehaviorTestStruct.farate(:,1,snb)*100,1),mean(PupilBehaviorTestStruct.farate(:,3,snb)*100,1)],[areablock(ListSubject(snb),1),areablock(ListSubject(snb),3)],'k-');
    h1(snb) = plot(mean(PupilBehaviorTestStruct.farate(:,1,snb)*100,1),areablock(ListSubject(snb),1),'.','MarkerSize',20,'Color',colorblockid{1})
    xlabel('Mean Fa rate [%]'); ylabel('Mean Pupil Size AUC [1-3s]');
    h2(snb) = plot(mean(PupilBehaviorTestStruct.farate(:,3,snb)*100,1),areablock(ListSubject(snb),3),'.','MarkerSize',20,'Color',colorblockid{3})
end
legend([h1(1),h2(1)],{'Early','Late'},'Location','SouthEast'); legend boxoff

ratiofa = (squeeze(mean(PupilBehaviorTestStruct.farate(:,1,:),1))./squeeze(mean(PupilBehaviorTestStruct.farate(:,3,:),1)));
ratiort = squeeze((PupilBehaviorTestStruct.MedianReactionTime(1,3,:)))-squeeze(PupilBehaviorTestStruct.MedianReactionTime(3,3,:));
diffpupilsize = (areablock(ListSubject,1))-((areablock(ListSubject,3)));

% add correlation fa and AUC
linfit = fitlm([ratiofa],diffpupilsize,'linear');
valuesfa =floor(min(ratiofa)):ceil(max(ratiofa));
subplot(1,2,2); hold on
h4 = plot(valuesfa,(valuesfa.*linfit.Coefficients.Estimate(2))+linfit.Coefficients.Estimate(1),'k','LineWidth',1.5); ylim([-20,50]); 
for snb= 1:length(ListSubject)
    h3(snb) = plot(ratiofa(snb),diffpupilsize(snb),'.','MarkerSize',20,'Color',[0.7,0.7,0.7]); %text((mean(PupilBehaviorTestStruct.farate(:,1,snb),1)./mean(PupilBehaviorTestStruct.farate(:,3,snb),1))+0.01,(areablock(ListSubject(snb),1)-areablock(ListSubject(snb),3))-1,num2str(snb));
end
ylim([-20,50]); xlim([valuesfa(1),valuesfa(end)])
xlabel('Early/Late Fa rate'); ylabel('Early - Late AUC');
legend([h3(1),h4(1)],{'Indiv. subject','Linear fit'},'Location','NorthWest'); legend boxoff

 saveas(h,[Rootaddress 'PupilPaper/PupilAUCCorrFA.pdf'])

%% save data trial area for regression model subject*outcome*bloc*areapupil % nan accorond gto max trial across dimensions
DataRegPupil = nan(5400,5400,5400,5400);
subjectlisttrial = (reshape(repmat(1:15,360,1),15*360,1));
arealisttrial = nan(5400,1);
bloclisttrial = nan(5400,1);
outcomelisttrial = nan(5400,1);
triallisttrial = nan(5400,1);

for snb = 1:15
    clear reconstructblockid
    reconstructblockid = nan(length(areatrial{ListSubject(snb)}),1)
    reconstructblockid(find(BlockId{ListSubject(snb)}==0),1) = 1;
    reconstructblockid(find(BlockId{ListSubject(snb)}==1.5),1) = 2;
    reconstructblockid(find(BlockId{ListSubject(snb)}==3),1) = 3;
   
    starttrialsubj = (360*(snb-1))+1;
    arealisttrial(starttrialsubj:(starttrialsubj+length(areatrial{ListSubject(snb)}))-1,1) = areatrial{ListSubject(snb)};
%     triallisttrial(starttrialsubj:(starttrialsubj+length(areatrial{ListSubject(snb)}))-1,1) = CorrectTrials{ListSubject(snb)};
    outcomelisttrial(starttrialsubj:(starttrialsubj+length(areatrial{ListSubject(snb)}))-1,1) = Outcome{ListSubject(snb)};
    bloclisttrial(starttrialsubj:(starttrialsubj+length(areatrial{ListSubject(snb)}))-1,1) = reconstructblockid;
end
% remove trials that don't have an area
outcomelisttrial(isnan(arealisttrial),1) = nan;
bloclisttrial(isnan(arealisttrial),1) = nan;
subjectlisttrial(isnan(arealisttrial),1) = nan;
DataPupilReg = [subjectlisttrial,arealisttrial,bloclisttrial,outcomelisttrial];
save([Rootaddress 'PupilPaper/110821_DataPupilReg.mat'],'DataPupilReg');
% % for slope
% h = figure(); hold on;
% for snb= 1:length(ListSubject)
%     subplot(1,2,1); hold on; plot([mean(PupilBehaviorTestStruct.farate(:,1,snb)*100,1),mean(PupilBehaviorTestStruct.farate(:,3,snb)*100,1)],[slopeblock(ListSubject(snb),1),slopeblock(ListSubject(snb),3)],'k-');
%     h1(snb) = plot(mean(PupilBehaviorTestStruct.farate(:,1,snb)*100,1),slopeblock(ListSubject(snb),1),'.','MarkerSize',20,'Color',colorblockid{1})
%     xlabel('Mean Fa rate [%]'); ylabel('Mean Pupil Size AUC [1-3s]');
%     h2(snb) = plot(mean(PupilBehaviorTestStruct.farate(:,3,snb)*100,1),slopeblock(ListSubject(snb),3),'.','MarkerSize',20,'Color',colorblockid{3})
% end
% legend([h1(1),h2(1)],{'Early','Late'},'Location','SouthEast'); legend boxoff
% 
% ratiofa = (squeeze(mean(PupilBehaviorTestStruct.farate(:,1,:),1))./squeeze(mean(PupilBehaviorTestStruct.farate(:,3,:),1)));
% ratiort = squeeze((PupilBehaviorTestStruct.MedianReactionTime(1,3,:)))-squeeze(PupilBehaviorTestStruct.MedianReactionTime(3,3,:));
% diffpupilsize = (slopeblock(ListSubject,1))-((slopeblock(ListSubject,3)));
% 
% % add correlation fa and AUC
% linfit = fitlm([ratiofa],diffpupilsize,'linear');
% valuesfa =floor(min(ratiofa)):ceil(max(ratiofa));
% subplot(1,2,2); hold on
% h4 = plot(valuesfa,(valuesfa.*linfit.Coefficients.Estimate(2))+linfit.Coefficients.Estimate(1),'k','LineWidth',1.5); ylim([-20,50]); 
% for snb= 1:length(ListSubject)
%     h3(snb) = plot(ratiofa(snb),diffpupilsize(snb),'.','MarkerSize',20,'Color',[0.7,0.7,0.7]); %text((mean(PupilBehaviorTestStruct.farate(:,1,snb),1)./mean(PupilBehaviorTestStruct.farate(:,3,snb),1))+0.01,(areablock(ListSubject(snb),1)-areablock(ListSubject(snb),3))-1,num2str(snb));
% end
% ylim([-0.002,0.02]); xlim([valuesfa(1),valuesfa(end)])
% xlabel('Early/Late Fa rate'); ylabel('Early - Late Slope');
% legend([h3(1),h4(1)],{'Indiv. subject','Linear fit'},'Location','NorthWest'); legend boxoff
% 
%  saveas(h,[Rootaddress 'PupilPaper/PupilSlopeCorrFA.pdf'])





%  close all
% 
%  h = figure(); hold on;
% for snb= 1:length(ListSubject)
%    hold on; plot([mean(PupilBehaviorTestStruct.farate(:,1,snb)*100,1),mean(PupilBehaviorTestStruct.farate(:,3,snb)*100,1)],[areablock(ListSubject(snb),1),areablock(ListSubject(snb),3)],'k-');
%     h1(snb) = plot(mean(PupilBehaviorTestStruct.farate(:,1,snb)*100,1),areablock(ListSubject(snb),1),'.','MarkerSize',20,'Color',colorblockid{1})
%     xlabel('Mean Fa rate [%]'); ylabel('Mean pupil size AUC [1-2s]');
%     h2(snb) = plot(mean(PupilBehaviorTestStruct.farate(:,3,snb)*100,1),areablock(ListSubject(snb),3),'.','MarkerSize',20,'Color',colorblockid{3})
% end
% legend([h1(1),h2(1)],{'Early','Late'},'Location','SouthEast'); legend boxoff
% 
% % ratiofa = squeeze(mean(PupilBehaviorTestStruct.farate(:,1,:),1))./squeeze(mean(PupilBehaviorTestStruct.farate(:,3,:),1));
% % diffpupilsize = areablock(ListSubject,1)-(areablock(ListSubject,3));
%  saveas(h,[Rootaddress 'PupilPaper/PupilAUCCorrFA1.fig'])
%  
%  h = figure(); hold on;
% h4 = plot(valuesfa,(valuesfa.*linfit.Coefficients.Estimate(2))+linfit.Coefficients.Estimate(1),'k','LineWidth',1.5); ylim([-20,50]); 
% for snb= 1:length(ListSubject)
%     h3(snb) = plot(ratiofa(snb),diffpupilsize(snb),'.','MarkerSize',20,'Color',[0.7,0.7,0.7]); %text((mean(PupilBehaviorTestStruct.farate(:,1,snb),1)./mean(PupilBehaviorTestStruct.farate(:,3,snb),1))+0.01,(areablock(ListSubject(snb),1)-areablock(ListSubject(snb),3))-1,num2str(snb));
% end
% ylim([-20,30]); xlim([valuesfa(1),valuesfa(end)])
% xlabel('Early/Late Fa rate'); ylabel('Early - Late AUC');
% legend([h3(1),h4(1)],{'Indiv. subject','Linear fit'},'Location','NorthWest'); legend boxoff
% saveas(h,[Rootaddress 'PupilPaper/PupilAUCCorrFA2.fig'])
 
% add correlation fa and Slope
diffpupilsize = (slopeblock(ListSubject,1))-((slopeblock(ListSubject,3)));
h = figure(); hold on;
for snb= 1:length(ListSubject)
    subplot(1,2,1); hold on; plot([mean(PupilBehaviorTestStruct.farate(:,1,snb)*100,1),mean(PupilBehaviorTestStruct.farate(:,3,snb)*100,1)],[slopeblock(ListSubject(snb),1),slopeblock(ListSubject(snb),3)],'k-');
    h1(snb) = plot(mean(PupilBehaviorTestStruct.farate(:,1,snb)*100,1),slopeblock(ListSubject(snb),1),'.','MarkerSize',20,'Color',colorblockid{1})
    xlabel('Mean Fa rate [%]'); ylabel('Mean Pupil Size AUC [1-2s]');
    h2(snb) = plot(mean(PupilBehaviorTestStruct.farate(:,3,snb)*100,1),slopeblock(ListSubject(snb),3),'.','MarkerSize',20,'Color',colorblockid{3})
end
legend([h1(1),h2(1)],{'Early','Late'},'Location','SouthEast'); legend boxoff
linfit = fitlm([ratiofa],diffpupilsize,'linear');
valuesfa =floor(min(ratiofa)):ceil(max(ratiofa));
subplot(1,2,2); hold on
h4 = plot(valuesfa,(valuesfa.*linfit.Coefficients.Estimate(2))+linfit.Coefficients.Estimate(1),'k','LineWidth',1.5); ylim([-20,50]); 
for snb= 1:length(ListSubject)
    h3(snb) = plot(ratiofa(snb),diffpupilsize(snb),'.','MarkerSize',20,'Color',[0.7,0.7,0.7]); %text((mean(PupilBehaviorTestStruct.farate(:,1,snb),1)./mean(PupilBehaviorTestStruct.farate(:,3,snb),1))+0.01,(areablock(ListSubject(snb),1)-areablock(ListSubject(snb),3))-1,num2str(snb));
end
ylim([min(diffpupilsize),max(diffpupilsize)]); xlim([valuesfa(1),valuesfa(end)])
xlabel('Early/Late Fa rate'); ylabel('Early - Late Slope');
legend([h3(1),h4(1)],{'Indiv. subject','Linear fit'},'Location','NorthWest'); legend boxoff

 saveas(h,[Rootaddress 'PupilPaper/PupilSlopeCorrFA.pdf'])

 % do the same plot only for hits
 h = figure(); hold on;
for snb= 1:length(ListSubject)
    subplot(1,2,1); hold on; plot([mean(PupilBehaviorTestStruct.farate(:,1,snb)*100,1),mean(PupilBehaviorTestStruct.farate(:,3,snb)*100,1)],[MeanAreaBlockOutcome(1,2,snb),MeanAreaBlockOutcome(3,2,snb)],'k-');
    h1(snb) = plot(mean(PupilBehaviorTestStruct.farate(:,1,snb)*100,1),MeanAreaBlockOutcome(1,2,snb),'.','MarkerSize',20,'Color',colorblockid{1})
    xlabel('Mean Fa rate [%]'); ylabel('Mean Pupil Size AUC [1-2s] only hits');
    h2(snb) = plot(mean(PupilBehaviorTestStruct.farate(:,3,snb)*100,1),MeanAreaBlockOutcome(3,2,snb),'.','MarkerSize',20,'Color',colorblockid{3})
end
legend([h1(1),h2(1)],{'Early','Late'},'Location','SouthEast'); legend boxoff

diffpupilsizehit = squeeze(MeanAreaBlockOutcome(1,2,:))-squeeze(MeanAreaBlockOutcome(3,2,:));

 % add correlation fa and AUC 
 linfit = fitlm([ratiofa],diffpupilsizehit,'linear');
 subplot(1,2,2); hold on
h4 = plot(valuesfa,(valuesfa.*linfit.Coefficients.Estimate(2))+linfit.Coefficients.Estimate(1),'k','LineWidth',1.5); ylim([-20,50]); 
for snb= 1:length(ListSubject)
    h3(snb) = plot(ratiofa(snb),diffpupilsizehit(snb),'.','MarkerSize',20,'Color',[0.7,0.7,0.7]);
end
ylim([-20,50]); xlim([valuesfa(1),valuesfa(end)])
xlabel('Early/Late Fa rate'); ylabel('Early - Late AUC (hit only)');
legend([h3(1),h4(1)],{'Indiv. subject','Linear fit'},'Location','NorthWest'); legend boxoff

 saveas(h,[Rootaddress 'PupilPaper/PupilAUCCorrFAOnlyHit.pdf'])
 
 
 % same but for RT instead of farate
%   h = figure(); hold on;
% for snb= 1:length(ListSubject)
%     subplot(1,2,1); hold on; plot([PupilBehaviorTestStruct.MedianReactionTime(1,3,snb),PupilBehaviorTestStruct.MedianReactionTime(3,3,snb)],[MeanAreaBlockOutcome(1,2,snb),MeanAreaBlockOutcome(3,2,snb)],'k-');
%     h1(snb) = plot(PupilBehaviorTestStruct.MedianReactionTime(1,3,snb),MeanAreaBlockOutcome(1,2,snb),'.','MarkerSize',20,'Color',colorblockid{1})
%     xlabel('Median RT [s]'); ylabel('Mean Pupil Size AUC [1-3s] only hits');
%     h2(snb) = plot(PupilBehaviorTestStruct.MedianReactionTime(3,3,snb),MeanAreaBlockOutcome(3,2,snb),'.','MarkerSize',20,'Color',colorblockid{3})
% end
% legend([h1(1),h2(1)],{'Early','Late'},'Location','SouthEast'); legend boxoff
% 
%  % add correlation rt and AUC 
%  linfit = fitlm([ratiort],diffpupilsizehit,'linear');
%  subplot(1,2,2); hold on
% for snb= 1:length(ListSubject)
%     h3(snb) = plot(ratiort(snb),diffpupilsizehit(snb),'.','MarkerSize',20,'Color',[0.7,0.7,0.7]);
% end
% h4 = plot([0:1],([0:1].*linfit.Coefficients.Estimate(2))+linfit.Coefficients.Estimate(1),'k','LineWidth',1.5); ylim([-20,50]); %xlim([0,4])
% 
% ylim([-20,50]); xlim([0,4])
% xlabel('Early Median RT (large)/Late Median RT (large)'); ylabel('Early AUC - Late AUC (hit only)');
% legend([h3(1),h4(1)],{'Indiv. subject','Linear fit'},'Location','NorthWest'); legend boxoff
% 
%  saveas(h,[Rootaddress 'PupilPaper/PupilAUCCorrRTOnlyHit.pdf'])
%  
% same for criterion
%  h = figure(); hold on;
% for snb= 1:length(ListSubject)
%     subplot(1,2,1); hold on; plot([nanmean(PupilBehaviorTestStruct.criterion(snb,1,:),3),nanmean(PupilBehaviorTestStruct.criterion(snb,3,:),3)],[areablock(ListSubject(snb),1),areablock(ListSubject(snb),3)],'k-');
%     h1(snb) = plot(nanmean(PupilBehaviorTestStruct.criterion(snb,1,:),3),areablock(ListSubject(snb),1),'.','MarkerSize',20,'Color',colorblockid{1})
%     xlabel('Criterion'); ylabel('Mean Pupil Size AUC [1-3s] only hits');
%     h2(snb) = plot(nanmean(PupilBehaviorTestStruct.criterion(snb,3,:),3),areablock(ListSubject(snb),3),'.','MarkerSize',20,'Color',colorblockid{3})
% end
% legend([h1(1),h2(1)],{'Early','Late'},'Location','SouthEast'); legend boxoff
% 
% 
% ratiocriterion = squeeze(nanmean(PupilBehaviorTestStruct.criterion(:,1,:),3))./squeeze(nanmean(PupilBehaviorTestStruct.criterion(:,3,:),3));
% 
% linfit = fitlm([ratiocriterion],diffpupilsize,'linear');
% subplot(1,2,2); hold on
% h4 = plot([-2:2],([-2:2].*linfit.Coefficients.Estimate(2))+linfit.Coefficients.Estimate(1),'k','LineWidth',1.5); ylim([-20,50]); xlim([-1.5,1.5])
% for snb= 1:length(ListSubject)
%     h3(snb) = plot(ratiocriterion(snb),areablock(ListSubject(snb),1)-areablock(ListSubject(snb),3),'.','MarkerSize',20,'Color',[0.7,0.7,0.7]);
% end
% ylim([-20,50]); xlim([-0.5,1.5])
% xlabel('Early criterion-Late criterion'); ylabel('Early AUC - Late AUC');
% legend([h3(1),h4(1)],{'Indiv. subject','Linear fit'},'Location','NorthWest'); legend boxoff
% % 
%  saveas(h,[Rootaddress 'PupilPaper/PupilAUCCorrCriterion.pdf'])
 
 
h = figure(); hold on;
for snb= 1:length(ListSubject)
    subplot(1,2,1); hold on; plot([mean(PupilBehaviorTestStruct.farate(:,1,snb),1),mean(PupilBehaviorTestStruct.farate(:,3,snb),1)],[slopeblock(ListSubject(snb),1,1),slopeblock(ListSubject(snb),3,1)],'k-');
    plot(mean(PupilBehaviorTestStruct.farate(:,1,snb),1),slopeblock(ListSubject(snb),1,1),'.','MarkerSize',20,'Color',colorblockid{1})
    xlabel('Mean Fa rate [%]'); ylabel('Mean Pupil Size slope [1-3s]');
    plot(mean(PupilBehaviorTestStruct.farate(:,3,snb),1),slopeblock(ListSubject(snb),3,1),'.','MarkerSize',20,'Color',colorblockid{3})
    subplot(1,2,2); hold on; plot(mean(PupilBehaviorTestStruct.farate(:,1,snb),1)/mean(PupilBehaviorTestStruct.farate(:,3,snb),1),slopeblock(ListSubject(snb),1,1)-slopeblock(ListSubject(snb),3,1),'.','MarkerSize',20,'Color',[0.5,0.5,0.5]);
    xlabel('Early Fa rate/Late Fa rate'); ylabel('Early slope - Late slope');
 end
 saveas(h,[Rootaddress 'PupilPaper/PupilSlopeCorr.pdf'])

%   corrcoef(squeeze(mean(PupilBehaviorTestStruct.farate(:,1,:),1))./squeeze(mean(PupilBehaviorTestStruct.farate(:,3,:),1)),slopeblock(ListSubject,1,1)-slopeblock(ListSubject,3,1))

% Plot variance

for b = 1:3
    for o =1:3
        Evokedvar(b,o,:) = cellfun(@(x,y,w,z,d) nanmean(nanvar(x(intersect(intersect(find(y==BlockValue(b)),find(w==o)),unique([find(z>110 & z<266),find(d>142)])),10:110))),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),Outcome(ListSubject),ChangeTime(ListSubject),ButtonTime(ListSubject));
%                 Evokedvar(b,o,:) = cellfun(@(x,y,w,z,d) nanmean(nanvar(x(intersect(intersect(find(y==BlockValue(b)),find(w==o)),unique([find(z>110 & z<266),find(d>142)])),10:110))),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),Outcome(ListSubject),ChangeTime(ListSubject),ButtonTime(ListSubject));
    end
end

h =figure(); hold on;
for b = 1:3
    for o = 1:3
        for snb= 1:length(ListSubject)
            plot((xblock(b)+xoutcome(o))+subjitt(snb),Evokedvar(b,o,snb),'.','Color',coutcome{o},'Markersize',10);
        end
        plot((xblock(b)+xoutcome(o)),nanmean(Evokedvar(b,o,:)),'.','Color',coutcome{o},'Markersize',30);
    end
end
xticks(xblock)
xticklabels({'Early','Interm.','Late'})
ylabel('Var (linear fit 1-3s)')
saveas(h,[Rootaddress 'PupilPaper/BlockOutcomeSlope.pdf'])


% % plot the same but only using subjects mean
% figure(); hold on;
% for b = 1:3
%     for o = 1:3
%         for snb= 1:length(ListSubject)
%             plot((xblock(b)+xoutcome(o))+subjitt(snb),slopeblockoutcome(snb,b,o,1),'.','Color',coutcome{o},'Markersize',10);
%         end
%         plot((xblock(b)+xoutcome(o)),nanmean(slopeblockoutcome(:,b,o,1)),'.','Color',coutcome{o},'Markersize',30);
%     end
% end
% xticks(xblock)
% xticklabels({'Early','Interm.','Late'})
% ylabel('Slope (linear fit 1-3s)')
% 
% figure(); hold on;
% for b = 1:3
%     for o = 1:3
%         for snb= 1:length(ListSubject)
%             plot((xblock(b)+xoutcome(o))+subjitt(snb),areablockoutcome(snb,b,o),'.','Color',coutcome{o},'Markersize',10);
%         end
%         plot((xblock(b)+xoutcome(o)),nanmean(areablockoutcome(:,b,o)),'.','Color',coutcome{o},'Markersize',30);
%     end
% end
% xticks(xblock)
% xticklabels({'Early','Interm.','Late'})
% ylabel('AUC (linear fit 1-3s)')


% % plot distribution slope trial per block per subject per change size
% for b = 1:3
%     for d = 1:3
%         TrialnBlockCT (b,d,:) = (cellfun(@(y,z,w) numel(intersect(find(y == BlockValue(b)),intersect(z,find(w == difflvlind(d))))),BlockId(ListSubject),CorrectTrials(ListSubject),ChangeSize(ListSubject)));
%         MeanSlopeBlockCT(b,d,:) = (cellfun(@(x,y,z,w) nanmean(x(intersect(find(y == BlockValue(b)),intersect(z,find(w == difflvlind(d)))),1)),slopetrial(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),ChangeSize(ListSubject)));
%         MeanAreaBlockCT(b,d,:) = (cellfun(@(x,y,z,w) nanmean(x(intersect(find(y == BlockValue(b)),intersect(z,find(w == difflvlind(d)))),1)),areatrial(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),ChangeSize(ListSubject)));
%     end
% end
% xoutcome = [-0.5 0 0.5];
% xblock = [2 5 8];
% subjitt = linspace(-0.1,0.1,length(ListSubject));
% figure(); hold on;
% for b = 1:3
%     for d = 1:3
%         for snb= 1:length(ListSubject)
%             plot((xblock(b)+xoutcome(d))+subjitt(snb),MeanAreaBlockCT(b,d,snb),'.','Color',colordiff{d},'Markersize',10);
%         end
%         plot((xblock(b)+xoutcome(d)),nanmean(MeanAreaBlockCT(b,d,:)),'.','Color',colordiff{d},'Markersize',30);
%     end
% end
% xticks(xblock)
% xticklabels({'Early','Interm.','Late'})
% ylabel('AUC 0-3s')
%
% figure(); hold on;
% for b = 1:3
%     for d = 1:3
%         for snb= 1:length(ListSubject)
%             plot((xblock(b)+xoutcome(d))+subjitt(snb),MeanSlopeBlockCT(b,d,snb),'.','Color',colordiff{d},'Markersize',10);
%         end
%         plot((xblock(b)+xoutcome(d)),nanmean(MeanSlopeBlockCT(b,d,:)),'.','Color',colordiff{d},'Markersize',30);
%     end
% end
% xticks(xblock)
% xticklabels({'Early','Interm.','Late'})
% ylabel('Slope (linear fit 1-3s)')

% plot hit as a function of difflvl
difflvlind = unique(DiffLvl{2});

for dd = 1:3
    HitEvokedResponseDiff{dd} = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z,w) x(intersect(find(y == difflvlind(dd)),intersect(z,find(w == 2))),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),DiffLvl(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
    changetimetarget{dd} = (cellfun(@(y,z,w,t) t(intersect(find(y == difflvlind(dd)),intersect(z,find(w == 2)))),DiffLvl(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject),ChangeTime(ListSubject),'UniformOutput',false));
    trialnb{dd}= cellfun(@(y,z,w) intersect(find(y == difflvlind(dd)),intersect(z,find(w == 2))),DiffLvl(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject),'UniformOutput',false);
    Pupilsizepriortochange{dd} = (cellfun(@(x,y,z) nanmean(x(y(find(z>30)),(z(find(z>30))+10)-30:z(find(z>30))+10),2),MatEvokedPupilNorm(ListSubject),trialnb{dd},changetimetarget{dd},'UniformOutput',false));
    xvaluesDiffCT{dd} = (cellfun(@(z) (z(find(z>100 & z<266))+10),changetimetarget{dd},'UniformOutput',false));
    trialDiffCT{dd} = (cellfun(@(z) (find(z>100 & z<266)),changetimetarget{dd},'UniformOutput',false));
    HitEvokedResponseDiffCT{dd} = cell2mat(cellfun(@nanmean,cellfun(@(x,y,z) x(y,((z-30):(z+10))),MatEvokedPupilNorm(ListSubject),trialDiffCT{dd},xvaluesDiffCT{dd},'UniformOutput',false),'UniformOutput',false)');
end
h = figure(); hold on;
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);
A = shadedErrorBar(xvalues,nanmean(HitEvokedResponseDiff{1}),(nanstd(HitEvokedResponseDiff{1})/sqrt(length(ListSubject))),{'Color',colordiff{1},'LineWidth',1.5},1); ...
    %     B = shadedErrorBar(xvalues,nanmean(HitEvokedResponseDiff{2}),(nanstd(HitEvokedResponseDiff{2})/sqrt(length(ListSubject))),{'Color',colordiff{2},'LineWidth',1.5},1); ...
C = shadedErrorBar(xvalues,nanmean(HitEvokedResponseDiff{3}),(nanstd(HitEvokedResponseDiff{3})/sqrt(length(ListSubject))),{'Color',colordiff{3},'LineWidth',1.5},1); ...
    xlim([xvalues(1),xvalues(length(xvalues))]); ylim([yvalues(1),yvalues(2)]);
title('Hit per difflvl');
xlabel('Time relative to sound onset [s]');
ylabel('Pupil size norm.');
plot([0,0],[-6,12],'--k');
legend([A.mainLine,C.mainLine],'Small','Large'); legend boxoff;


% plot size hit prior to change as a function of diff lvl and ct
h = figure(); hold on;
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);
subplot(1,3,[1,2]); hold on;
for snb = 1:15
    for dd = [1,3]
        clear sortchange sortchangeind
        A = plot(changetimetarget{dd}{snb}(changetimetarget{dd}{snb}>30)+10,(Pupilsizepriortochange{dd}{snb}),'.','Color',colordiff{dd},'LineWidth',1.5);
        [sortchange,sortchangeind] = sort(changetimetarget{dd}{snb}(changetimetarget{dd}{snb}>30));
        Matchangediff{dd}{snb}(:,1) = Pupilsizepriortochange{dd}{snb}(sortchangeind);
        Matchangediff{dd}{snb}(:,2) = changetimetarget{dd}{snb}(sortchangeind);
    end
end
title('Hit per difflvl');
xlabel('Time relative to sound onset [s]');
ylabel('Pupil size norm.');
plot([0,0],[-6,12],'--k');
% % plot diftribution easy vs hard
[ydistrihard,xdistrihard] = hist(vertcat(Pupilsizepriortochange{1}{:}),20);
[ydistrieasy,xdistrieasy] = hist(vertcat(Pupilsizepriortochange{3}{:}),20);
subplot(1,3,[3]); hold on; plot(ydistrieasy/max(ydistrieasy),xdistrieasy,'Color',colordiff{3}); plot(ydistrihard/max(ydistrihard),xdistrihard,'Color',colordiff{1})

% Plot block as a function of change time
blockind = unique(BlockId{2});
for dd = 1:3
    HitEvokedResponseblock{dd} = cell2mat(cellfun(@nanmean,(cellfun(@(x,y,z,w) x(intersect(find(y == blockind(dd)),intersect(z,find(w == 2))),1:length(xvalues)),MatEvokedPupilNorm(ListSubject),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject),'UniformOutput',false)),'UniformOutput',false)');
    changetimetargetblock{dd} = (cellfun(@(y,z,w,t) t(intersect(find(y == blockind(dd)),intersect(z,find(w == 2)))),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject),ChangeTime(ListSubject),'UniformOutput',false));
    trialnbblock{dd}= cellfun(@(y,z,w) intersect(find(y == blockind(dd)),intersect(z,find(w == 2))),BlockId(ListSubject),CorrectTrials(ListSubject),Outcome(ListSubject),'UniformOutput',false);
    Pupilsizepriortochangeblock{dd} = (cellfun(@(x,y,z) nanmean(x(y(find(z>30)),(z(find(z>30))+10)-30:z(find(z>30))+10),2),MatEvokedPupilNorm(ListSubject),trialnbblock{dd},changetimetargetblock{dd},'UniformOutput',false));
end

% h = figure(); hold on;
% set(gca, 'FontName', 'Arial');
% set(gca, 'FontSize', 12);
% subplot(1,3,[1,2]); hold on;
% for snb = 1:15
%     for dd = [1,3]
%         clear sortchange sortchangeind ctind
%         ctind = find(changetimetargetblock{dd}{snb}>100 & changetimetargetblock{dd}{snb}<=266);
%         A = plot((changetimetargetblock{dd}{snb}(ctind)*0.03),(Pupilsizepriortochangeblock{dd}{snb}(ctind)),'.','Color',colorblockid{dd},'LineWidth',1.5);
%         %         [sortchange,sortchangeind] = sort(changetimetargetblock{dd}{snb}(changetimetargetblock{dd}{snb}>30));
%         %         Matchangediff{dd}{snb}(:,1) = Pupilsizepriortochange{dd}{snb}(sortchangeind);
%         %         Matchangediff{dd}{snb}(:,2) = changetimetargetblock{dd}{snb}(sortchangeind);
%     end
% end
% xlabel('Change time [s]'); ylabel('Pupil size prior to change for hits (900ms average)')
% subplot(1,3,[3]); hold on;
% [ydistriearly,xdistriearly] = hist(vertcat(Pupilsizepriortochangeblock{1}{:}),20);
% [ydistrilate,xdistrilate] = hist(vertcat(Pupilsizepriortochangeblock{3}{:}),20);
% plot(ydistriearly/max(ydistriearly),xdistriearly,'Color',colorblockid{1}); plot(ydistrilate/max(ydistrilate),xdistrilate,'Color',colorblockid{3})
% set(gca,'YTickLabel',[]);
% legend('early','late'); legend boxoff
% saveas(h,[Rootaddress 'PupilPaper/PupilsizeDistriPriorChange.pdf'])


% % Plot hit reate as a function of block diff and change time
% bins = round((0:2:12)/0.03);
% for b = 1:3
%     for diff = 1:3
%         for bb = 1:length(bins)-1
%             for o = 1:2
%                 trialinddiffblock(b,diff,bb,o,:) = cellfun(@(x,y,w,z) intersect(intersect(find(x == difflvlind(diff)),intersect(find(y == blockind(b)),find(w == (o+1)))),intersect(find(z>bins(bb)),find(z<bins(bb+1)))),DiffLvl(ListSubject),BlockId(ListSubject),Outcome(ListSubject),ChangeTime(ListSubject),'UniformOutput',false);
%                 trialnumeldiffblock(b,diff,bb,o,:) = squeeze(cellfun(@numel,trialinddiffblock(b,diff,bb,o,:)))';
%             end
%             perfblockdiff(b,diff,bb,:) = squeeze(trialnumeldiffblock(b,diff,bb,1,:))./(squeeze(trialnumeldiffblock(b,diff,bb,1,:))+squeeze(trialnumeldiffblock(b,diff,bb,2,:)));
%         end
%     end
% end
% figure;
%
% for b = 1:3
%     subplot(1,3,b); hold on;
%     for diff = 1:3
%         A(diff) = plot(bins(2:end),squeeze(nanmean(perfblockdiff(b,diff,:,:),4))','Color',colordiff{diff},'Markersize',30);
%     end
% end
% xticks(xblock)
% xticklabels({'Early','Interm.','Late'})
% ylabel('Perf [%]')
% legend([A(1),A(2),A(3)],{'small','interm.','large'},'Location','SouthWest'); legend boxoff

% % plot correlation slope and diff effect
% diffslope = squeeze((slopeblock(ListSubject,1,1)./slopeblockoutcome(ListSubject,3,1)));
% diffarea = squeeze((areablock(ListSubject,1)./areablock(ListSubject,3)));
% diffperf = squeeze(perfblockdiff(1,1,:)-perfblockdiff(3,1,:));
% figure();
% plot(diffarea,diffperf,'.');
% figure();
% plot(diffslope,diffperf,'.');

% % RT for trial 3s to 8 s as a function of time per block
% figure(); hold on;
% for b = 1:3
%     subplot(1,3,b); hold on; title(['Block ' num2str(blockind(b))]); xlabel('Change time bin (end of bin 1s)'); ylabel('Reaction time')
%     for diff = 1:3
%         clear sortctbinind meanretallct
%         binbctind = [3:8];
%         for ctbin = 2:length(binbctind)
%             [sortctbinind{ctbin-1}] = find((ctall{b,diff}*0.03)>=binbctind(ctbin-1) &(ctall{b,diff}*0.03)<binbctind(ctbin));
%             medianretallct(b,diff,ctbin-1) = median(rtall{b,diff}(sortctbinind{ctbin-1}));
%             semrtallct(b,diff,ctbin-1) = std(rtall{b,diff}(sortctbinind{ctbin-1}))./sqrt(length((sortctbinind{ctbin-1})));
%         end
%         shadedErrorBar((binbctind(2:end)),medianretallct(b,diff,:),semrtallct(b,diff,:),{'Color',colordiff{diff}})
%     end
% end

% baseline per block

baselinemat = (cellfun(@(x) (x),BaselineTrialNorm(ListSubject),'UniformOutput',0));
for b =1:3
    baselineblock(b,:) = cellfun(@(x,y) nanmean(x(1,intersect(find(y==blockind(b)),[30:120,150:240,270:360]))),baselinemat,BlockId(ListSubject));
    varbaselineblock(b,:) = cellfun(@(x,y) nanvar(x(1,intersect(find(y==blockind(b)),[30:120,150:240,270:360]))),baselinemat,BlockId(ListSubject));
    baselineoutcome(b,:) = cellfun(@(x,y) nanmean(x(1,find(y==b))),baselinemat,Outcome(ListSubject));
end

figure(); hold on;
for snb= 1:length(ListSubject)
    for b = 1:3
        plot((xblock(:))+subjitt(snb)',baselineblock(:,snb)','Color',colorblockid{b},'Markersize',10);
        plot((xblock(b))+subjitt(snb)',baselineblock(b,snb)','.','Color',colorblockid{b},'Markersize',10);
    end
end
for b = 1:3
    G(b) = plot((xblock(b)),nanmean(baselineblock(b,2:end)),'.','Color',colorblockid{b},'Markersize',30);
end
xticks(xblock)
xticklabels({'Early','Interm.','Late'})
ylabel('Baseline pupil size')
legend([G(1),G(2),G(3)],{'Early','Interm.','Late'},'Location','SouthWest'); legend boxoff
indexbaseline = baselineblock(1,:)-baselineblock(3,:);



figure(); hold on;
for snb= 1:length(ListSubject)
    for b = 1:3
        plot((xblock(:))+subjitt(snb)',baselineoutcome(:,snb)','Color',coutcome{b},'Markersize',10);
        plot((xblock(b))+subjitt(snb)',baselineoutcome(b,snb)','.','Color',coutcome{b},'Markersize',10);
    end
end
for b = 1:3
    G(b) = plot((xblock(b)),nanmean(baselineoutcome(b,2:end)),'.','Color',coutcome{b},'Markersize',30);
end
xticks(xblock)
xticklabels({'Fa','Hit.','Miss'})
ylabel('Baseline pupil size')
legend([G(1),G(2),G(3)],{'Early','Interm.','Late'},'Location','SouthWest'); legend boxoff
indexbaselineoutcome = baselineoutcome(1,:)-baselineoutcome(3,:);
% figure(); hold on;
% for snb= 1:length(ListSubject)
%     for b = [1,3]
%         plot((xblock(:))+subjitt(snb)',varbaselineblock(:,snb)','Color',colorblockid{b},'Markersize',10);
%         plot((xblock(b))+subjitt(snb)',varbaselineblock(b,snb)','.','Color',colorblockid{b},'Markersize',10);
%     end
% end
% for b = [1,3]
%     G(b) = plot((xblock(b)),nanmean(varbaselineblock(b,2:end)),'.','Color',colorblockid{b},'Markersize',30);
% end
% xticks(xblock)
% xticklabels({'Early','Interm.','Late'})
% ylabel('Baseline pupil size')
% legend([G(1),G(2),G(3)],{'Early','Interm.','Late'},'Location','SouthWest'); legend boxoff
%% Decoding attempt from CT
% MatEvokedCTLocked{SubjectNum}(CorrectTrials{SubjectNum},:);
for SubjectNum = ListSubject
    %     largechanges{SubjectNum} = find(ismember(ChangeSize{SubjectNum},[130,95]));
    trialpool{SubjectNum} = intersect(intersect(CorrectTrials{SubjectNum},find(Outcome{SubjectNum}==2)),find(ChangeTime{SubjectNum}>100 & ChangeTime{SubjectNum}<276));
    blockidtrialpool{SubjectNum} = BlockId{SubjectNum}(trialpool{SubjectNum});
    cttrialpool{SubjectNum} = ChangeTime{SubjectNum}(trialpool{SubjectNum});
    changesizetrialpool{SubjectNum} = ChangeSize{SubjectNum}(trialpool{SubjectNum});
    bptrialpool{SubjectNum} = ButtonTime{SubjectNum}(trialpool{SubjectNum});
    outcometrialpool{SubjectNum} = Outcome{SubjectNum}(trialpool{SubjectNum});
    for ii = 1:length(cttrialpool{SubjectNum})
        decomat{SubjectNum}(ii,:) = (matrawsubject{SubjectNum}(trialpool{SubjectNum}(ii),(cttrialpool{SubjectNum}(ii)+10)-34:cttrialpool{SubjectNum}(ii)+10)-nanmean(matrawsubject{SubjectNum}(trialpool{SubjectNum}(ii),(cttrialpool{SubjectNum}(ii):(cttrialpool{SubjectNum}(ii)+7)+10))))./nanmean(matrawsubject{SubjectNum}(trialpool{SubjectNum}(ii),(cttrialpool{SubjectNum}(ii):(cttrialpool{SubjectNum}(ii)+7)+10)));
    end
    [sortct{SubjectNum},sortctind{SubjectNum}] = sort(cttrialpool{SubjectNum});
    [sortbd{SubjectNum},sortbdind{SubjectNum}] = sort(blockidtrialpool{SubjectNum});
    [sortdiff{SubjectNum},sortdiffind{SubjectNum}] = sort(changesizetrialpool{SubjectNum});
    %     sortbddiffind{SubjectNum} =
end

% ListSubject([12,13]) = [];
% plot diff prior to change time
smalldiffmean = nanmean(cell2mat(cellfun(@(x,y) nanmean(x(find(y==60),:)),decomat(ListSubject),sortdiff(ListSubject),'UniformOutput',0)'));
stdsmalldiffmean = nanstd(cell2mat(cellfun(@(x,y) nanmean(x(find(y==60),:)),decomat(ListSubject),sortdiff(ListSubject),'UniformOutput',0)'));
largediffmean = nanmean(cell2mat(cellfun(@(x,y) nanmean(x(find(y==130),:)),decomat(ListSubject),sortdiff(ListSubject),'UniformOutput',0)'));
stdlargediffmean = nanstd(cell2mat(cellfun(@(x,y) nanmean(x(find(y==130),:)),decomat(ListSubject),sortdiff(ListSubject),'UniformOutput',0)'));

h = figure(); hold on;
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);

A = shadedErrorBar([(-34:0)*0.03],smalldiffmean,(stdsmalldiffmean./sqrt(length(ListSubject))),{'Color',colordiff{1},'LineWidth',1.5},1); ...
    C = shadedErrorBar([(-34:0)*0.03],largediffmean,(stdlargediffmean./sqrt(length(ListSubject))),{'Color',colordiff{3},'LineWidth',1.5},1);
%     xlim([xvalues(1),xvalues(length(xvalues))]); ylim([yvalues(1),yvalues(2)]);
title('Hit per difflvl');
xlabel('Time relative to Change [s]');
ylabel('Pupil size norm.');
xlim([-34*0.03,0]);
legend([A.mainLine,C.mainLine],'Hard','Easy'); legend boxoff;


% plot change lock block id
earlymean = nanmean(cell2mat(cellfun(@(x,y) nanmean(x(find(y==0),:)),decomat(ListSubject),sortbd(ListSubject),'UniformOutput',0)'));
stdsearlymean = nanstd(cell2mat(cellfun(@(x,y) nanmean(x(find(y==0),:)),decomat(ListSubject),sortbd(ListSubject),'UniformOutput',0)'));
intermmean = nanmean(cell2mat(cellfun(@(x,y) nanmean(x(find(y==1.5),:)),decomat(ListSubject),sortbd(ListSubject),'UniformOutput',0)'));
stdintermmean = nanstd(cell2mat(cellfun(@(x,y) nanmean(x(find(y==1.5),:)),decomat(ListSubject),sortbd(ListSubject),'UniformOutput',0)'));
latemean = nanmean(cell2mat(cellfun(@(x,y) nanmean(x(find(y==3),:)),decomat(ListSubject),sortbd(ListSubject),'UniformOutput',0)'));
stdlatemean = nanstd(cell2mat(cellfun(@(x,y) nanmean(x(find(y==3),:)),decomat(ListSubject),sortbd(ListSubject),'UniformOutput',0)'));

h = figure(); hold on;
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 12);

A = shadedErrorBar([(-34:0)*0.03],earlymean,(stdsearlymean./sqrt(length(ListSubject))),{'Color',colorblockid{1},'LineWidth',1.5},1); ...
    B = shadedErrorBar([(-34:0)*0.03],intermmean,(stdintermmean./sqrt(length(ListSubject))),{'Color',colorblockid{2},'LineWidth',1.5},1); ...
    C = shadedErrorBar([(-34:0)*0.03],latemean,(stdlatemean./sqrt(length(ListSubject))),{'Color',colorblockid{3},'LineWidth',1.5},1);
%     xlim([xvalues(1),xvalues(length(xvalues))]); ylim([yvalues(1),yvalues(2)]);
title('Hit per difflvl');
xlabel('Time relative to Change [s]');
ylabel('Pupil size norm.');
xlim([-34*0.03,0]);
legend([A.mainLine,B.mainLine,C.mainLine],'Early','Interm.','Late'); legend boxoff;



trainingpool = cellfun(@(x) find(ismember(x,datasample(x,round(length(x)/2),'Replace',false))),trialpool(ListSubject),'UniformOutput',0);
trainingpoollabel = cellfun(@(x,y) x(y),blockidtrialpool(ListSubject),trainingpool,'UniformOutput',0);
% cttrainingpool = cellfun(@(x,y,z) x(ismember(z,y)),cttrialpool(ListSubject),trainingpool,trialpool(ListSubject),'UniformOutput',0);
% bptrainingpool = cellfun(@(x,y,z) x(ismember(z,y)),bptrialpool(ListSubject),trainingpool,trialpool(ListSubject),'UniformOutput',0);;


testpool= cellfun(@(x,y) find(ismember([1:length(x)],y)==0),trialpool(ListSubject),trainingpool,'UniformOutput',0);
testpoollabel = cellfun(@(x,y) x(y),blockidtrialpool(ListSubject),testpool,'UniformOutput',0);
% cttestpool = cellfun(@(x,y,z) x(ismember(z,y)),cttrialpool(ListSubject),testpool,trialpool(ListSubject),'UniformOutput',0);;
% bptestpool = cellfun(@(x,y,z) x(ismember(z,y)),bptrialpool(ListSubject),testpool,trialpool(ListSubject),'UniformOutput',0);;

% for bs = 1:10
%   fa(bs) = cellfun(@(x) y(x<bs*50),trialpool(ListSubject),bptrialpool(ListSubject))
% end

vec1 = cellfun(@(x,y) nanmean(x(find(y==0),:),2),decomat(ListSubject),trainingpoollabel,'UniformOutput',0);
vec2 = cellfun(@(x,y) nanmean(x(find(y==1.5),:),2),decomat(ListSubject),trainingpoollabel,'UniformOutput',0);
vec3 = cellfun(@(x,y) nanmean(x(find(y==3),:),2),decomat(ListSubject),trainingpoollabel,'UniformOutput',0);
minnbtrial = cellfun(@(x,y,z) min([size(x,1),size(y,1),size(z,1)]),vec1,vec2,vec3,'UniformOutput',0);
w = cellfun(@(x,y,z) (x(1:z)-y(1:z)),vec1,vec3,minnbtrial,'UniformOutput',0);
b = -(dot(w{1},vec1{1}(1:17)) + dot(w{1},vec3{1}(1:17)))/2;

for i = 1:length(testpool{1})
    classif(i) = dot(mean(decomat{2}(testpool{1}(i),:),2),mean(w{1}))+b;
end
idvec1{1} = numel(intersect(find(testpoollabel{1}==0),find(classif(:)<b)))/numel(find(testpoollabel{1}==0));

