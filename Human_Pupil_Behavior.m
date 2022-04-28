function Human_Pupil_Behavior()
% load('/home/jennifer/Dropbox/Pilote12/ElectrophyStruct/20_07_03_DataGLMlatestPsychophysicsPupilHuman.mat')
% load '/home/jennifer/Dropbox/Pilote12/ElectrophyStruct/20_07_03_DataBehaviorGLMlatestPsychophysicsPupilHuman.mat'
% load('C:\Users\jlawlor3\Desktop\PsychophysicsJen\Pilote12\ElectrophyStruct\DataBehaviorGLMlatestPsychophysicsPupilHuman_190209.mat')
% load 'C:\Users\jlawlor3\Desktop\PsychophysicsJen\Pilote12\ElectrophyStruct\DataGLMlatestPsychophysicsPupilHuman_190209.mat'
Rootaddress = 'C:/Users/jlawlor3/Dropbox/Pilote12/';%' '/home/jennifer/Dropbox/Pilote12/';%%
load([Rootaddress 'ElectrophyStruct/20_07_03_DataBehaviorGLMlatestPsychophysicsPupilHuman.mat'])
load([Rootaddress 'ElectrophyStruct/20_07_03_DataGLMlatestPsychophysicsPupilHuman.mat'])



% general parameters
BinSize = data(1).BinSize;
ListSubject = [2:9 12:18];
BlockValue =[0,1.5,3];
RemoveMotor = 32;
BaselineLength = round(0.79/BinSize); % 800ms of baseline
xwindowslope = 40:110;% 50 bins is s +10 bins for baseline until 3s That's 1.2s of signal
xwindowpca = 10:110;% 50 bins is s +10 bins for baseline until 3s
% general tools
coutcome = {[0.8,0.1,0.4],[0.3,0.8,0.4],[0.1,0.4,0.8]};
colorblockid = cellfun(@(x) x/255,{[102,178,255],[0,128,255],[0,76,153]},'UniformOutput',0);
yvalues =[-100,400];
baselinewindow = 10;
xvalues = [-0.3+BinSize:BinSize:(178*BinSize)-0.3]; % Up to 5.2s
colordiff = cellfun(@(x) x/255,{[255,180,80],[255,108,0],[255,50,0]},'UniformOutput',0);%cellfun(@(x) x/255,{[255,153,51],[255,108,0],[255,70,0]},'UniformOutput',0);
diffvalue = [60,95,130]
% process pupil data
for SubjectNum = ListSubject
    clear StartTrial OnsetTrial weirdpointsline weirdpointscol MatRawPupiltmp MatRawPupilChange MatRawPupilFa RealBaselineLength BaselineSignal weirdpoints TrialLength MatRawPupil vectorstart weirdpointsline weirdpointscol allrawpupil stdthreshold
    % behavioral and stim events
    StartTrial = find(data(SubjectNum).t == 0);
    OnsetTrial = find(data(SubjectNum).s(1,:) == 1);
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
    fatrialnb{SubjectNum} = ListTrials{SubjectNum}(find(Outcome{SubjectNum}==1));
end 

blockind = unique(BlockId{2});
difflvlind = unique(DiffLvl{2});


% Plot behavior performance diff as a function of blocks
for b = 1:3
    for diff = 1:3
        for o = 1:2
            trialinddiffblock(b,diff,o,:) = cellfun(@(x,y,w,z) numel(intersect(intersect(find(x == difflvlind(diff)),intersect(find(y == blockind(b)),find(w == (o+1)))),find(z>=100 & z<267))),DiffLvl(ListSubject),BlockId(ListSubject),Outcome(ListSubject),ChangeTime(ListSubject));
        end
        perfblockdiff(b,diff,:) = squeeze(trialinddiffblock(b,diff,1,:))./(squeeze(trialinddiffblock(b,diff,1,:))+squeeze(trialinddiffblock(b,diff,2,:)));
    end
end


xoutcome = [-0.5 0 0.5];
xblock = [2 5 8];
subjitt = linspace(-0.1,0.1,length(ListSubject));
h = figure(); hold on;
for snb= 1:length(ListSubject)
    for diff = 1:3
        for b = 1:3
            %             plot((xblock(b)+xoutcome(diff))+subjitt(snb),perfblockdiff(:,diff,snb)','-','Color',colordiff{diff},'Markersize',10);
            %             plot((xblock(:)+xoutcome(diff))+subjitt(snb)',perfblockdiff(:,diff,snb)','Color',colordiff{diff},'LineWidth',0.5,'Markersize',10);
            plot((xblock(:)+xoutcome(diff))+subjitt(snb)',100*(perfblockdiff(:,diff,snb)'),'.','Color',colordiff{diff},'Markersize',10);
        end
    end
end
for b = 1:3
    for diff = 1:3
        A(diff) = plot((xblock(b)+xoutcome(diff)),nanmean(perfblockdiff(b,diff,:))*100,'.','Color',colordiff{diff},'Markersize',30);
    end
end
xticks(xblock)
xticklabels({'Early','Interm.','Late'})
ylabel('Perf [%]')
legend([A(1),A(2),A(3)],{'small','medium','large'},'Location','SouthWest'); legend boxoff
saveas(h,[Rootaddress 'PupilPaper/PerformanceBlockDiff.pdf'])
saveas(h,[Rootaddress 'PupilPaper/PerformanceBlockDiff.fig'])

% add perf as a function of time per diff per block
% hitct = cellfun(@(x,y) y(find(x==2))*0.03,Outcome(ListSubject),ChangeTime(ListSubject),'UniformOutput',0);
% missct = cellfun(@(x,y) y(find(x==2))*0.03,Outcome(ListSubject),ChangeTime(ListSubject),'UniformOutput',0);
% hitcs = cellfun(@(x,y) y(find(x==2)),Outcome(ListSubject),ChangeSize(ListSubject),'UniformOutput',0);
% misscs = cellfun(@(x,y) y(find(x==2)),Outcome(ListSubject),ChangeSize(ListSubject),'UniformOutput',0);
% for tb = 1:10
%     for diff = 1:3
%     hitnbtb(tb,diff,:) = cellfun(@(x,y) numel(intersect(find(y == difflvlind(diff)),find(x>=tb-1 & x<tb))),hitct,hitcs);
%     missnbtb(tb,diff,:) = cellfun(@(x,y) numel(intersect(find(y == difflvlind(diff)),find(x>=tb-1 & x<tb))),missct,misscs);  
%     perftb(tb,diff,:)  = ;
%     end
% end
% plot perf per as a fucntion of time per block per diff lvl npo normalise
% by trialnb
for ct = [4:8]
    for d = 1:3
        for b =1:3
            hitblockct(ct-3,d,b,:) = cellfun(@(x,y,z,w) numel(intersect(intersect(intersect(find(x==2),find(y==difflvlind(d))),find(z==BlockValue(b))),find(w>=floor((ct-1)/0.03) & w<floor(ct/0.03)))),Outcome(ListSubject),DiffLvl(ListSubject),BlockId(ListSubject),ChangeTime(ListSubject));
            missblockct(ct-3,d,b,:) = cellfun(@(x,y,z,w) numel(intersect(intersect(intersect(find(x==3),find(y==difflvlind(d))),find(z==BlockValue(b))),find(w>=floor((ct-1)/0.03) & w<floor(ct/0.03)))),Outcome(ListSubject),DiffLvl(ListSubject),BlockId(ListSubject),ChangeTime(ListSubject));
        end
    end
end
perfblockct = hitblockct./(hitblockct+missblockct);
changetrialnb = (hitblockct+missblockct);

h = figure(); hold on;
for d = 1:3
    subplot(1,3,d); hold on; title(['Diff: ' num2str(diffvalue(d))])
    for b = 1:3
        shadedErrorBar(3:7,nanmean(perfblockct(:,d,b,:)*100,4),nanstd(squeeze(perfblockct(:,d,b,:)*100)')./sqrt(length(ListSubject)),{'Color',colorblockid{b}},1);
    end
    xlim([3,7]); ylim([0 100]); xlabel('Change time [s]'); ylabel('Perf [%]')
end

saveas(h,[Rootaddress 'PupilPaper/PerformanceBlockDiffChangetime.pdf'])
saveas(h,[Rootaddress 'PupilPaper/PerformanceBlockDiffChangetime.fig'])
%test per diff lvl block effect
[smallpval,~,smallstat]=friedman([squeeze(perfblockdiff(1,1,:)),squeeze(perfblockdiff(2,1,:)),squeeze(perfblockdiff(3,1,:))])
[middlepval,~,middlestat]=friedman([squeeze(perfblockdiff(1,2,:)),squeeze(perfblockdiff(2,2,:)),squeeze(perfblockdiff(3,2,:))])
[largepval,~,largestat]=friedman([squeeze(perfblockdiff(1,3,:)),squeeze(perfblockdiff(2,3,:)),squeeze(perfblockdiff(3,3,:))])

% cumsum hit per diff lvl
figure(); hold on; 
for d = 1:3
    for b =1:3
        for snb = 1:15
            hitcdf(d,b,snb,:) = (cumsum(hitblockct(:,d,b,snb))./cumsum(changetrialnb(:,d,b,snb)))*100;
        end
        subplot(1,3,d); hold on; plot(3:7,squeeze(nanmean(hitcdf(d,b,:,:),3))','Color',colorblockid{b}); ylim([0 100])
    end
end
    



% plot reaction time per diff per block

for b = 1:3
    for diff = 1:3
        RTdiffblock(b,diff,:) = cellfun(@median,cellfun(@(x,y,w,z,bp) bp(intersect(intersect(find(x == difflvlind(diff)),intersect(find(y == blockind(b)),find(w == 2))),find(z>=100 & z<267)))-z(intersect(intersect(find(x == difflvlind(diff)),intersect(find(y == blockind(b)),find(w == 2))),find(z>=100 & z<267))),DiffLvl(ListSubject),BlockId(ListSubject),Outcome(ListSubject),ChangeTime(ListSubject),ButtonTime(ListSubject),'UniformOutput',false))*0.03;
        varRTdiffblock(b,diff,:) = cellfun(@var,cellfun(@(x,y,w,z,bp) bp(intersect(intersect(find(x == difflvlind(diff)),intersect(find(y == blockind(b)),find(w == 2))),find(z>=100 & z<267)))-z(intersect(intersect(find(x == difflvlind(diff)),intersect(find(y == blockind(b)),find(w == 2))),find(z>=100 & z<267))),DiffLvl(ListSubject),BlockId(ListSubject),Outcome(ListSubject),ChangeTime(ListSubject),ButtonTime(ListSubject),'UniformOutput',false))*0.03;
        allRTdiffblock(b,diff,:) = cellfun(@(x,y,w,z,bp) bp(intersect(intersect(find(x == difflvlind(diff)),intersect(find(y == blockind(b)),find(w == 2))),find(z>=100 & z<267)))-z(intersect(intersect(find(x == difflvlind(diff)),intersect(find(y == blockind(b)),find(w == 2))),find(z>=100 & z<267))),DiffLvl(ListSubject),BlockId(ListSubject),Outcome(ListSubject),ChangeTime(ListSubject),ButtonTime(ListSubject),'UniformOutput',false);
        allCTdiffblock(b,diff,:) = cellfun(@(x,y,w,z) z(intersect(intersect(find(x == difflvlind(diff)),intersect(find(y == blockind(b)),find(w == 2))),find(z>=100 & z<267))),DiffLvl(ListSubject),BlockId(ListSubject),Outcome(ListSubject),ChangeTime(ListSubject),'UniformOutput',false);
    end
end
% plot median reaction time per block per diff
h = figure(); hold on;
for snb= 1:length(ListSubject)
    for diff = 1:3
        for b = 1:3
            % plot((xblock(:)+xoutcome(diff))+subjitt(snb)',RTdiffblock(:,diff,snb)','Color',colordiff{diff},'LineWidth',0.5,'Markersize',10);
            plot((xblock(:)+xoutcome(diff))+subjitt(snb)',RTdiffblock(:,diff,snb)','.','Color',colordiff{diff},'Markersize',10);
        end
    end
end
for b = 1:3
    for diff = 1:3
        A(diff) = plot((xblock(b)+xoutcome(diff)),nanmean(RTdiffblock(b,diff,:)),'.','Color',colordiff{diff},'Markersize',30); ylim([0, 2])
    end
end
xticks(xblock)
xticklabels({'Early','Interm.','Late'})
ylabel('Median Reaction Time [s]')
legend([A(1),A(2),A(3)],{'small','interm.','large'},'Location','SouthWest'); legend boxoff
saveas(h,[Rootaddress 'PupilPaper/ReactionTimeBlockDiff.pdf'])
saveas(h,[Rootaddress 'PupilPaper/ReactionTimeBlockDiff.fig'])

%test per diff lvl block effect
[smallrtpval,~,smallrtstat]=friedman([squeeze(RTdiffblock(1,1,:)),squeeze(RTdiffblock(2,1,:)),squeeze(RTdiffblock(3,1,:))])
[middlertpval,~,middlertstat]=friedman([squeeze(RTdiffblock(1,2,:)),squeeze(RTdiffblock(2,2,:)),squeeze(RTdiffblock(3,2,:))])
[largertpval,~,largertstat]=friedman([squeeze(RTdiffblock(1,3,:)),squeeze(RTdiffblock(2,3,:)),squeeze(RTdiffblock(3,3,:))])

% plot var rt per block per diff
% h = figure(); hold on;
% for snb= 1:length(ListSubject)
%     for diff = 1:3
%         for b = 1:3
%             % plot((xblock(:)+xoutcome(diff))+subjitt(snb)',RTdiffblock(:,diff,snb)','Color',colordiff{diff},'LineWidth',0.5,'Markersize',10);
%             plot((xblock(:)+xoutcome(diff))+subjitt(snb)',varRTdiffblock(:,diff,snb)','.','Color',colordiff{diff},'Markersize',10);
%         end
%     end
% end
% for b = 1:3
%     for diff = 1:3
%         A(diff) = plot((xblock(b)+xoutcome(diff)),nanmean(varRTdiffblock(b,diff,:)),'.','Color',colordiff{diff},'Markersize',30); 
%     end
% end
% xticks(xblock)
% xticklabels({'Early','Interm.','Late'})
% ylabel('Variance Reaction Time')
% legend([A(1),A(2),A(3)],{'small','interm.','large'},'Location','SouthWest'); legend boxoff
% saveas(h,[Rootaddress 'PupilPaper/VarReactionTimeBlockDiff.pdf'])
%test per diff lvl block effect
% [smallrtpval,~,smallrtstat]=friedman([squeeze(RTdiffblock(1,1,:)),squeeze(RTdiffblock(2,1,:)),squeeze(RTdiffblock(3,1,:))])
% [middlertpval,~,middlertstat]=friedman([squeeze(RTdiffblock(1,2,:)),squeeze(RTdiffblock(2,2,:)),squeeze(RTdiffblock(3,2,:))])
% [largertpval,~,largertstat]=friedman([squeeze(RTdiffblock(1,3,:)),squeeze(RTdiffblock(2,3,:)),squeeze(RTdiffblock(3,3,:))])



% PULL ALL TRIALS ALL SUBJECTS TOGETHER PER CONDITION: DIFF*BLOCK
for snb= 1:length(ListSubject)
    for diff = 1:3
        allRTdiffblocknorm(:,diff,snb) = cellfun(@(x) ((x*0.03)-median(allRTdiffblock{1,3,snb}*0.03)),allRTdiffblock(:,diff,snb),'UniformOutput',false);%./std(horzcat(allRTdiffblock{:,diff,snb})*0.03)
%         medianallRTdiffblocknorm(:,diff,snb) = cellfun(@(x) ((median(x*0.03)-median(horzcat(allRTdiffblock{:,:,snb})*0.03))./max(horzcat(allRTdiffblock{:,:,snb})*0.03)),allRTdiffblock(:,diff,snb));
%         allRTdiffblocknormonlysubj(:,diff,snb) = cellfun(@(x) ((x*0.03)-median(horzcat(allRTdiffblock{:,:,snb})*0.03)),allRTdiffblock(:,diff,snb),'UniformOutput',false);%./std(horzcat(allRTdiffblock{:,diff,snb})*0.03)
    end
end

for b = 1:3
    for diff = 1:3
        normrtall{b,diff} = horzcat(allRTdiffblocknorm{b,diff,:});
        rtall{b,diff} = horzcat(allRTdiffblock{b,diff,:})*0.03;
        ctall{b,diff} = horzcat(allCTdiffblock{b,diff,:});
    end
end
figure(); hold on;
for b = 1:3
    for diff = 1:3
        A(b) = plot((xblock(diff)+xoutcome(b)),nanmean(normrtall{b,diff}),'.','Color',colorblockid{b},'Markersize',30);
        errorbar((xblock(diff)+xoutcome(b)),nanmean(normrtall{b,diff}),nanstd(normrtall{b,diff})./sqrt(length(normrtall{b,diff})),'Color',colorblockid{b},'Markersize',30);
    end
end
xticks(xblock)
xticklabels({'Small','Medium','Large'})
ylabel('Mean reaction time norm.')
legend([A(1),A(2),A(3)],{'Early','Interm.','Late'},'Location','SouthWest'); legend boxoff

% stat test for reaction time deviation from the diff effect per block overall
bla1 = [normrtall{1,3},normrtall{2,3},normrtall{3,3}];
grp = [ones(1,length(normrtall{1,3})),ones(1,length(normrtall{2,3}))*2,ones(1,length(normrtall{3,3}))*3];
[a,b,c]=anovan(bla1',{grp})
multcompare(c)

%plot rt as a function of changetime per diff
figure(); hold on;
for diff = 1:3
    for snb = 1:15
        clear poolrt ctsort indsort
        poolrt = horzcat(allRTdiffblock{:,diff,snb})*0.03;
        [ctsort,indsort]=sort(horzcat(allCTdiffblock{:,diff,snb})*0.03);
        %         plot(ctsort,poolrt(indsort),'.','Color',colordiff{diff},'Markersize',5);
        for bt = 4:8
            if ~isempty(find(ctsort>=bt-1 & ctsort<bt))
                plot(bt-0.05,median(poolrt(indsort(find(ctsort>=bt-1 & ctsort<bt)))),'.','Color',colordiff{diff},'Markersize',10);
                medianrttime(diff,snb,bt-3) = median(poolrt(indsort(find(ctsort>=bt-1 & ctsort<bt))));
            else
                medianrttime(diff,snb,bt-3) = nan;
            end
        end
%         plot([4:8]-0.05,squeeze(nanmean(medianrttime(diff,:,:),2))','Color',colordiff{diff});        
    end
            plot([4:8]-0.05,squeeze(nanmean(medianrttime(diff,:,:),2))','.','Color',colordiff{diff},'Markersize',30);
end
% test effect of time of rt per diff
[rttimepvalsmall] = friedman(squeeze(medianrttime(1,:,1:4)));
[rttimepvalmedium] = friedman(squeeze(medianrttime(2,:,:)));
[rttimepvallarge] = friedman(squeeze(medianrttime(3,:,:)));
% RT for trial 3s to 8 s as a function of time per block
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

% compute instantaneous fa and hit and miss per block
timingfa = cellfun(@(x,y) y(find(x==1)),Outcome(ListSubject),ButtonTime(ListSubject),'UniformOutput',0);
timingfablock = cellfun(@(x,y) (y(find(x==1))),Outcome(ListSubject),BlockId(ListSubject),'UniformOutput',0);
timingfaseq = cellfun(@(x,y) (y(find(x==1))),Outcome(ListSubject),ListTrials(ListSubject),'UniformOutput',0);
blockstart = [0,120,240];
figure(15); hold on;
figure(16); hold on; 
% for snb =1:15
% for b = 1:3
%     clear a bli
%  [a,bli] = hist(timingfa{snb}(timingfablock{snb}==BlockValue(b)),5)   
%  figure(15); hold on; subplot(2,8,snb); hold on; plot(bli*0.03,a/max(a),'Color',colorblockid{b});
%  [~,blockstarind] = sort(abs(blockstart-min(timingfaseq{snb}(timingfablock{snb}==BlockValue(b)))));
%  figure(16); hold on; subplot(2,8,snb); hold on; plot(1:size(timingfaseq{snb}(timingfablock{snb}==BlockValue(b)),2),timingfaseq{snb}(timingfablock{snb}==BlockValue(b))-blockstart(blockstarind(1)),'Color',colorblockid{b});
%  figure(17); hold on; subplot(2,8,snb); hold on; plot(timingfaseq{snb}(timingfablock{snb}==BlockValue(b))-blockstart(blockstarind(1)),timingfa{snb}(timingfablock{snb}==BlockValue(b))*0.03,'.','Color',colorblockid{b},'MarkerSize',20);
%  trialnbblock3{snb,b} = timingfaseq{snb}(((timingfa{snb}(timingfablock{snb}==BlockValue(b)))*0.03)<=3);
% end
% end

binboundaries = [0.5:0.5:7];
sizewindow = 1
for b = 1:3
    for tb = 1:length(binboundaries)
     
        fabin(tb,b,:) = cellfun(@(x,y,z) numel(intersect(intersect(find(x == 1),find((y*0.03)>=binboundaries(tb)-sizewindow/2 & (y*0.03)<(binboundaries(tb)+sizewindow/2))),find(z == BlockValue(b)))),Outcome(ListSubject),ButtonTime(ListSubject),BlockId(ListSubject));
        hitbin(tb,b,:) = cellfun(@(x,y,z) numel(intersect(intersect(find(x == 2),find((y*0.03)>=binboundaries(tb)-sizewindow/2 & (y*0.03)<(binboundaries(tb)+sizewindow/2))),find(z == BlockValue(b)))),Outcome(ListSubject),ChangeTime(ListSubject),BlockId(ListSubject));
        missbin(tb,b,:) = cellfun(@(x,y,z) numel(intersect(intersect(find(x == 3),find((y*0.03)>=binboundaries(tb)-sizewindow/2 & (y*0.03)<(binboundaries(tb)+sizewindow/2))),find(z == BlockValue(b)))),Outcome(ListSubject),ChangeTime(ListSubject),BlockId(ListSubject));
        crbin(tb,b,:) = cellfun(@(x,y,z) numel(intersect(unique([find((x*0.03)>=binboundaries(tb)+sizewindow/2),find((y*0.03)>=binboundaries(tb)+sizewindow/2)]),find(z == BlockValue(b)))),ButtonTime(ListSubject),ChangeTime(ListSubject),BlockId(ListSubject));
        farate(tb,b,:) = fabin(tb,b,:)./(fabin(tb,b,:)+crbin(tb,b,:));
        hitrate(tb,b,:) = hitbin(tb,b,:)./(hitbin(tb,b,:)+missbin(tb,b,:));
        
        farate(tb,b,(find(farate(tb,b,:)==0))) = 1/(2*((fabin(tb,b,find(farate(tb,b,:)==0)))+crbin(tb,b,find(farate(tb,b,:)==0))));
        farate(tb,b,(find(farate(tb,b,:)==1))) = 1-1/(2*((fabin(tb,b,find(farate(tb,b,:)==1)))+crbin(tb,b,find(farate(tb,b,:)==1))));
        hitrate(tb,b,(find(hitrate(tb,b,:)==0))) = 1/(2*((hitbin(tb,b,find(hitrate(tb,b,:)==0)))+missbin(tb,b,find(hitrate(tb,b,:)==0))));
        hitrate(tb,b,(find(hitrate(tb,b,:)==1))) = 1-1/(2*((hitbin(tb,b,find(hitrate(tb,b,:)==1)))+missbin(tb,b,find(hitrate(tb,b,:)==1))));
           
        % for quarter block
        faindtrialnb(tb,b,:) = cellfun(@(x,y,z,w) (w(intersect(intersect(find(x == 1),find((y*0.03)>=binboundaries(tb)-sizewindow/2 & (y*0.03)<(binboundaries(tb)+sizewindow/2))),find(z == BlockValue(b))))-find(z == BlockValue(b),1))+1,Outcome(ListSubject),ButtonTime(ListSubject),BlockId(ListSubject),ListTrials(ListSubject),'UniformOutput',0);
        crindtrialnb(tb,b,:) = cellfun(@(x,y,z,w) (w(intersect([find((x*0.03)>=binboundaries(tb)-sizewindow/2),find((y*0.03)>=binboundaries(tb)-sizewindow/2)],find(z == BlockValue(b))))-find(z == BlockValue(b),1))+1,ButtonTime(ListSubject),ChangeTime(ListSubject),BlockId(ListSubject),ListTrials(ListSubject),'UniformOutput',0);
        hitindtrialnb(tb,b,:) = cellfun(@(x,y,z,w) (w(intersect(intersect(find(x == 2),find((y*0.03)>=binboundaries(tb)-sizewindow/2 & (y*0.03)<(binboundaries(tb)+sizewindow/2))),find(z == BlockValue(b))))-find(z == BlockValue(b),1))+1,Outcome(ListSubject),ButtonTime(ListSubject),BlockId(ListSubject),ListTrials(ListSubject),'UniformOutput',0);
        missindtrialnb(tb,b,:) = cellfun(@(x,y,z,w) (w(intersect(intersect(find(x == 3),find((y*0.03)>=binboundaries(tb)-sizewindow/2 & (y*0.03)<(binboundaries(tb)+sizewindow/2))),find(z == BlockValue(b))))-find(z == BlockValue(b),1))+1,Outcome(ListSubject),ButtonTime(ListSubject),BlockId(ListSubject),ListTrials(ListSubject),'UniformOutput',0);
        for q = 1:4
            fanumel(tb,b,q,:)=cell2mat(squeeze(cellfun(@(x) numel(find(x>((q-1)*10)+1 & x<=q*10)),faindtrialnb(tb,b,:),'UniformOutput',false)));
            crnumel(tb,b,q,:)=cell2mat(squeeze(cellfun(@(x) numel(find(x>((q-1)*10)+1 & x<=q*10)),crindtrialnb(tb,b,:),'UniformOutput',false)));
            hitnumel(tb,b,q,:)=cell2mat(squeeze(cellfun(@(x) numel(find(x>((q-1)*10)+1 & x<=q*10)),hitindtrialnb(tb,b,:),'UniformOutput',false)));
            missnumel(tb,b,q,:)=cell2mat(squeeze(cellfun(@(x) numel(find(x>((q-1)*10)+1 & x<=q*10)),missindtrialnb(tb,b,:),'UniformOutput',false)));
            % compute the rates
            faratequarter(tb,b,q,:) = (fanumel(tb,b,q,:)./(fanumel(tb,b,q,:)+crnumel(tb,b,q,:)));
            hitratequarter(tb,b,q,:) = (hitnumel(tb,b,q,:)./(hitnumel(tb,b,q,:)+missnumel(tb,b,q,:)));
            % correct for either 0 fa or hit and 100% hit or fa
            faratequarter(tb,b,q,find(faratequarter(tb,b,q,:)==0)) =  1/(2*crnumel(tb,b,q,find(faratequarter(tb,b,q,:)==0)));
            faratequarter(tb,b,q,find(faratequarter(tb,b,q,:)==1)) =  1 - 1/(2*fanumel(tb,b,q,find(faratequarter(tb,b,q,:)==1)));
            hitratequarter(tb,b,q,find(hitratequarter(tb,b,q,:)==0)) =  1/(2*missnumel(tb,b,q,find(hitratequarter(tb,b,q,:)==0)));
            hitratequarter(tb,b,q,find(hitratequarter(tb,b,q,:)==1)) =  1 - 1/(2*hitnumel(tb,b,q,find(hitratequarter(tb,b,q,:)==1)));
            dprimeratequarter(tb,b,q,:) = dprime(squeeze(hitratequarter(tb,b,q,:))',squeeze(faratequarter(tb,b,q,:))');
            criterionquarter(tb,b,q,:) = -0.5*(norminv(squeeze(hitratequarter(tb,b,q,:))')+norminv(squeeze(faratequarter(tb,b,q,:))'))';
        end
        %         pnoise(tb,b,:) = cellfun(@(x,z) numel(intersect(find((x*0.03)>=binboundaries(tb)+sizewindow/2),find(z == BlockValue(b)))),ChangeTime(ListSubject),BlockId(ListSubject));
%         psignal(tb,b,:) = cellfun(@(x,z) numel(intersect(find((x*0.03)>=binboundaries(tb) &(x*0.03)<binboundaries(tb)+sizewindow/2),find(z == BlockValue(b)))),ChangeTime(ListSubject),BlockId(ListSubject));
        %     likelihood(tb,b,:)=pnoise(tb,b,:)./psignal(tb,b,:);
    
    % do it for difflvl
    for diff = 1:3
        
       
%         fabindiff(diff,tb,b,:) = cellfun(@(x,y,z,r) numel(intersect(intersect(intersect(find(x == 1),find((y*0.03)>=binboundaries(tb)-sizewindow/2 & (y*0.03)<(binboundaries(tb)+sizewindow/2))),find(z == BlockValue(b))),find(r==difflvlind(diff)))),Outcome(ListSubject),ButtonTime(ListSubject),BlockId(ListSubject),DiffLvl(ListSubject));
        hitbindiff(diff,tb,b,:) = cellfun(@(x,y,z,r) numel(intersect(intersect(intersect(find(x == 2),find((y*0.03)>=binboundaries(tb)-sizewindow/2 & (y*0.03)<(binboundaries(tb)+sizewindow/2))),find(z == BlockValue(b))),find(r==difflvlind(diff)))),Outcome(ListSubject),ChangeTime(ListSubject),BlockId(ListSubject),DiffLvl(ListSubject));
        missbindiff(diff,tb,b,:) = cellfun(@(x,y,z,r) numel(intersect(intersect(intersect(find(x == 3),find((y*0.03)>=binboundaries(tb)-sizewindow/2 & (y*0.03)<(binboundaries(tb)+sizewindow/2))),find(z == BlockValue(b))),find(r==difflvlind(diff)))),Outcome(ListSubject),ChangeTime(ListSubject),BlockId(ListSubject),DiffLvl(ListSubject));
%         crbindiff(diff,tb,b,:) = cellfun(@(x,y,z,r) numel(intersect(intersect(unique([find((x*0.03)>=binboundaries(tb)+sizewindow/2),find((y*0.03)>=binboundaries(tb)+sizewindow/2)]),find(z == BlockValue(b))),find(r==difflvlind(diff)))),ButtonTime(ListSubject),ChangeTime(ListSubject),BlockId(ListSubject),DiffLvl(ListSubject));
%         faratediff(diff,tb,b,:) = fabindiff(diff,tb,b,:)./(fabindiff(diff,tb,b,:)+crbindiff(diff,tb,b,:));
        hitratediff(diff,tb,b,:) = hitbindiff(diff,tb,b,:)./(hitbindiff(diff,tb,b,:)+missbindiff(diff,tb,b,:));
        
%         faratediff(diff,tb,b,(find(faratediff(diff,tb,b,:)==0))) = 1/(2*((fabindiff(diff,tb,b,find(faratediff(diff,tb,b,:)==0)))+crbindiff(diff,tb,b,find(faratediff(diff,tb,b,:)==0))));
%         faratediff(diff,tb,b,(find(faratediff(diff,tb,b,:)==1))) = 1-1/(2*((fabindiff(diff,tb,b,find(faratediff(diff,tb,b,:)==1)))+crbindiff(diff,tb,b,find(faratediff(diff,tb,b,:)==1))));
        hitratediff(diff,tb,b,(find(hitratediff(diff,tb,b,:)==0))) = 1/(2*((hitbindiff(diff,tb,b,find(hitratediff(diff,tb,b,:)==0)))+missbindiff(diff,tb,b,find(hitratediff(diff,tb,b,:)==0))));
        hitratediff(diff,tb,b,(find(hitratediff(diff,tb,b,:)==1))) = 1-1/(2*((hitbindiff(diff,tb,b,find(hitratediff(diff,tb,b,:)==1)))+missbindiff(diff,tb,b,find(hitratediff(diff,tb,b,:)==1))));
    end
   end 
end



for snb = 1:15
    for b = 1:3
        dprimerate(snb,b,:) = dprime(squeeze(hitrate(:,b,snb))',squeeze(farate(:,b,snb))');
        iocriterionrate(snb,b,:) = 0.5*dprimerate(snb,b,:);
        criterion(snb,b,:) = -0.5*(norminv(squeeze(hitrate(:,b,snb))')+norminv(squeeze(farate(:,b,snb))'))';
        lnlikelihood(snb,b,:) = squeeze(dprimerate(snb,b,:))'.*squeeze(criterion(snb,b,:))';
        for diff = 1:3
            dprimeratediff(diff,snb,b,:) = dprime(squeeze(hitratediff(diff,:,b,snb))',squeeze(farate(:,b,snb)));
            iocriterionratediff(diff,snb,b,:) = 0.5*dprimeratediff(diff,snb,b,:);
            criteriondiff(diff,snb,b,:) = -0.5*(norminv(squeeze(hitratediff(diff,:,b,snb))')+norminv(squeeze(farate(:,b,snb))))';
            lnlikelihooddiff(diff,snb,b,:) = squeeze(dprimeratediff(diff,snb,b,:))'.*squeeze(criteriondiff(diff,snb,b,:))';
        end
    end
end

% plot farate as function of time
h = figure(); hold on;
for b = 1:3
    shadedErrorBar(binboundaries(1:end),squeeze(nanmean(hitrate(:,b,:),3))'*100,nanstd(squeeze(hitrate(:,b,:))')*100/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b},'LineWidth',1.5},1); ylim([0,100]); ylabel('Instant. Hit [%]'); xlabel('Time from onset [s]');
end
saveas(h,[Rootaddress 'PupilPaper/InstantaneaousHitRate.pdf'])
saveas(h,[Rootaddress 'PupilPaper/InstantaneaousHitRate.fig'])


% plot farate as function of time
h = figure(); hold on;
for b = 1:3
    shadedErrorBar(binboundaries(1:end),squeeze(nanmean(farate(:,b,:),3))'*100,nanstd(squeeze(farate(:,b,:))')*100/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b},'LineWidth',1.5},1); ylim([0,50]); ylabel('Instant. Fa [%]'); xlabel('Time from onset [s]'); xlim([binboundaries(1),binboundaries(end)])
end
saveas(h,[Rootaddress 'PupilPaper/InstantaneaousFARate.pdf'])
saveas(h,[Rootaddress 'PupilPaper/InstantaneaousFARate.fig'])

% plot dprime and criterion
% 
h = figure(); hold on;
for b = 1:3
    subplot(1,2,1); hold on; shadedErrorBar(binboundaries(1:end),squeeze(nanmean(dprimerate(:,b,:),1))',nanstd(squeeze(dprimerate(:,b,:)))/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b},'LineWidth',1.5},1); ylim([0,3]); ylabel('Instant. d'''); xlabel('Time from onset [s]'); xlim([3,7])
    subplot(1,2,2); hold on; shadedErrorBar(binboundaries(1:end),squeeze(nanmean(criterion(:,b,:),1))',nanstd(squeeze(criterion(:,b,:)))/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b},'LineWidth',1.5},1); ylim([-0.5,1]); ylabel('Instant. critterion'); xlabel('Time from onset [s]'); xlim([3,7])
end
saveas(h,[Rootaddress 'PupilPaper/DprimeCriterionBlock.pdf'])

% print d' and criterion figure independently for the paper

h = figure(); hold on;
for b = 1:3
    hold on; shadedErrorBar(binboundaries(1:end),squeeze(nanmean(dprimerate(:,b,:),1))',nanstd(squeeze(dprimerate(:,b,:)))/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b},'LineWidth',1.5},1); ylim([1,3]); ylabel('Instantaneous d'''); xlabel('Time from onset [s]'); xlim([3,7])
end
saveas(h,[Rootaddress 'PupilPaper/DprimeBlock.fig'])


h = figure(); hold on;
for b = 1:3
    hold on; shadedErrorBar(binboundaries(1:end),squeeze(nanmean(criterion(:,b,:),1))',nanstd(squeeze(criterion(:,b,:)))/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b},'LineWidth',1.5},1); ylim([-0.5,1]); ylabel('Instantaneous critterion'); xlabel('Time from onset [s]'); xlim([3,7])
end
saveas(h,[Rootaddress 'PupilPaper/CriterionBlock.fig'])

% plot dprime and criterion per block per diff
h = figure(); hold on;
for snb= 1:length(ListSubject)
    for diff = 1:3
            % plot((xblock(:)+xoutcome(diff))+subjitt(snb)',RTdiffblock(:,diff,snb)','Color',colordiff{diff},'LineWidth',0.5,'Markersize',10);
            plot((xblock(:)+xoutcome(diff))+subjitt(snb)',squeeze(nanmean(dprimeratediff(diff,snb,:,6:14),4)),'.','Color',colordiff{diff},'Markersize',10);
    end
end
for b = 1:3
    for diff = 1:3
        A(diff) = plot((xblock(b)+xoutcome(diff)),nanmean(squeeze(nanmean(dprimeratediff(diff,:,b,6:14),4)),2),'.','Color',colordiff{diff},'Markersize',30); ylim([0, 4])
    end
end
xticks(xblock)
xticklabels({'Early','Interm.','Late'})
ylabel('d''')
legend([A(1),A(2),A(3)],{'small','interm.','large'},'Location','SouthWest'); legend boxoff
saveas(h,[Rootaddress 'PupilPaper/DprimeBlockDiff.pdf'])
saveas(h,[Rootaddress 'PupilPaper/DprimeBlockDiff.fig'])

% same but for criterion (per block per diff)
h = figure(); hold on;
for snb= 1:length(ListSubject)
    for diff = 1:3
            % plot((xblock(:)+xoutcome(diff))+subjitt(snb)',RTdiffblock(:,diff,snb)','Color',colordiff{diff},'LineWidth',0.5,'Markersize',10);
            plot((xblock(:)+xoutcome(diff))+subjitt(snb)',squeeze(nanmean(criteriondiff(diff,snb,:,6:14),4)),'.','Color',colordiff{diff},'Markersize',10);
    end
end
for b = 1:3
    for diff = 1:3
        A(diff) = plot((xblock(b)+xoutcome(diff)),nanmean(squeeze(nanmean(criteriondiff(diff,:,b,6:14),4)),2),'.','Color',colordiff{diff},'Markersize',30); ylim([0, 1.5])
    end
end
xticks(xblock)
xticklabels({'Early','Interm.','Late'})
ylabel('Criterion')
legend([A(1),A(2),A(3)],{'small','interm.','large'},'Location','SouthWest'); legend boxoff
saveas(h,[Rootaddress 'PupilPaper/CriterionBlockDiff.pdf'])
saveas(h,[Rootaddress 'PupilPaper/CriterionBlockDiff.fig'])


PupilBehaviorTestStruct.MedianReactionTime = RTdiffblock;
PupilBehaviorTestStruct.HitRate = perfblockdiff;
PupilBehaviorTestStruct.hitrate = hitrate;
PupilBehaviorTestStruct.farate = farate;
PupilBehaviorTestStruct.dprimerate = dprimerate;
PupilBehaviorTestStruct.criterion = criterion;
PupilBehaviorTestStruct.binboundaries = binboundaries;

% plot criterion, hit rate, per block per quarter block
h = figure(); hold on;
for q = 1:4
    for b = 1:3
        subplot(1,4,q);hold on; shadedErrorBar(binboundaries(1:end),squeeze(nanmean(faratequarter(:,b,q,:),4)*100)',nanstd(squeeze(faratequarter(:,b,q,:)*100)')/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b}},1); ylim([0,50]); xlim ([binboundaries(1), binboundaries(end)]); ylabel('Fa rate [%]'); xlabel('Time from onset [s]'); title([num2str(q) '/4'])
%         subplot(2,4,q+4);hold on;  shadedErrorBar(binboundaries(1:end),squeeze(nanmean(criterionquarter(:,b,q,:),4))',nanstd(squeeze(criterionquarter(:,b,q,:))')/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b}},1); ylim([-0.7,0.5]); ylabel('Criterion'); xlabel('Time from onset [s]');
    end
end
saveas(h,[Rootaddress 'PupilPaper/FaRateBlockPerquarter.pdf'])

h = figure(); hold on;
for q = 1:4
    for b = 1:3
        subplot(1,4,q);hold on; shadedErrorBar(binboundaries(1:end),squeeze(nanmean(dprimeratequarter(:,b,q,:),4))',nanstd(squeeze(dprimeratequarter(:,b,q,:))')/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b}},1); ylim([0,2]); xlim ([binboundaries(1), binboundaries(end)]); ylabel('d'''); xlabel('Time from onset [s]'); title([num2str(q) '/4'])
%         subplot(2,4,q+4);hold on;  shadedErrorBar(binboundaries(1:end),squeeze(nanmean(criterionquarter(:,b,q,:),4))',nanstd(squeeze(criterionquarter(:,b,q,:))')/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b}},1); ylim([-0.7,0.5]); ylabel('Criterion'); xlabel('Time from onset [s]');
    end
end
saveas(h,[Rootaddress 'PupilPaper/dprimeBlockPerquarter.pdf'])
h = figure(); hold on;
for q = 1:4
    for b = 1:3
        subplot(1,4,q);hold on; shadedErrorBar(binboundaries(1:end),squeeze(nanmean(criterionquarter(:,b,q,:),4))',nanstd(squeeze(criterionquarter(:,b,q,:))')/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b}},1); ylim([0,1]); xlim ([binboundaries(1), binboundaries(end)]); ylabel('Criterion'); xlabel('Time from onset [s]'); title([num2str(q) '/4'])
%         subplot(2,4,q+4);hold on;  shadedErrorBar(binboundaries(1:end),squeeze(nanmean(criterionquarter(:,b,q,:),4))',nanstd(squeeze(criterionquarter(:,b,q,:))')/(sqrt(length(ListSubject)))*2,{'Color',colorblockid{b}},1); ylim([-0.7,0.5]); ylabel('Criterion'); xlabel('Time from onset [s]');
    end
end
saveas(h,[Rootaddress 'PupilPaper/criterionBlockPerquarter.pdf'])

% fa test vector
fatestvector = reshape(faratequarter,1,size(faratequarter,1)*size(faratequarter,2)*size(faratequarter,3)*size(faratequarter,4));

% hitrateblockquarter = cellfun(@(x) numel(find(x == 2))/(numel(find(x == 3))+ (numel(find(x == 2)))),outcomeblock) 
% farateblockquarter = cellfun(@(x) numel(find(x == 2))/(numel(find(x == 3))+ (numel(find(x == 2)))),outcomeblock) 

PupilBehaviorTestStruct.fatrialnb = fatrialnb(ListSubject);
PupilBehaviorTestStruct.timingfa = timingfa;
PupilBehaviorTestStruct.timingfablock = timingfablock;
PupilBehaviorTestStruct.changetime = ChangeTime(ListSubject);
PupilBehaviorTestStruct.blocid = BlockId(ListSubject);


save([Rootaddress '/PupilPaper/PupilBehaviorTestStruct_031121.mat'],'PupilBehaviorTestStruct')


% plot cdf for false alarm
h = figure(); hold on;
for b = 1:3
    
fablocktiming{b} = cellfun(@(x,y) x(intersect(find(y==BlockValue(b)), find(x<266))),timingfa,timingfablock,'UniformOutput',0);
for snb = 1:15
for tt = 1:10
    fatimebinind(b,snb,tt) = cellfun(@(x) numel(find(x>=tt-1 & x<tt)),{[fablocktiming{b}{snb}]*0.03});
end
cumsumfasubject(b,snb,:) = cumsum(squeeze(fatimebinind(b,snb,:)))./max(cumsum(squeeze(fatimebinind(b,snb,:))));
AUCfacdf(b,snb) = trapz(squeeze(cumsumfasubject(b,snb,:)));
subplot(1,3,[1,2]); hold on; plot((0:tt),[0;squeeze(cumsumfasubject(b,snb,:))],'LineWidth',0.1,'Color',colorblockid{b});
% subplot(1,4,1); hold on; 
% plot(b+subjitt(snb),sum(squeeze(fatimebinind(b,snb,:))'),'.','Color',colorblockid{b},'Markersize',10);
subplot(1,3,3); hold on; 
plot(b+subjitt(snb),AUCfacdf(b,snb),'.','Color',colorblockid{b},'Markersize',10);
end
% subplot(1,3,[1,2]); hold on; plot((0:tt),[0;squeeze(mean(cumsumfasubject(b,:,:),2))],'LineWidth',4,'Color',colorblockid{b});
xlabel('Time within trial [s]'); ylabel('CDF false alarm'); xlim([0,tt-1])

subplot(1,3,3); hold on;
N(b) = plot(b,mean(AUCfacdf(b,:)),'.','Color',colorblockid{b},'Markersize',30); ylim([3,8]);


% subplot(1,4,1); hold on; 
% M(b) = plot(b,mean(sum(squeeze(fatimebinind(b,:,:))')),'.','Color',colorblockid{b},'Markersize',30);

end

% xticks(1:3)
% xticklabels({'Early','Interm.','Late'})
% ylabel('Average nb of false alarm [0-8s]')
% legend([M(1),M(2),M(3)],{'Early','Interm.','Late'},'Location','SouthWest'); legend boxoff
% only for figure purpose

subplot(1,3,3); hold on;
xticks(1:3)
xticklabels({'Early','Interm.','Late'})
ylabel('Average AUC false alarm CDF [0-8s]')
legend([N(1),N(2),N(3)],{'Early','Interm.','Late'},'Location','SouthWest'); legend boxoff
[pfablock,fatbl,fatimestat] = friedman(AUCfacdf');
str = ['p = ', num2str(pfablock)];
text(1.75,7.75,str)
print([Rootaddress 'PupilPaper/FABlockCDF.pdf'],'-dpdf','-bestfit')
saveas(h,[Rootaddress 'PupilPaper/FABlockCDF.fig'])
% 
% [fatimestatcorrected] = multcompare(fatimestat);
% [pfablock,stats,fatimestat] = anova1(AUCfacdf');
% multcompare(fatimestat)
% test repeated measure model sphericity
% taucfa =  table([1:15]',squeeze(AUCfacdf(1,:)'),squeeze(AUCfacdf(2,:)'),squeeze(AUCfacdf(3,:)'),'VariableNames',{'Subject','b1','b2','b3'});
% WithinStructure = table(([1 2 3]'),'VariableNames',{'bloc'});
% rm = fitrm(taucfa,'b1-b3~Subject','WithinDesign',WithinStructure);
% ranovatable = ranova(rm,'WithinModel','bloc');
% mauchly(rm)

% plot increase in hit rate as a function of increase in false alarm rate
% figure(); hold on;
% for i = [1,3]
% plot(squeeze(nanmean(farate(:,1,:)*100))-squeeze(nanmean(farate(:,3,:)*100)),(squeeze(perfblockdiff(1,i,:)*100)-squeeze(perfblockdiff(3,i,:)*100)),'.','MarkerSize',25,'Color',colordiff{i})
% end
% xlabel('Early-Late fa rate [%] (3-8s)'); ylabel('Early-Late hit rate [%] (3-8s)')
% 
% linfit = fitlm(squeeze(mean(farate(:,1,:)*100))-squeeze(mean(farate(:,3,:)*100)),(squeeze(perfblockdiff(1,1,:)*100)-squeeze(perfblockdiff(3,1,:)*100)),'linear');
% valuesfa =floor(min(ratiofa)):ceil(max(ratiofa));
% subplot(1,2,2); hold on
% h4 = plot(valuesfa,(valuesfa.*linfit.Coefficients.Estimate(2))+linfit.Coefficients.Estimate(1),'k','LineWidth',1.5); ylim([-20,50]); 
% for snb= 1:length(ListSubject)
%     h3(snb) = plot(ratiofa(snb),diffpupilsize(snb),'.','MarkerSize',20,'Color',[0.7,0.7,0.7]); %text((mean(PupilBehaviorTestStruct.farate(:,1,snb),1)./mean(PupilBehaviorTestStruct.farate(:,3,snb),1))+0.01,(areablock(ListSubject(snb),1)-areablock(ListSubject(snb),3))-1,num2str(snb));
% end


% Plot chronometric plots * diff
medianrtcond = cellfun(@median,rtall);
semrtcond = cellfun(@(x) std(x)./sqrt(length(x)),rtall);
meanmedianrt = squeeze(mean(RTdiffblock,3));
semmedianrt =  (squeeze(std(permute(RTdiffblock(:,:,:),[3,1,2]))./sqrt(15)));
meanrtnorm = mean(cellfun(@median, allRTdiffblocknorm),3);
stdrtnorm = squeeze(std(permute(cellfun(@median, allRTdiffblocknorm),[3,1,2]))./(sqrt(15)));
for b = 1:3
    for snb = 1:15
        for diff = 1:3
            difffromlargestdiffwithinblock(b,diff,snb) = mean(allRTdiffblock{b,diff,snb}-mean(allRTdiffblock{1,3,snb}))*0.03;
            perfdifffromlargestdiffwithinblock(b,diff,snb) = (perfblockdiff(b,diff,snb)-perfblockdiff(1,3,snb));
        end
    end
end

h = figure(); hold on;
for diff = [1,2,3]
     plot(squeeze(mean(difffromlargestdiffwithinblock(:,diff,:),3)),squeeze(mean(perfdifffromlargestdiffwithinblock(:,diff,:),3))*100,'LineWidth',2,'Color', colordiff{diff});
    %plot(squeeze(mean(RTdiffblock(b,:,:),3)),squeeze(mean(perfblockdiff(b,:,:),3))*100,'LineWidth',2,'Color', colorblockid{b});
end

for b = 1:3
    plot(squeeze(mean(difffromlargestdiffwithinblock(b,:,:),3)),squeeze(mean(perfdifffromlargestdiffwithinblock(b,:,:),3))*100,'.','MarkerSize', 30,'Color', colorblockid{b});
%     plot(squeeze(median(RTdiffblock(3,diff,:),3)),squeeze(mean(perfblockdiff(3,diff,:),3))*100,'.','MarkerSize', 20,'Color', colordiff{diff});
end
xlabel({'Mean RT -  RT_e_a_r_l_y_/_l_a_r_g_e [s]'}); ylabel({'Mean Perf -  Perf_e_a_r_l_y_/_l_a_r_g_e [%]'})

semperfsubj = (squeeze(std(permute(perfdifffromlargestdiffwithinblock(:,:,:)*100,[3,1,2])))./sqrt(15));
meanperfsubj = squeeze(mean(permute(perfdifffromlargestdiffwithinblock(:,:,:)*100,[3,1,2])));
semrtsubj = (squeeze(std(permute(difffromlargestdiffwithinblock(:,:,:),[3,1,2])))./sqrt(15));
meanrtsubj = squeeze(mean(difffromlargestdiffwithinblock(:,:,:),3));
for b = 1:3
    for diff = 1:3
        plot([squeeze(meanrtsubj(b,diff)),squeeze(meanrtsubj(b,diff))],[meanperfsubj(b,diff)-semperfsubj(b,diff),meanperfsubj(b,diff)+semperfsubj(b,diff)],'Color', colorblockid{b},'LineWidth',1.5);
        plot([meanrtsubj(b,diff)-semrtsubj(b,diff),meanrtsubj(b,diff)+semrtsubj(b,diff)],[meanperfsubj(b,diff),meanperfsubj(b,diff)],'Color', colordiff{diff},'LineWidth',1.5);
        %     plot(squeeze(median(RTdiffblock(3,diff,:),3)),squeeze(mean(perfblockdiff(3,diff,:),3))*100,'.','MarkerSize', 20,'Color', colordiff{diff});
    end
end

print([Rootaddress 'PupilPaper/PerfRTSummaryplot.pdf'],'-dpdf','-bestfit')
saveas(h,[Rootaddress 'PupilPaper/PerfRTSummaryplot.fig'])




% plot perf per diffuclty as a function of change time and rt per diff as a
% function of change time for Supplementary figure 1
binboundariessupp = [0:2:10];
sizewindow = 2;

for tb = 1:length(binboundariessupp)
    for diff = 1:3
        hitbindiffsupp(tb,diff,:) = cellfun(@(x,y,z) numel(intersect(intersect(find(x == 2),find((y*0.03)>=binboundariessupp(tb)-sizewindow/2 & (y*0.03)<(binboundariessupp(tb)+sizewindow/2))),find(z == diffvalue(diff)))),Outcome(ListSubject),ChangeTime(ListSubject),DiffLvl(ListSubject)); % b is replace by diff here to compute the hit rate
        missbindiffsupp(tb,diff,:) = cellfun(@(x,y,z) numel(intersect(intersect(find(x == 3),find((y*0.03)>=binboundariessupp(tb)-sizewindow/2 & (y*0.03)<(binboundariessupp(tb)+sizewindow/2))),find(z == diffvalue(diff)))),Outcome(ListSubject),ChangeTime(ListSubject),DiffLvl(ListSubject)); % b is replace by diff here to compute the hit rate
        hitratediffsupp(tb,diff,:) = hitbindiffsupp(tb,diff,:)./(hitbindiffsupp(tb,diff,:)+missbindiffsupp(tb,diff,:));
        rtbindiffsupp(tb,diff,:) =  cellfun(@nanmedian,cellfun(@(x,w,y,bp) bp(intersect(intersect(find(x == difflvlind(diff)),find(w == 2)),find((y*0.03)>=binboundariessupp(tb)-sizewindow/2 & (y*0.03)<(binboundariessupp(tb)+sizewindow/2))))-y(intersect(intersect(find(x == difflvlind(diff)),find(w == 2)),find((y*0.03)>=binboundariessupp(tb)-sizewindow/2 & (y*0.03)<(binboundariessupp(tb)+sizewindow/2)))),DiffLvl(ListSubject),Outcome(ListSubject),ChangeTime(ListSubject),ButtonTime(ListSubject),'UniformOutput',false))*0.03;
    end
end
h = figure(); hold on;
for diff = 1:3
    shadedErrorBar(binboundariessupp(1:end),squeeze(nanmean(hitratediffsupp(:,diff,:),3))'*100,nanstd(squeeze(hitratediffsupp(:,diff,:))')*100/(sqrt(length(ListSubject)))*2,{'Color',colordiff{diff},'LineWidth',1.5},1); ylim([0,100]); ylabel('Detection rate [%]'); xlabel('Time from onset [s]');
end
xlim([binboundariessupp(1),binboundariessupp(end)]) % half of the bin
% plot([3,3],[0,100],'k--'); plot([8,8],[0,100],'k--');
saveas(h,[Rootaddress 'PupilPaper/Supphitrate.pdf'])
saveas(h,[Rootaddress 'PupilPaper/Supphitrate.fig'])
h = figure(); hold on;
for diff = 1:3
    shadedErrorBar(binboundariessupp(1:end),squeeze(nanmean(rtbindiffsupp(:,diff,:),3))',nanstd(squeeze(rtbindiffsupp(:,diff,:))')/(sqrt(length(ListSubject)))*2,{'Color',colordiff{diff},'LineWidth',1.5},1); ylim([0.5,1.5]); ylabel('Median reaction time [s]'); xlabel('Time from onset [s]');
end
xlim([binboundariessupp(1),binboundariessupp(end)]) % half of the bin
% plot([3,3],[0,2],'k--'); plot([8,8],[0,2],'k--');
saveas(h,[Rootaddress 'PupilPaper/Supprt.pdf'])
saveas(h,[Rootaddress 'PupilPaper/Supprt.fig'])


% plot rt dirstibution for large changes as afunction of perforamce
% figure(); hold on;
% bla = [prctile(rtall{1,3},[20,40,60,80]);prctile(rtall{2,3},[20,40,60,80]);prctile(rtall{3,3},[20,40,60,80])];
% for b = [1,2,3]
%     plot(bla(b,:),squeeze(mean(perfblockdiff(b,3,:),3))*100,'LineWidth',2,'Color', colorblockid{b});
% end
%%%%%%%only for figure paper purposes
% close all
% h = figure(1); hold on;
% h2 = figure(2); hold on;
% for b = 1:3
% fablocktiming{b} = cellfun(@(x,y) x(intersect(find(y==BlockValue(b)), find(x<266))),timingfa,timingfablock,'UniformOutput',0);
% for snb = 1:15
% for tt = 1:10
%     fatimebinind(b,snb,tt) = cellfun(@(x) numel(find(x>=tt-1 & x<tt)),{[fablocktiming{b}{snb}]*0.03});
% end
% cumsumfasubject(b,snb,:) = cumsum(squeeze(fatimebinind(b,snb,:)))./max(cumsum(squeeze(fatimebinind(b,snb,:))));
% AUCfacdf(b,snb) = trapz(squeeze(cumsumfasubject(b,snb,:)));
% figure(1); hold on; plot((0:tt),[0;squeeze(cumsumfasubject(b,snb,:))],'LineWidth',0.1,'Color',colorblockid{b});
% figure(2); hold on; 
% plot(b+subjitt(snb),AUCfacdf(b,snb),'.','Color',colorblockid{b},'Markersize',10);
% end
% figure(1); hold on;% plot((0:tt),[0;squeeze(mean(cumsumfasubject(b,:,:),2))],'LineWidth',4,'Color',colorblockid{b});
% xlabel('Time within trial [s]'); ylabel('CDF false alarm'); xlim([0,tt-1])
% figure(2); hold on; 
% hold on;
% N(b) = plot(b,mean(AUCfacdf(b,:)),'.','Color',colorblockid{b},'Markersize',30); ylim([3,8]);xlim([0,4]);
% end
% hold on;
% xticks(1:3)
% xticklabels({'Early','Interm.','Late'})
% ylabel('Average AUC false alarm CDF [0-8s]')
% legend([N(1),N(2),N(3)],{'Early','Interm.','Late'},'Location','SouthWest'); legend boxoff
% [pfablock,~,fatimestat] = friedman(AUCfacdf');
% str = ['p = ', num2str(pfablock)];
% text(1.75,7.75,str)
% saveas(h,[Rootaddress 'PupilPaper/FABlockCDF1.fig'])
% saveas(h2,[Rootaddress 'PupilPaper/FABlockCDF2.fig'])

% 
% % check perf per block as a function of change time [since criterion lowers over time compare 3-5s and 5 to 8s]
% hit35 = cellfun(@(x,y,z) numel(intersect(intersect(find(x == 2),find(y>=100 & y<167)),find(z == 60))),Outcome(ListSubject),ChangeTime(ListSubject),ChangeSize(ListSubject))
% miss35 = cellfun(@(x,y,z) numel(intersect(intersect(find(x == 3),find(y>=100 & y<167)),find(z == 60))),Outcome(ListSubject),ChangeTime(ListSubject),ChangeSize(ListSubject))
% perf35 = hit35./(hit35+miss35);
% hit58 = cellfun(@(x,y,z) numel(intersect(intersect(find(x == 2),find(y>=167 & y<266)),find(z == 60))),Outcome(ListSubject),ChangeTime(ListSubject),ChangeSize(ListSubject))
% miss58 = cellfun(@(x,y,z) numel(intersect(intersect(find(x == 3),find(y>=167 & y<266)),find(z == 60))),Outcome(ListSubject),ChangeTime(ListSubject),ChangeSize(ListSubject))
% perf58 = hit58./(hit58+miss58);

% % Stats for behavior 
%     % first test perf per diff per block with paired anova2 (block*diff)
%     sizedata = size(perfblockdiff,1)*size(perfblockdiff,2)*size(perfblockdiff,3);
%     vectperfblockdiff = reshape(perfblockdiff,1,sizedata);
%     kstest(vectperfblockdiff); % data not normal use friedman per diff
%     blocknb = repmat([1:3]',sizedata/3,1)';
%     diffnb = repmat([ones(1,3)*difflvlind(1),ones(1,3)*difflvlind(2),ones(1,3)*difflvlind(3)],1,sizedata/9);
%     subjectind = reshape(repmat([1:15],sizedata/15,1),1,sizedata);
%    
% %     [a,b,c] = anovan(vectperfblockdiff,{blocknb diffnb},'model',2,'varnames',{'block','diff'})
% 
%     repeatperfblockdiffmat = reshape(vectperfblockdiff,sizedata/9,9); 
%     tperf = table([1:15]',squeeze(perfblockdiff(1,1,:)),squeeze(perfblockdiff(1,2,:)),squeeze(perfblockdiff(1,3,:)),squeeze(perfblockdiff(2,1,:)),squeeze(perfblockdiff(2,2,:)),squeeze(perfblockdiff(2,3,:)),squeeze(perfblockdiff(3,1,:)),squeeze(perfblockdiff(3,2,:)),squeeze(perfblockdiff(3,3,:))...
%         ,'VariableNames',{'Subject','b1d1','b1d2','b1d3','b2d1','b2d2','b2d3','b3d1','b3d2','b3d3'});
%     trt = table([1:15]',squeeze(RTdiffblock(1,1,:)),squeeze(RTdiffblock(1,2,:)),squeeze(RTdiffblock(1,3,:)),squeeze(RTdiffblock(2,1,:)),squeeze(RTdiffblock(2,2,:)),squeeze(RTdiffblock(2,3,:)),squeeze(RTdiffblock(3,1,:)),squeeze(RTdiffblock(3,2,:)),squeeze(RTdiffblock(3,3,:))...
%         ,'VariableNames',{'Subject','b1d1','b1d2','b1d3','b2d1','b2d2','b2d3','b3d1','b3d2','b3d3'});
%     datatable.Properties.VariableNames = {'diff','bloc'};
%     trtnorm = table([1:15]',squeeze(medianallRTdiffblocknorm(1,1,:)),squeeze(medianallRTdiffblocknorm(1,2,:)),squeeze(medianallRTdiffblocknorm(1,3,:)),squeeze(medianallRTdiffblocknorm(2,1,:)),squeeze(medianallRTdiffblocknorm(2,2,:)),squeeze(medianallRTdiffblocknorm(2,3,:)),squeeze(medianallRTdiffblocknorm(3,1,:)),squeeze(medianallRTdiffblocknorm(3,2,:)),squeeze(medianallRTdiffblocknorm(3,3,:))...
%         ,'VariableNames',{'Subject','b1d1','b1d2','b1d3','b2d1','b2d2','b2d3','b3d1','b3d2','b3d3'});
%     datatable.Properties.VariableNames = {'diff','bloc'};
% % % When you have more than one repeated-measures factor, you must set up a table
% % % to indicate the levels on each factor for each of your different variables.
% % % Here is the command you need for this case:
% WithinStructure = table(repmat([1 2 3],1,3)',reshape(repmat([1 2 3],3,1),9,1),'VariableNames',{'diff','bloc'});
% % % The 4 different rows of the WithinStructure table correspond to the 4 different
% % % columns, 'pre_A','pre_B','post_A','post_B', respectively, in your data table.
% % % Each 'pre_A','pre_B','post_A','post_B' column is coded as 1/2 on the PrePost factor
% % % and as 1/2 on the TreatAB factor.
% % % (Added based on later comments:)
% % % Indicate that the 1's and 2's of WithinStructure are category labels
% % % rather than regression-type numerical covariates:
% WithinStructure.bloc = categorical(WithinStructure.bloc);
% WithinStructure.diff = categorical(WithinStructure.diff);    
% rm = fitrm(tperf,'b1d1,b1d2,b1d3,b2d1,b2d2,b2d3,b3d1,b3d2,b3d3~Subject','WithinDesign',WithinStructure);
% rm = fitrm(trt,'b1d1,b1d2,b1d3,b2d1,b2d2,b2d3,b3d1,b3d2,b3d3~Subject','WithinDesign',WithinStructure);
% [ranovatable,b,c] = ranova(rm,'WithinModel','diff*bloc');
% % test for sphericity of the model
% mauchly(rm,c)
% % 
% % % test ranova dprime block time subject
% dtest = table([1:15]',squeeze(dprimerate(:,1,6)),squeeze(dprimerate(:,2,6)),squeeze(dprimerate(:,3,6)),squeeze(dprimerate(:,1,7)),squeeze(dprimerate(:,2,7)),squeeze(dprimerate(:,3,7)),squeeze(dprimerate(:,1,8)),squeeze(dprimerate(:,2,8)),squeeze(dprimerate(:,3,8))...
%     ,squeeze(dprimerate(:,1,9)),squeeze(dprimerate(:,2,9)),squeeze(dprimerate(:,3,9)),squeeze(dprimerate(:,1,10)),squeeze(dprimerate(:,2,10)),squeeze(dprimerate(:,3,10)),squeeze(dprimerate(:,1,11)),squeeze(dprimerate(:,2,11)),squeeze(dprimerate(:,3,11))...
%         ,squeeze(dprimerate(:,1,12)),squeeze(dprimerate(:,2,12)),squeeze(dprimerate(:,3,12)),squeeze(dprimerate(:,1,13)),squeeze(dprimerate(:,2,13)),squeeze(dprimerate(:,3,13)),squeeze(dprimerate(:,1,14)),squeeze(dprimerate(:,2,14)),squeeze(dprimerate(:,3,14))...
%     ,'VariableNames',{'Subject','b1t6','b2t6','b3t6','b1t7','b2t7','b3t7','b1t8','b2t8','b3t8','b1t9','b2t9','b3t9','b1t10','b2t10','b3t10','b1t11','b2t11','b3t11','b1t12','b2t12','b3t12','b1t13','b2t13','b3t13','b1t14','b2t14','b3t14'});
% datatable.Properties.VariableNames = {'bloc','time'};
% WithinStructure = table(repmat([1 2 3],1,9)',reshape(repmat([1:9],3,1),27,1),'VariableNames',{'bloc','time'});
% WithinStructure.bloc = categorical(WithinStructure.bloc);
% WithinStructure.diff = categorical(WithinStructure.time);    
% rm = fitrm(dtest,'b1t6,b2t6,b3t6,b1t7,b2t7,b3t7,b1t8,b2t8,b3t8,b1t9,b2t9,b3t9,b1t10,b2t10,b3t10,b1t11,b2t11,b3t11,b1t12,b2t12,b3t12,b1t13,b2t13,b3t13,b1t14,b2t14,b3t14~Subject','WithinDesign',WithinStructure);
% ranovatable = ranova(rm,'WithinModel','bloc*time');
% 
% ctest = table([1:15]',squeeze(criterion(:,1,6)),squeeze(criterion(:,2,6)),squeeze(criterion(:,3,6)),squeeze(criterion(:,1,7)),squeeze(criterion(:,2,7)),squeeze(criterion(:,3,7)),squeeze(criterion(:,1,8)),squeeze(criterion(:,2,8)),squeeze(criterion(:,3,8))...
%     ,squeeze(criterion(:,1,9)),squeeze(criterion(:,2,9)),squeeze(criterion(:,3,9)),squeeze(criterion(:,1,10)),squeeze(criterion(:,2,10)),squeeze(criterion(:,3,10)),squeeze(criterion(:,1,11)),squeeze(criterion(:,2,11)),squeeze(criterion(:,3,11))...
%         ,squeeze(criterion(:,1,12)),squeeze(criterion(:,2,12)),squeeze(criterion(:,3,12)),squeeze(criterion(:,1,13)),squeeze(criterion(:,2,13)),squeeze(criterion(:,3,13)),squeeze(criterion(:,1,14)),squeeze(criterion(:,2,14)),squeeze(criterion(:,3,14))...
%     ,'VariableNames',{'Subject','b1t6','b2t6','b3t6','b1t7','b2t7','b3t7','b1t8','b2t8','b3t8','b1t9','b2t9','b3t9','b1t10','b2t10','b3t10','b1t11','b2t11','b3t11','b1t12','b2t12','b3t12','b1t13','b2t13','b3t13','b1t14','b2t14','b3t14'});
% datatable.Properties.VariableNames = {'bloc','time'};
% WithinStructure = table(repmat([1 2 3],1,9)',reshape(repmat([1:9],3,1),27,1),'VariableNames',{'bloc','time'});
% WithinStructure.bloc = categorical(WithinStructure.bloc);
% WithinStructure.diff = categorical(WithinStructure.time);    
% rm = fitrm(ctest,'b1t6,b2t6,b3t6,b1t7,b2t7,b3t7,b1t8,b2t8,b3t8,b1t9,b2t9,b3t9,b1t10,b2t10,b3t10,b1t11,b2t11,b3t11,b1t12,b2t12,b3t12,b1t13,b2t13,b3t13,b1t14,b2t14,b3t14~Subject','WithinDesign',WithinStructure);
% ranovatable = ranova(rm,'WithinModel','bloc*time');



end