%runEmpiricalSurnameStats.m
%EW Tekwa May 24, 2022

set(0,'DefaultAxesFontSize',20)
scrsz = get(0,'ScreenSize');
set(0,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);

% %set up data matrics (rows: entrance LR, advancement ratio, publication
% %ratio relative to the global population)
% HighIncomeUS=[];
% LowIncomeUS=[];
%
% HarvardEnvSci=[];
% HarvardImmun=[];
% HarvardSocio=[];
% HarvardEvol=[];
% HarvardEthni=[];
% HarvardComSci=[];

% HarvardSci=[];
% HarvardSocioSci=[];
% HarvardMed=[];
% HarvardBus=[];
% HarvardPubHea=[];
% HarvardArt=[];

%load data:
load('EmpiricalSurnameData.mat') %ratios with global population as reference

USEntranceLR=10.31;
USAdvancementR=0.45;
USPubR=4.65;

EntLRNorm=[1 USEntranceLR];
AdvRNorm=[1 USAdvancementR];
PubNorm=[1 USPubR];

NormType=2; %1 for global reference, 2 for US reference

%overwrite ratio data with chosen normalization:
HighIncomeUS=HighIncomeUS.*[EntLRNorm(3-NormType) AdvRNorm(3-NormType) PubNorm(3-NormType)];
LowIncomeUS=LowIncomeUS.*[EntLRNorm(3-NormType) AdvRNorm(3-NormType) PubNorm(3-NormType)];

HarvardEnvSci=HarvardEnvSci./[EntLRNorm(NormType) AdvRNorm(NormType) PubNorm(NormType) 1];
HarvardImmun=HarvardImmun./[EntLRNorm(NormType) AdvRNorm(NormType) PubNorm(NormType) 1];
HarvardSocio=HarvardSocio./[EntLRNorm(NormType) AdvRNorm(NormType) PubNorm(NormType) 1];
HarvardEvol=HarvardEvol./[EntLRNorm(NormType) AdvRNorm(NormType) PubNorm(NormType) 1];
HarvardEthni=HarvardEthni./[EntLRNorm(NormType) AdvRNorm(NormType) PubNorm(NormType) 1];
HarvardComSci=HarvardComSci./[EntLRNorm(NormType) AdvRNorm(NormType) PubNorm(NormType) 1];

%get regressions of publications & rel pub on LR:
HarvardFocusData={HarvardComSci,HarvardEnvSci,HarvardEthni,HarvardEvol,HarvardImmun,HarvardSocio};
HarvardFocusUnit={'ComSci','EnvSci','Ethni','Evol','Immun','Socio'};
for plotType=[4]
    figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.5 scrsz(4)/2]);
    for i=1:6
        mdlFocus=fitlm(log10(HarvardFocusData{i}(:,1)+1),log10(HarvardFocusData{i}(:,plotType)+1));
        subplot(2,3,i)
        pl=plot(mdlFocus);
        set(pl,'LineWidth',2,'Color','k')
        title(HarvardFocusUnit{i},'fontweight','normal')
        if i>3
            xlabel 'log(LR+1)'
        else
            xlabel ''
        end
        if i==1 || i==4
            label_h=ylabel('log(pub+1)');
        else
            label_h=ylabel('');
        end
        legend off
        getX=xlim;
        text((getX(2)-getX(1))*0.3+getX(1),0.5,{['\beta=' num2str(mdlFocus.Coefficients.Estimate(2),2) '\pm' num2str(mdlFocus.Coefficients.SE(2),2)];['p=' num2str(mdlFocus.Coefficients.pValue(2),2) ' (n=' num2str(mdlFocus.NumObservations) ')']},'Fontsize',16);
        ylim([0,3])
        ylims=ylim;
        text(label_h.Position(1),ylims(2)+diff(ylims)*0.11,char(64+i),'Fontsize',20)
    end
end

HarvardArt=HarvardArt*EntLRNorm(3-NormType);
HarvardSocioSci=HarvardSocioSci*EntLRNorm(3-NormType);
HarvardSci=HarvardSci*EntLRNorm(3-NormType);
HarvardMed=HarvardMed*EntLRNorm(3-NormType);
HarvardPubHea=HarvardPubHea*EntLRNorm(3-NormType);
HarvardBus=HarvardBus*EntLRNorm(3-NormType);

%compute entrance LR and advancement ratio for groups
% entRange=linspace(0,maxEntAll,101);
% advRange=linspace(0,maxAdvAll,101);
if NormType==1
    maxEntAll=60; % 5
    maxAdvAll=3.5; % 6
    maxEntPAll=0.6; % 0.35
    maxAdvPAll=0.6; % 0.35
else
    maxEntAll=5;
    maxAdvAll=6;
    maxEntPAll=0.35;
    maxAdvPAll=0.35;
end
maxPubAll=12;
entRange=[linspace(0,maxEntAll,21) Inf];
advRange=[linspace(0,maxAdvAll,21) Inf];

LowIncomeUSEnt=histcounts(LowIncomeUS(:,1),entRange)/length(LowIncomeUS);
HighIncomeUSEnt=histcounts(HighIncomeUS(:,1),entRange)/length(HighIncomeUS);
LowIncomeUSAdv=histcounts(LowIncomeUS(:,2),advRange)/length(LowIncomeUS);
HighIncomeUSAdv=histcounts(HighIncomeUS(:,2),advRange)/length(HighIncomeUS);
minEntP=zeros(1,6);
maxEntP=zeros(1,6);
minEnt=zeros(1,6);
maxAdv=zeros(1,6);
minAdvP=zeros(1,6);
maxAdvP=zeros(1,6);
minAdv=zeros(1,6);
maxAdv=zeros(1,6);

%figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/5 scrsz(4)/3]);
figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)/2.5]);
colorSeq=['y','g','k','w','c','m','b','r']; %colour for ComSci, EnvSci, Ethi, Evol, Immun, Socio, HighIncome, LowIncome
%HarvardData={HarvardComSci,HarvardEnvSci,HarvardEthni,HarvardEvol,HarvardImmun,HarvardSocio};
%HarvardUnit={'ComSci','EnvSci','Ethni','Evol','Immun','Socio'};
HarvardData={HarvardArt,HarvardSocioSci,HarvardSci,HarvardMed,HarvardPubHea,HarvardBus};
HarvardUnit={'Art','Soc','Sci','Med','Pub','Bus'};
for i=1:length(HarvardData)
    if ~isempty(HarvardData{i})
        subplot(1,11,i+4)
        HarvardEnt=histcounts(HarvardData{i}(:,1),entRange)/length(HarvardData{i});
        [h,p,ci,stats] = ttest(HarvardData{i}(:,1),1,'Tail','right');
        [hHi,pHi,ciHi,statsHi] = ttest2(HarvardData{i}(:,1),HighIncomeUS(:,1));
        hold on
        yline(mean(HighIncomeUS(:,1)),'r')
        yline(1)
        barh(entRange(1:end-1)+0.5*entRange(2),HarvardEnt,1,colorSeq(i),'FaceAlpha',0.2)
        scatter(0.055,mean(HarvardData{i}(:,1)),150,'>','filled','MarkerFaceColor',colorSeq(i),'MarkerEdgeColor','k','LineWidth',2)
        ylabel 'entrance LR'
        box off
        minEntP(i)=min(HarvardEnt);
        maxEntP(i)=max(HarvardEnt);
        minEnt(i)=min(HarvardData{i}(:,1));
        maxEnt(i)=max(HarvardData{i}(:,1));
%         if i==1
%             set(get(gca,'XAxis'),'visible','off')
%             curLabels=yticklabels;
%             curLabels{end}=['>' num2str(curLabels{end})];
%             yticklabels(curLabels)
%             xlim([-maxEntPAll*1.7,0])
%             ax = gca;
%             ax.Position(1) = ax.Position(1)*(1-1.1/1.7*0.25);
%             ax.Position(3) = (1.7/1.1)*ax.Position(3);
%         else
%             axis off
%             xlim([-maxEntPAll*1.1,0])
%         end
        axis off
        title(HarvardUnit{i},'fontweight','normal')
        xlim([0,maxEntPAll*1.1])
        ylim([-0.5,maxEntAll*1.1])
        text(0.09,maxEntAll*0.55,['p_1=' num2str(p,0)],'Fontsize',14)
        text(0.02,maxEntAll*0.75,['p_{hi}=' num2str(pHi,0)],'Fontsize',14,'Color','r')
    end
end

subplot(1,11,1:4)
[hLo,pLo,ciLo,statsLo] = ttest2(HighIncomeUS(:,1),LowIncomeUS(:,1));
[h,p,ci,stats] = ttest(HighIncomeUS(:,1),1,'Tail','right'); %p value for HighIncome LR greater than one
% yyaxis left
% yticks([])
% yticklabels([])
% yyaxis right
hold on
% [f_HighIncomeUS_ent,xi_HighIncomeUS_ent]=ksdensity(HighIncomeUS(:,1)); %entrance
% [f_HighIncomeUS_adv,xi_HighIncomeUS_adv]=ksdensity(HighIncomeUS(:,2)); %advancement
% [f_LowIncomeUS_ent,xi_LowIncomeUS_ent]=ksdensity(LowIncomeUS(:,1)); %entrance
% [f_LowIncomeUS_adv,xi_LowIncomeUS_adv]=ksdensity(LowIncomeUS(:,2)); %advancement
% plot(xi_HighIncomeUS_ent,f_HighIncomeUS_ent,'-r','LineWidth',2);
% plot(xi_LowIncomeUS_ent,-f_LowIncomeUS_ent,'-b','LineWidth',2);
yline(mean(HighIncomeUS(:,1)),'r')
yline(1)
barh(entRange(1:end-1)+0.5*entRange(2),LowIncomeUSEnt,1,colorSeq(7),'FaceAlpha',0.5)
barh(entRange(1:end-1)+0.5*entRange(2),-HighIncomeUSEnt,1,colorSeq(8),'FaceAlpha',0.5)
scatter(-0.07,mean(HighIncomeUS(:,1)),150,'>','filled','MarkerEdgeColor','k','MarkerFaceColor',colorSeq(8),'LineWidth',2)
scatter(0.06,mean(LowIncomeUS(:,1)),150,'>','filled','MarkerEdgeColor','k','MarkerFaceColor',colorSeq(7),'LineWidth',2)
%ylabel 'entrance LR'
%xlabel 'probability'
%box on
curLabels=yticklabels;
curLabels{end}=['>' num2str(curLabels{end})];
yticklabels(curLabels)
xlim([-max(HighIncomeUSEnt)*1.1,max(LowIncomeUSEnt)])
ylim([-0.5,maxEntAll*1.1])
set(get(gca,'XAxis'),'visible','off')
ylabel 'entrance LR'
title('High/low income','fontweight','normal')
text(0.2,maxEntAll*0.7,{['hi/lo ratio=' num2str(mean(HighIncomeUS(:,1))/mean(LowIncomeUS(:,1)),2)];['p_{hi/lo}=' num2str(pLo,0)];['p_{hi1}=' num2str(p,0)]},'Fontsize',14)
% 
% subplot(2,8,9:10)
% [hLo,pLo,ciLo,statsLo] = ttest2(HighIncomeUS(:,2),LowIncomeUS(:,2));
% [h,p,ci,stats] = ttest(HighIncomeUS(:,2),1,'Tail','right'); %p value for HighIncome LR greater than one
% hold on
% barh(advRange(1:end-1)+0.5*advRange(2),LowIncomeUSAdv,1,colorSeq(7),'FaceAlpha',0.5)
% barh(advRange(1:end-1)+0.5*advRange(2),-HighIncomeUSAdv,1,colorSeq(8),'FaceAlpha',0.5)
% yline(1)
% scatter(-0.035,mean(HighIncomeUS(:,2)),150,'>','filled','MarkerEdgeColor','k','MarkerFaceColor',colorSeq(8),'LineWidth',2)
% scatter(0.025,mean(LowIncomeUS(:,2)),150,'>','filled','MarkerEdgeColor','k','MarkerFaceColor',colorSeq(7),'LineWidth',2)
% %ylabel 'advancement ratio'
% %xlabel 'probability'
% %box on
% curLabels=yticklabels;
% curLabels{end}=['>' num2str(curLabels{end})];
% yticklabels(curLabels)
% xlim([-max(HighIncomeUSAdv)*1.1,max(LowIncomeUSAdv)*1.1])
% ylim([0,maxAdvAll*1.1])
% % curLabels=yticklabels;
% % curLabels{end+1}=['>' num2str(maxAdvAll)];
% % yticklabels(curLabels)
% h = gca;
% h.XAxis.Visible = 'off';
% ylabel 'advancement ratio'
% %title('High/low income','fontweight','normal')
% text(0.04,maxAdvAll*0.8,{['hi/lo ratio=' num2str(mean(HighIncomeUS(:,2))/mean(LowIncomeUS(:,2)),2)];['p_{hi/lo}=' num2str(pLo,0)];['p_{hi1}=' num2str(p,0)]},'Fontsize',14)
% 
% for i=1:length(HarvardData)
%     if ~isempty(HarvardData{i})
%         subplot(2,8,10+i)
%         HarvardAdv=histcounts(HarvardData{i}(:,2),advRange)/length(HarvardData{i});
%         [h,p,ci,stats] = ttest(HarvardData{i}(:,2),1,'Tail','right');
%         [hHi,pHi,ciHi,statsHi] = ttest2(HarvardData{i}(:,2),HighIncomeUS(:,2));
%         hold on
%         barh(advRange(1:end-1)+0.5*advRange(2),HarvardAdv,1,colorSeq(i),'FaceAlpha',0.2)
%         yline(1)
%         scatter(0.03,mean(HarvardData{i}(:,2)),150,'>','filled','MarkerFaceColor',colorSeq(i),'MarkerEdgeColor','k','LineWidth',2)
%         box off
%         minAdvP(i)=min(HarvardAdv);
%         maxAdvP(i)=max(HarvardAdv);
%         minAdv(i)=min(HarvardData{i}(:,2));
%         maxAdv(i)=max(HarvardData{i}(:,2));
%         %title(HarvardUnit{i},'fontweight','normal')
% %         if i==1
% %             yyaxis left
% %             yticks([])
% %             yticklabels([])
% %             yyaxis right
% %             set(get(gca,'XAxis'),'visible','off')
% %             curLabels=yticklabels;
% %             curLabels{end}=['>' num2str(curLabels{end})];
% %             yticklabels(curLabels)
% %             ylabel 'advancement ratio'
% %             xlim([0,maxAdvPAll*1.7])
% %             ax = gca;
% %             ax.Position(3) = (1.7/1.1)*ax.Position(3);
% %         else
% %             axis off
% %             xlim([0,maxAdvPAll*1.1])
% %         end
%         axis off
%         xlim([0,maxAdvPAll*1.1])
%         ylim([0,maxAdvAll*1.1])
%         text(0.05,maxAdvAll*0.2,['p_1=' num2str(p,0)],'Fontsize',14)
%         text(0.12,maxAdvAll*0.7,['p_{hi}=' num2str(pHi,0)],'Fontsize',14)
%     end
% end

% maxEntAll=max(maxEnt)
% maxAdvAll=max(maxAdv)
minEntPAll=min(minEntP)
maxEntPAll=max(maxEntP)
% minAdvPAll=min(minAdvP)
% maxAdvPAll=max(maxAdvP)


%figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2.8 scrsz(4)/5.5]);
allHarvardAcademics=[];
for i=1:length(HarvardData)
    allHarvardAcademics=[allHarvardAcademics; HarvardData{i}];
end
%subplot(1,3,1)
subplot(1,11,i+5)
HarvardEnt=histcounts(allHarvardAcademics(:,1),entRange)/length(allHarvardAcademics(:,1));
[h,p,ci,stats] = ttest(allHarvardAcademics(:,1),1,'Tail','right');
[hHi,pHi,ciHi,statsHi] = ttest2(allHarvardAcademics(:,1),HighIncomeUS(:,1));
hold on
yline(mean(HighIncomeUS(:,1)),'r')
yline(1)
barh(entRange(1:end-1)+0.5*entRange(2),HarvardEnt,1,'r','FaceAlpha',0.9)
scatter(0.055,mean(allHarvardAcademics(:,1)),150,'>','filled','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2)
axis off
title('Harvard (all)','fontweight','normal')
xlim([0,maxEntPAll*1.1])
ylim([-0.5,maxEntAll*1.1])
text(0.09,maxEntAll*0.55,['p_1=' num2str(p,0)],'Fontsize',14)
text(0.02,maxEntAll*0.75,['p_{hi}=' num2str(pHi,0)],'Fontsize',14,'Color','r')
        
%histogram(allHarvardAcademics(:,1),[linspace(0,maxEntAll,21) Inf],'FaceColor','b','FaceAlpha',0.2)
% xline(1)
% xline(mean(allHarvardAcademics(:,1)),'linewidth',2)
% curX=xlim;
% curY=ylim;
% text(curX(2)*0.6,curY(2)*0.8,['p_1=' num2str(p,0)],'Fontsize',14);
% box off
% ylabel 'frequency'
% xlabel 'entrance LR'
% subplot(1,3,2)
% [h,p,ci,stats] = ttest(allHarvardAcademics(:,2),1,'Tail','right');
% hold on
% histogram(allHarvardAcademics(:,2),[linspace(0,maxAdvAll,21) Inf],'FaceColor','c','FaceAlpha',0.2)
% xline(1)
% xline(mean(allHarvardAcademics(:,2)),'linewidth',2)
% curX=xlim;
% curY=ylim;
% text(curX(2)*0.6,curY(2)*0.8,['p_1=' num2str(p,0)],'Fontsize',14);
% box off
% xlabel 'advancement ratio'
% subplot(1,3,3)
% [h,p,ci,stats] = ttest(allHarvardAcademics(:,3),1,'Tail','right');
% hold on
% %bar(edges(1:end-1)+0.5*edges(2),HarvardPub,1,'k','FaceAlpha',0.8)
% histogram(allHarvardAcademics(:,3),[linspace(0,maxPubAll,21) Inf],'FaceColor','k','FaceAlpha',0.2)
% xline(1)
% xline(mean(allHarvardAcademics(:,3)),'linewidth',2)
% curX=xlim;
% curY=ylim;
% text(curX(2)*0.6,curY(2)*0.8,['p_1=' num2str(p,0)],'Fontsize',14);
% %scatter(0.03,mean(HarvardData{i}(:,2)),150,'>','filled','MarkerFaceColor',colorSeq(i),'MarkerEdgeColor','k','LineWidth',2)
% box off
% xlabel 'publication ratio'

means=mean(allHarvardAcademics)
stds=std(allHarvardAcademics)
numEntBelow1=sum(allHarvardAcademics(:,1)<1)
% numAdvBelow1=sum(allHarvardAcademics(:,2)<1)
% numPubBelow1=sum(allHarvardAcademics(:,3)<1)

%power analysis based on aggregate Harvard data:
powerHarvard = sampsizepwr('t',[1 stds(1)],means(1),[],length(allHarvardAcademics),'Tail','right');
betaHarvard = 1 - powerHarvard %probability of false negative (falsely rejecting alternative hypothesis)
nHarvard = sampsizepwr('t',[1 stds(1)],means(1),0.9,[],'Tail','right') %estimated sample size required to reject null at alpha=0.05 with a power of 0.9


meanHiIncome=mean(HighIncomeUS(:,1))
stdHiIncome=std(HighIncomeUS(:,1))
powerHiIncome = sampsizepwr('t',[1 stdHiIncome],meanHiIncome,[],length(HighIncomeUS),'Tail','right');
betaHiIncome = 1 - powerHiIncome %probability of false negative (falsely rejecting alternative hypothesis)
nHiIncome = sampsizepwr('t',[1 stdHiIncome],meanHiIncome,0.9,[],'Tail','right')  %estimated sample size required to reject null at alpha=0.05 with a power of 0.9

meanLowIncome=mean(LowIncomeUS(:,1))
stdLowIncome=std(LowIncomeUS(:,1))