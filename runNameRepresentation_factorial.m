%run and plot a series of nameRepresentation scenarios
set(0,'DefaultAxesFontSize',14)
scrsz = get(0,'ScreenSize');
fontS=14; %font size
%figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)]);
%rngSeed=50;
%rng(rngSeed); %set random number generator seed (1,10,100,1000)
rngSeed=rng('shuffle')
time=datetime;
%baseline values for parameters that will vary across scenarios:

numReps=1
evaluations=[0,1/3,2/3,1]; %capital evaluation weight
mateChoice_all=[0:0.4:0.8]; %capital mating weight
merit_herit_all=[0:0.4:0.8]; %merit heritability
capital_herit_all=[0:0.4:0.8]; %capital heritability
capital_added_all=[0:25:50]; %capital added to academic children
MeritCapCausal_all=[0:0.4:0.8]; %merit-capital slope (from outside of academia)
CV_all=[0.5 1]; %variance over mean ratio (for merit, capital, and capital added)
AcademiaPorp_all=[0.02]; %percent of people in academia
% evaluations=[0:0.25:1];
% mateChoice_all=[0:0.25:0.75];
% merit_herit_all=[0:0.25:1];
% capital_herit_all=[0:0.25:1];
% capital_added_all=[0:25:50];
numScenarios=length(evaluations)*length(mateChoice_all)*length(merit_herit_all)*length(capital_herit_all)*length(capital_added_all)*length(MeritCapCausal_all)*length(CV_all)*length(AcademiaPorp_all)
numScenarios=numScenarios*numReps;


%CV=1; %variance over mean ratio
numGen=20; %number of generations to simulate (20)
numNames=1000; %initial number of names (2000)
initPopPerName=50; %per name mean population (238 globally, 63 in US)
nameDistr=2; %name frequency distribution is Poisson (1) or Exponential (2)
%AcademiaPorp=0.002; %percent of people in academia, 0.0022 global, 0.002, 0.0004
%initAcademiaPop=round(numNames*initPopPerName*AcademiaPorp); %number of people in academia
reprodRate=2; %mean children per pair of parents per generation, actual is Poisson random
reprodCost1=0; %reduction in numbeer of children per generation for couples with one academics
reprodCost2=reprodCost1; %reduction in numbeer of children per generation for couples with two academics
%evaluation=1; %evaluation for academic entrance based on only merit (0) to only capital (1)
merit_mean=100;
%merit_var=(merit_mean*CV)*(merit_mean>0);
capital_mean=100;
%capital_var=(capital_mean*CV)*(capital_mean>0);
% mateChoice=1; %mate choice random (0) or ranked by capital (1) or ranked by merit (2)
% mateChoiceRelVar=0.1*(mateChoice>0); %mate choice ranking variability scaled to score of 100 (this variability is added to true capital or capital+merit scores before ranking)
%merit_herit=0.5; %0 to 1
%capital_herit=0.9;
%capital_added=0; %mean capital added to children of academics
%capital_added_var=(capital_added/2)*(capital_added>0); %variance in captial added to children of academics
child_surname=0; %Child inherits surname of parent with higher capital (1), or with higher capital+merit (2), or random (0)
child_surnameRelVar=0.1*(child_surname>0); %name convention ranking variability scaled to score of 100 (this variability is added to true capital or capital+merit scores before ranking)
nameMutation=0.002; %probability of surname mutation per child per generation

%vectors to store outcomes across scenarios:
authorLikeRatio=zeros(1,numScenarios);
topAcadLikeRatio=zeros(1,numScenarios);
botAcadLikeRatio=zeros(1,numScenarios);
topCapLikeRatio=zeros(1,numScenarios);
botCapLikeRatio=zeros(1,numScenarios);
Port_top=zeros(1,numScenarios);
rankMeritSlope=zeros(1,numScenarios);
capMeritSlope=zeros(1,numScenarios);
capMeritCor=zeros(1,numScenarios);
LRMeritSlope=zeros(1,numScenarios);
Merit_pop=zeros(1,numScenarios);
Merit_top=zeros(1,numScenarios);
Merit_acad=zeros(1,numScenarios);
Capital_pop=zeros(1,numScenarios);
Capital_top=zeros(1,numScenarios);
Capital_acad=zeros(1,numScenarios);
ParamValues=zeros(8,numScenarios);

i=1; %scenario/rep index
for j=1:length(evaluations)
    for k=1:length(mateChoice_all)
        for l=1:length(merit_herit_all)
            for m=1:length(capital_herit_all)
                for n=1:length(capital_added_all)
                    for o=1:length(MeritCapCausal_all)
                        for p=1:length(CV_all)
                            for q=1:length(AcademiaPorp_all)
                                for r=1:numReps
                                    evaluation=evaluations(j);
                                    mateChoice=mateChoice_all(k);
                                    merit_herit=merit_herit_all(l);
                                    capital_herit=capital_herit_all(m);
                                    capital_added=capital_added_all(n);
                                    MeritCapCausal=MeritCapCausal_all(o);
                                    CV=CV_all(p);
                                    AcademiaPorp=AcademiaPorp_all(q);
                                    ParamValues(:,i)=[round(evaluation,2),mateChoice,merit_herit,capital_herit,capital_added,MeritCapCausal,CV,AcademiaPorp]';
                                    
                                    capital_added_var=(capital_added*CV)*(capital_added>0); %variance in captial added to children of academics
                                    merit_var=(merit_mean*CV)*(merit_mean>0);
                                    capital_var=(capital_mean*CV)*(capital_mean>0);
                                    initAcademiaPop=round(numNames*initPopPerName*AcademiaPorp); %number of people in academia
                                    
                                    display(['running scenario ' num2str(ceil(i/numReps)) ' of ' num2str(numScenarios/numReps) ', rep # ' num2str(r)])
                                    %subplot(3,numScenarios/3,i)
                                    [authorLikeRatio(i),topAcadLikeRatio(i),botAcadLikeRatio(i),topCapLikeRatio(i),botCapLikeRatio(i),Port_top(i),rankMeritSlope(i),capMeritSlope(i),capMeritCor(i),LRMeritSlope(i),Merit_pop(i),Merit_top(i),Merit_acad(i),Capital_pop(i),Capital_top(i),Capital_acad(i)]=nameRepresentation_noPlot(numGen,numNames,initPopPerName,nameDistr,initAcademiaPop,reprodRate,reprodCost1,reprodCost2,evaluation,merit_mean,merit_var,capital_mean,capital_var,MeritCapCausal,mateChoice,merit_herit,capital_herit,capital_added,capital_added_var,child_surname,child_surnameRelVar,nameMutation);
                                    i=i+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

save(['allStats_capitalSurname' char(time)]);

% %merge data:
% %-------------
% % % %first, copy data from current set
% authorLikeRatio2=authorLikeRatio;
% topAcadLikeRatio2=topAcadLikeRatio;
% botAcadLikeRatio2=botAcadLikeRatio;
% topCapLikeRatio2=topCapLikeRatio;
% botCapLikeRatio2=botCapLikeRatio;
% Port_top2=Port_top;
% rankMeritSlope2=rankMeritSlope;
% capMeritSlope2=capMeritSlope;
% capMeritCor2=capMeritCor;
% LRMeritSlope2=LRMeritSlope;
% Merit_pop2=Merit_pop;
% Merit_top2=Merit_top;
% Merit_acad2=Merit_acad;
% Capital_pop2=Capital_pop;
% Capital_top2=Capital_top;
% Capital_acad2=Capital_acad;
% ParamValues2=ParamValues;
% 
% authorLikeRatio3=authorLikeRatio2;
% topAcadLikeRatio3=topAcadLikeRatio2;
% botAcadLikeRatio3=botAcadLikeRatio2;
% topCapLikeRatio3=topCapLikeRatio2;
% botCapLikeRatio3=botCapLikeRatio2;
% Port_top3=Port_top2;
% rankMeritSlope3=rankMeritSlope2;
% capMeritSlope3=capMeritSlope2;
% capMeritCor3=capMeritCor2;
% LRMeritSlope3=LRMeritSlope2;
% Merit_pop3=Merit_pop2;
% Merit_top3=Merit_top2;
% Merit_acad3=Merit_acad2;
% Capital_pop3=Capital_pop2;
% Capital_top3=Capital_top2;
% Capital_acad3=Capital_acad2;
% ParamValues3=ParamValues2;
% % %ParamValues2=[ParamValues;[AcademiaPorp_all(1)*ones(1,length(ParamValues)/2) AcademiaPorp_all(2)*ones(1,length(ParamValues)/2)]]; %add AcademiaPorp
% % % %
% % % % %second, load the other data set and merge
% authorLikeRatio=[authorLikeRatio authorLikeRatio2];
% topAcadLikeRatio=[topAcadLikeRatio topAcadLikeRatio2];
% botAcadLikeRatio=[botAcadLikeRatio botAcadLikeRatio2];
% topCapLikeRatio=[topCapLikeRatio topCapLikeRatio2];
% botCapLikeRatio=[botCapLikeRatio botCapLikeRatio2];
% Port_top=[Port_top Port_top2];
% rankMeritSlope=[rankMeritSlope rankMeritSlope2];
% capMeritSlope=[capMeritSlope capMeritSlope2];
% capMeritCor=[capMeritCor capMeritCor2];
% LRMeritSlope=[LRMeritSlope LRMeritSlope2];
% Merit_pop=[Merit_pop Merit_pop2];
% Merit_top=[Merit_top Merit_top2];
% Merit_acad=[Merit_acad Merit_acad2];
% Capital_pop=[Capital_pop Capital_pop2];
% Capital_top=[Capital_top Capital_top2];
% Capital_acad=[Capital_acad Capital_acad2];
% %ParamValues=[ParamValues;repmat(AcademiaPorp_all,1,length(ParamValues)/2)]; %add AcademiaPorp
% ParamValues=[ParamValues ParamValues2];
% % % %
% 
% authorLikeRatio=[authorLikeRatio authorLikeRatio3];
% topAcadLikeRatio=[topAcadLikeRatio topAcadLikeRatio3];
% botAcadLikeRatio=[botAcadLikeRatio botAcadLikeRatio3];
% topCapLikeRatio=[topCapLikeRatio topCapLikeRatio3];
% botCapLikeRatio=[botCapLikeRatio botCapLikeRatio3];
% Port_top=[Port_top Port_top3];
% rankMeritSlope=[rankMeritSlope rankMeritSlope3];
% capMeritSlope=[capMeritSlope capMeritSlope3];
% capMeritCor=[capMeritCor capMeritCor3];
% LRMeritSlope=[LRMeritSlope LRMeritSlope3];
% Merit_pop=[Merit_pop Merit_pop3];
% Merit_top=[Merit_top Merit_top3];
% Merit_acad=[Merit_acad Merit_acad3];
% Capital_pop=[Capital_pop Capital_pop3];
% Capital_top=[Capital_top Capital_top3];
% Capital_acad=[Capital_acad Capital_acad3];
% %ParamValues=[ParamValues;repmat(AcademiaPorp_all,1,length(ParamValues)/2)]; %add AcademiaPorp
% ParamValues=[ParamValues ParamValues3];
% % %-----------


% predictors=[authorLikeRatio;authorLikeRatio;authorLikeRatio;authorLikeRatio;topAcadLikeRatio./botAcadLikeRatio;topCapLikeRatio./botCapLikeRatio];
% responses=[Merit_acad;Port_top;rankMeritSlope;capMeritSlope;Merit_acad;Merit_acad];
effectGroup=find(ParamValues(8,:)>=0.002);
predictors=[topAcadLikeRatio(effectGroup);topAcadLikeRatio(effectGroup)];
responses=[Merit_acad(effectGroup);Port_top(effectGroup)];
% predictors=[topAcadLikeRatio];
% responses=[Merit_acad];
colors=['r','b','g','m'];
for type=1:size(responses,1) %plot academic merit or portion of top merit people in academia
    %figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)/1.5]);
    figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/2.7]);
    %titles={{'0.4 capital mating weight &';'0.4 capital evaluation weight'},'capital evaluation weight','capital mating weight','merit heritability','capital heritability','academic children capital gain'};
    titles={'capital evaluation weight','capital mating weight','merit heritability','capital heritability','academic children capital gain','general merit-capital slope','coeff. var','portion in academia'};
    %intermediateGroup=find(ParamValues(8,:)==0.02);
    for pt=1:8
        subplot(2,4,pt)
        hold on
        %         if pt==1
        %             scatter(predictors(type,:),responses(type,:),14,'k','filled');
        %             [x,I]=sort(predictors(type,intermediateGroup));
        %             y=responses(type,intermediateGroup);
        %             y=y(I);
        %             scatter(x,y,20,'r','filled');
        %             fitobj=fit(x',y','exp2');
        %             %predobj=predint(fitobj,x');
        %             pl=plot(x',feval(fitobj,x'),'r');
        %             set(pl,'LineWidth',2);
        %             %plot(x,predobj,'r--');
        %         else
        groupVals=unique(ParamValues(pt,effectGroup));
        %colors = lines(length(groupVals));
        %colors=colormap(hsv(length(groupVals)));
        gscatter(predictors(type,:),responses(type,:),ParamValues(pt,effectGroup),colors,'',1);
        for group=1:length(groupVals)
            groupIDs=find(ParamValues(pt,effectGroup)==groupVals(group));
            [x,I]=sort(predictors(type,groupIDs));
            I0=find(x>0);
            x=x(I0);
            I=I(I0);
            trim=0; %a priori trim extreme values in predictor
            x=x(trim+1:end-trim);
            I=I(trim+1:end-trim);
            y=responses(type,groupIDs);
            y=y(I);
            miny=min(y);
            maxy=max(y);
            [fitobj,gof,output]=fit(x',y','exp2');
            predy=feval(fitobj,x');
            inty=predint(fitobj,x',0.95,'functional','off');
            if isnan(max(diff(inty')))
                [fitobj,gof,output]=fit(x',y','exp1');
                predy=feval(fitobj,x');
                inty=predint(fitobj,x',0.95,'functional','off');
            end
            if min(diff(inty'))<eps
                [fitobj,gof,output]=fit(x',y','poly1');
                predy=feval(fitobj,x');
                inty=predint(fitobj,x',0.95,'functional','off');
            end
            [diffSize goodI]=find(diff(inty')<(maxy-miny)/2); %find x or y indices where interval size is sufficiently small
            %plot(x',feval(fitobj,x'),'LineWidth',2,'Color',colors(group,:));
            %set(pl,'LineWidth',2,'Color',colors(group,:));
            %plot(x',inty','--','LineWidth',1,'Color',colors(group,:));
            bl1=boundedline(x(goodI)',predy(goodI),[predy(goodI)-inty(goodI,1) inty(goodI,2)-predy(goodI)],colors(group),'alpha'); set(bl1,'linewidth',2);
        end
        hleg = legend ('show');
        hleg.String(length(groupVals)+1:end)=[];
        %         end
        ax=gca;
        ax.XAxis.Scale='log';
        if type<6
            xlabel 'top academia^\primes author likelihood (LR)'
            %xticks([1+eps,2,4,8,16,32])
            %xlim([1,32])
        else
            xlabel 'top academia^\primes author likelihood (LR)'
            %xlim([0.125,160])
        end
%         xticks([0.25 1 1.414 2 2.828 4 16 64 256])
%         ax.XMinorTick='off';
%         xticklabels([0.25 1 1.4 2 2.8 4 16 64 256])
        xticks([0.25 1 2 4 8 16 32 64 256])
        ax.XMinorTick='off';
        xticklabels([0.25 1 2 4 8 16 32 64 256])

        title(titles{pt},'fontweight','normal')
        %         if pt==1
        %             legend({'all','focal group'},'Location','northeast')
        %             %         elseif pt==3
        %             %             legend({'random','90% by capital'},'Location','northeast')
        %         end
        legend('AutoUpdate', 'off')
        if type==1 || type==6
            ylim([97,148]) %95 145
            yline(100);
            label_h=ylabel('academia^\primes mean merit');
        elseif type==2 || type==7
            ylim([0,1])
            yline(0);
            yline(AcademiaPorp);
            label_h=ylabel('academia^\primes fairness');
        elseif type==3
            %ylim([-2 2])
            yline(0);
            yline(1);
            label_h=ylabel('academia^\primes score-merit slope');
        elseif type==4
            %ylim([-2 2])
            yline(1);
            label_h=ylabel('academia^\primes capital-merit slope');
        elseif type==5
            ylim([-0.1 0.8])
            yline(0);
            label_h=ylabel('overall capital-merit corr');
        end
        ylims=ylim;
        text(label_h.Position(1),ylims(2)+diff(ylims)*0.07,char(64+pt),'Fontsize',fontS)
    end
    %     subplot(2,4,1)
    %     label_h=ylabel('');
    %     ylims=ylim;
    %     text(label_h.Position(1),ylims(2)+diff(ylims)*0.07,char(65),'Fontsize',fontS)
end


predictors=[topAcadLikeRatio;topCapLikeRatio];%[authorLikeRatio;topAcadLikeRatio;topCapLikeRatio];
responses=[Merit_acad;Port_top;rankMeritSlope;capMeritSlope];
for pred=1:size(predictors,1)
    %figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/1.5]);
    figure('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/4 scrsz(4)/2.7]);
    intermediateGroup=find(ParamValues(1,:)>0 & ParamValues(2,:)>0 & ParamValues(3,:)>0  & ParamValues(4,:)>0 & ParamValues(5,:)>=0  & ParamValues(6,:)>0 & ParamValues(8,:)==0.02); %define group with 1/3 evaluation and 0.4 mate choice and get indices
    intPoint1=find(ParamValues(1,:)==0.33 & ParamValues(2,:)==0.4 & ParamValues(3,:)==0.4  & ParamValues(4,:)==0.8 & ParamValues(5,:)==25  & ParamValues(6,:)==0.4 & ParamValues(7,:)==1 & ParamValues(8,:)==0.02);
    intPoint2=find(ParamValues(1,:)==0.67 & ParamValues(2,:)==0.4 & ParamValues(3,:)==0.4  & ParamValues(4,:)==0.8 & ParamValues(5,:)==25  & ParamValues(6,:)==0.4 & ParamValues(7,:)==1 & ParamValues(8,:)==0.02);
    intPoint3=find(ParamValues(1,:)==0.0 & ParamValues(2,:)==0.8 & ParamValues(3,:)==0.8  & ParamValues(4,:)==0.8 & ParamValues(5,:)==50  & ParamValues(6,:)==0.8 & ParamValues(7,:)==1 & ParamValues(8,:)==0.02);
    for pt=1:size(responses,1)
        subplot(2,size(responses,1)/2,pt)
        %subplot(1,3,pt)
        hold on
        %scatter(predictors(pred,:),responses(pt,:),14,'k','filled');
        [x,I]=sort(predictors(pred,intermediateGroup));
        I0=find(x>0 & ~isnan(x) & imag(x)==0);
        x=x(I0);
        I=I(I0);
        trim=0; %a priori trim extreme values in predictor
        x=x(trim+1:end-trim);
        I=I(trim+1:end-trim);
        y=responses(pt,intermediateGroup);
        y=y(I);
        scatter(x,y,6,'k','filled');
        miny=min(y);
        maxy=max(y);
        [fitobj,gof,output]=fit(x',y','exp2');
        predy=feval(fitobj,x');
        inty=predint(fitobj,x',0.95,'functional','off');
        if isnan(max(diff(inty')))
            [fitobj,gof,output]=fit(x',y','exp1');
            predy=feval(fitobj,x');
            inty=predint(fitobj,x',0.95,'functional','off');
        end
        if min(diff(inty'))<eps
            [fitobj,gof,output]=fit(x',y','poly1');
            predy=feval(fitobj,x');
            inty=predint(fitobj,x',0.95,'functional','off');
        end
        [diffSize goodI]=find(diff(inty')<(maxy-miny)/2); %find x or y indices where interval size is sufficiently small
        bl1=boundedline(x(goodI)',predy(goodI),[predy(goodI)-inty(goodI,1) inty(goodI,2)-predy(goodI)],'k','alpha'); set(bl1,'linewidth',2);
        ax=gca;
        ax.XAxis.Scale='log';
        if pred==0
            xlabel 'academia^\primes author likelihood (LR)'
        elseif pred==1
            xlabel 'top academia^\primes author likelihood (LR)'
        else
            xlabel 'top capital^\primes author likelihood (LR)'
        end
        xticks([0.25 1 1.414 2 2.828 4 16 64 256])
        ax.XMinorTick='off';
        xticklabels([0.25 1 1.4 2 2.8 4 16 64 256])
%         xticks([0.25 1 2 4 8 16 32 64 256])
%         ax.XMinorTick='off';
%         xticklabels([0.25 1 2 4 8 16 32 64 256])
        scatter(predictors(pred,intPoint1),responses(pt,intPoint1),50,'db','filled')
        %kscontour([predictors(pred,intPoint1);responses(pt,intPoint1)]','color','blue','Nlevels', 3);
        scatter(predictors(pred,intPoint2),responses(pt,intPoint2),50,'dr','filled')
        %kscontour([predictors(pred,intPoint2);responses(pt,intPoint2)]','color','vermillion','Nlevels', 3);
        scatter(predictors(pred,intPoint3),responses(pt,intPoint3),50,'dg','filled')
        legend('off')
        if pt==1
            label_h=ylabel('academics^\primes mean merit');
            ylim([99,133]) %90 155
            %ylim([98,143])
            yline(100);
        elseif pt==2
            ylabel 'academia^\primes fairness'
            ylim([0,1]) %0.65
            yline(0);
            yline(AcademiaPorp);
        elseif pt==4
            ylabel 'academia^\primes capital-merit slope'
            ylim([-1 0.6])
            %ylim([-1.7 1.3])
            yline(0);
            yline(1);
        elseif pt==3
            ylabel 'academia^\primes score-merit slope'
            ylim([-0.5 1.25])
            %ylim([-1.2 1.7])
            yline(0);
            yline(1);
        end
        ylims=ylim;
        text(label_h.Position(1),ylims(2)+diff(ylims)*0.07,char(64+pt),'Fontsize',fontS)
        xlim([0.94 6])
        %xlim([0.94 8.5])
    end
end
