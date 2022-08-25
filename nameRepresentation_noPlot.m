function [authorLikeRatio,topAcadLikeRatio,botAcadLikeRatio,topCapLikeRatio,botCapLikeRatio,Port_top,rankMeritSlope,capMeritSlope,capMeritCor,LRMeritSlope,Merit_pop,Merit_top,Merit_acad,Capital_pop,Capital_top,Capital_acad]=nameRepresentation_noPlot(numGen,numNames,initPopPerName,nameDistr,initAcademiaPop,reprodRate,reprodCost1,reprodCost2,evaluation,merit_mean,merit_var,capital_mean,capital_var,MeritCapCausal,mateChoice,merit_herit,capital_herit,capital_added,capital_added_var,child_surname,child_surnameRelVar,nameMutation)

%nameRepresentation.m
%E.W. Tekwa Apr 28, 2022
%simulate name-publication associations

fs=11; %font size
% set(0,'DefaultAxesFontSize',fs)
% scrsz = get(0,'ScreenSize');
% figure('Color', [1 1 1]);

%clear
% rng(1); %set random number generator seed

%Create general population with sexual reproduction, heritable ability, and
%heritable academic capital. This population feeds into a subpopulation of
%academia, which can become reproductively isolated to a degree through
%mate choice

% numGen=20; %number of generations to simulate
% numNames=2000; %initial number of names
% initPopPerName=50; %per name mean population
% nameDistr=2; %name frequency distribution is Poisson (1) or Exponential (2)
% AcademiaPorp=0.01; %percent of people in academia, 0.0022 global
% initAcademiaPop=round(numNames*initPopPerName*AcademiaPorp); %number of people in academia
% reprodRate=2; %mean children per pair of parents per generation, actual is Poisson random
% reprodCost1=0.1; %reduction in numbeer of children per generation for couples with one academics
% reprodCost2=0.2; %reduction in numbeer of children per generation for couples with two academics
% evaluation=2; %evaluation for academic entrance based on merit (0), capital (1), or capital+merit(2)
% merit_mean=100;
% merit_var=(merit_mean/2)*(merit_mean>0);
% capital_mean=100;
% capital_var=(capital_mean/2)*(capital_mean>0);
% mateChoice=1; %mate choice random (0) or ranked by capital (1) or ranked by merit (2)
% mateChoiceRelVar=0.1*(mateChoice>0); %mate choice ranking variability scaled to score of 100 (this variability is added to true capital or capital+merit scores before ranking)
% merit_herit=0.5; %0 to 1
% capital_herit=0.9;
% capital_added=0; %mean capital added to children of academics
% capital_added_var=(capital_added/2)*(capital_added>0); %variance in captial added to children of academics
% child_surname=1; %Child inherits surname of parent with higher capital (1), or with higher capital+merit (2), or random (0)
% child_surnameRelVar=0.1*(child_surname>0); %name convention ranking variability scaled to score of 100 (this variability is added to true capital or capital+merit scores before ranking)
% nameMutation=0.002; %probability of surname mutation per child per generation

%initialize population variables
PopPerNameEmpty=60; %60*numGen/20; %multiplicative extra room over initial population per name to store population growth
NameEmpty=numGen; %5*numGen/20; %multiplicative extra room over initial names to store name growth
GenAuthors=zeros(numNames*NameEmpty,numGen); %authors per surname over generations
GenPop=zeros(numNames*NameEmpty,numGen); %people per surname over generations
GenIndiv_ID=repmat([1:numNames*NameEmpty]',1,PopPerNameEmpty*initPopPerName); %assign surname ID to individuals (20x for potential growth of individuals in a surname)
GenIndiv_merit=nan(numNames*NameEmpty,PopPerNameEmpty*initPopPerName,numGen); %merit of each individual with a surname over generations
GenIndiv_capital=nan(numNames*NameEmpty,PopPerNameEmpty*initPopPerName,numGen); %capital of each individual with a surname over generations
AcademicIndiv=nan(initAcademiaPop,numGen); %individuals
AcademicOutput=nan(1,numGen); %academic merit per capita
AcademicCapital=nan(1,numGen); %academic capital per capita
GenMerit=zeros(1,numGen); %merit per capita
GenCapital=zeros(1,numGen); %capital per capita
GenTopMeritInAcad=zeros(1,numGen); %portion of people with merits ranked within the academia size who are in academia
GenMerit_TopMeritInPop=zeros(1,numGen); %mean merit of top-merit people in the general population
GenCapital_TopMeritInPop=zeros(1,numGen); %mean capital of top-merit people in the general population
rankMeritSlopes=zeros(1,numGen); %regression slope of merit on academic score within academia
capMeritSlopes=zeros(1,numGen); %regression slope of merit on capital within academia
capMeritCors=zeros(1,numGen); %Pearson correlation of merit and capital in population
numNames_new=numNames; %current number of names including extinct and existing (with new mutations)

%initial population values
if nameDistr==1
    GenIndiv=poissrnd(initPopPerName,numNames,1);
else
    GenIndiv=round(random('exponential',initPopPerName,numNames,1));
end
GenPop(1:numNames,1)=GenIndiv;
for name=1:numNames %only assign names to initiated individuals among initially existing surnames
    GenIndiv_merit(name,1:GenIndiv(name),1)=merit_mean+sqrt(merit_var)*randn(1,GenIndiv(name));
    %GenIndiv_capital(name,1:GenIndiv(name),1)=capital_mean+sqrt(capital_var)*randn(1,GenIndiv(name));
    GenIndiv_capital(name,1:GenIndiv(name),1)=(1-MeritCapCausal)*capital_mean+MeritCapCausal*(GenIndiv_merit(name,1:GenIndiv(name),1))+sqrt(capital_var)*randn(1,GenIndiv(name)); %additional causal relationship between capital and merit
end

%run dynamics
for t=1:numGen
    %academic selection:
    GenIndiv_merit_t=GenIndiv_merit(:,:,t);
    GenIndiv_capital_t=GenIndiv_capital(:,:,t);
    GenIndiv_total_t=GenIndiv_merit_t+GenIndiv_capital_t;
    GenIndiv_score_t=(1-evaluation)*GenIndiv_merit_t+evaluation*GenIndiv_capital_t; %weighted evaluation score
    GenPop(:,t)=sum(GenIndiv_total_t>0,2); %track population size, only count those with total scores above 0 (otherwise death)
    GenMerit(:,t)=nanmean(GenIndiv_merit_t,'all');
    GenCapital(:,t)=nanmean(GenIndiv_capital_t,'all');
%    if evaluation==1 %academic entrance based on capital
    [rank, rankID]=sort(GenIndiv_score_t(:),'descend','MissingPlacement','last'); %sort all individuals by evaluation score
%     elseif evaluation==2 %academic entrance based on capital+merit
%         [rank, rankID]=sort(GenIndiv_total_t(:),'descend','MissingPlacement','last'); %sort all individuals by merit+capital
%     else %academic entrance based on merit
%         [rank, rankID]=sort(GenIndiv_merit_t(:),'descend','MissingPlacement','last'); %sort all individuals by merit
%     end
    mdl_merit_score=fitlm(rank(1:initAcademiaPop),GenIndiv_merit_t(rankID(1:initAcademiaPop))); %perform linear regression of merit (response) on academic score (predictor) within academia
    rankMeritSlopes(t)=mdl_merit_score.Coefficients.Estimate(2); %record regression slope within academia
    mdl_merit_cap=fitlm(GenIndiv_capital_t(rankID(1:initAcademiaPop)),GenIndiv_merit_t(rankID(1:initAcademiaPop))); %perform linear regression of merit (response) on capital (predictor) within academia
    capMeritSlopes(t)=mdl_merit_cap.Coefficients.Estimate(2); %record regression slope within academia
    capMeritCors(t)=corr(GenIndiv_merit_t(:),GenIndiv_capital_t(:),'Rows','complete'); %record correlation between merit and capital in population
    AcademicOutput(t)=mean(GenIndiv_merit_t(rankID(1:initAcademiaPop))); %record mean academic merit
    AcademicCapital(t)=mean(GenIndiv_capital_t(rankID(1:initAcademiaPop))); %record mean academic capital
    AcademicIndiv(:,t)=GenIndiv_ID(rankID(1:initAcademiaPop)); %record surnames of academic individuals
    GenIndiv_capital_t(AcademicIndiv(:,t))=GenIndiv_capital_t(AcademicIndiv(:,t))+capital_added+sqrt(capital_added_var)*randn(initAcademiaPop,1); %add capital to academics
    %GenIndiv_merit(:,:,t)=GenIndiv_merit_t; %replace current generation merit
    for i=1:initAcademiaPop
        GenAuthors(AcademicIndiv(i,t),t)=GenAuthors(AcademicIndiv(i,t),t)+1; %add authorship record by rank
    end
    [rankMerit, rankIDMerit]=sort(GenIndiv_merit_t(:),'descend','MissingPlacement','last'); %sort all individuals by merit
    GenTopMeritInAcad(t)=numel(intersect(rankIDMerit(1:initAcademiaPop),rankID(1:initAcademiaPop)))/initAcademiaPop; %portion of top merit people in academia
    GenMerit_TopMeritInPop(t)=mean(rankMerit(1:initAcademiaPop)); %mean merit of top merit people in population
    GenCapital_TopMeritInPop(t)=mean(GenIndiv_capital_t(rankIDMerit(1:initAcademiaPop))); %mean capital of top merit people in population
    
    %reproduction and inheritance
    if t<numGen
        IndivPos=find(GenIndiv_total_t>0); %get current generation's individual positions
        %pick order of mating parents from current generation (each consecutive
        %pair eg. 1&2 3&4 try for children):
        GenIndiv_capital_t_mate=mateChoice*GenIndiv_capital_t+(1-mateChoice)*capital_mean+sqrt(capital_var)*randn; %mate assortment is the sum of true capital (weighted by mateChoice), population mean capital (weighted by 1-mateChoice), and noise
%         if mateChoice==1 %pair by capital:
%             GenIndiv_capital_t_mate=GenIndiv_capital_t+sqrt(mateChoiceRelVar*capital_mean)*randn(size(GenIndiv_capital_t)); %noise added to capital data before ranking
            [rank_c, rankID_c]=sort(GenIndiv_capital_t_mate(:),'descend','MissingPlacement','last'); %sort all individuals by noisy capital
            ParentOrder=rankID_c;
%         elseif mateChoice==2 %pair by capital+merit:
%             GenIndiv_capital_t_mate=GenIndiv_total_t+sqrt(mateChoiceRelVar*(capital_mean+merit_mean))*randn(size(GenIndiv_capital_t)); %noise added to capital+merit data before ranking
%             [rank_c, rankID_c]=sort(GenIndiv_capital_t_mate(:),'descend','MissingPlacement','last'); %sort all individuals by noisy capital+merit
%             ParentOrder=rankID_c;
%         else %random
%             ParentOrder=randsample(IndivPos(:),sum(GenPop(:,t))); %order individuals randomly
%         end
        %pick children's surname from parents:
        if child_surname==1 %surname inherited by capital:
            GenIndiv_capital_t_name=GenIndiv_capital_t+sqrt(child_surnameRelVar*capital_mean)*randn(size(GenIndiv_capital_t)); %noise added to capital data before ranking
        elseif child_surname==2 %surname inherited by capital+merit:
            GenIndiv_capital_t_name=GenIndiv_capital_t+sqrt(child_surnameRelVar*(capital_mean+merit_mean))*randn(size(GenIndiv_capital_t)); %noise added to capital+merit data before ranking
        end
        for i=1:floor(sum(GenPop(:,t),1)/2) %for each living pair
            if GenIndiv_capital_t(ParentOrder(2*i-1))>0 && GenIndiv_capital_t(ParentOrder(2*i))>0 %check parents not dead
                numAcadInPair=sum(rankID(1:initAcademiaPop)==ParentOrder(2*i-1) + rankID(1:initAcademiaPop)==ParentOrder(2*i)); %number of parents in pair who is an academic
                if numAcadInPair==1
                    children=poissrnd(reprodRate-reprodCost1); %number of children born to two academics
                elseif numAcadInPair==2
                    children=poissrnd(reprodRate-reprodCost2); %number of children born to one academics
                else
                    children=poissrnd(reprodRate); %number of children born to non-academics
                end
                if children>0
                    surname1=GenIndiv_ID(ParentOrder(2*i-1)); %get surname of parent1, which children inherit if egalitarian
                    surname2=GenIndiv_ID(ParentOrder(2*i)); %get surname of parent2
                    if child_surname==1 && GenIndiv_capital_t_name(ParentOrder(2*i-1))<GenIndiv_capital_t_name(ParentOrder(2*i)) %children inherit parent surname with highest capital
                        temp=surname1;
                        surname1=surname2; %swap parents order
                        surname2=temp; %swap parents order
                    elseif child_surname==2 && GenIndiv_total_t_name(ParentOrder(2*i-1))<GenIndiv_total_t_name(ParentOrder(2*i)) %children inherit parent surname with highest capital+merit
                        temp=surname1;
                        surname1=surname2; %swap parents order
                        surname2=temp; %swap parents order
                    end
                    for child=1:children
                        if rand<nameMutation %mutate to new child surname (independent for each child)
                            surname_child=numNames_new+1;
                            numNames_new=numNames_new+1;
                        else
                            surname_child=surname1; %else keep parent's name
                        end
                        surname_cur=sum(~isnan(GenIndiv_merit(surname_child,:,t+1))); %get number of people already with the children's surname in next generation
                        %add children's merits to next generation:
                        %GenIndiv_merit(surname_child,surname_cur+1,t+1)=(1-merit_herit)*(merit_mean+sqrt(merit_var)*randn)+merit_herit*mean(GenIndiv_merit_t(ParentOrder(2*i-1:2*i)));
                        GenIndiv_merit(surname_child,surname_cur+1,t+1)=merit_herit*mean(GenIndiv_merit_t(ParentOrder(2*i-1:2*i)))+(1-merit_herit)*(merit_mean)+sqrt(merit_var)*randn;
                        %add children's inherited capitals to next generation:
                        %GenIndiv_capital(surname_child,surname_cur+1,t+1)=(1-capital_herit)*(capital_mean+sqrt(capital_var)*randn)+capital_herit*mean(GenIndiv_capital_t(ParentOrder(2*i-1:2*i)));
                        %GenIndiv_capital(surname_child,surname_cur+1,t+1)=capital_herit*mean(GenIndiv_capital_t(ParentOrder(2*i-1:2*i)))+(1-capital_herit)*capital_mean+sqrt(capital_var)*randn; %no additional causal relationship between capital and merit
                        GenIndiv_capital(surname_child,surname_cur+1,t+1)=capital_herit*mean(GenIndiv_capital_t(ParentOrder(2*i-1:2*i)))+(1-capital_herit)*((1-MeritCapCausal)*capital_mean+MeritCapCausal*GenIndiv_merit(surname_child,surname_cur+1,t+1))+sqrt(capital_var)*randn; %additional causal relationship between capital and merit
                    end
                end
            end
        end
    end
end

%AuthorsPerPerson=sum(GenAuthors,2)./sum(GenPop,2); %each surname's
%probability of authorship per person measured across all generations
AuthorsPerPerson=sum(max(GenAuthors-1,0),2)./sum(max(GenPop-1,1),2); %bias-corrected
populationTrend=sum(GenPop)

LastGenPop_AuthorsPerPerson=AuthorsPerPerson(GenIndiv_ID(GenIndiv_capital(:,:,end)>0)); %get author/person for all individuals in last generation
LastGenAcademics_AuthorsPerPerson=AuthorsPerPerson(AcademicIndiv(:,end)); %get author/person for all individuals in academia in last generation
LastGenTopAcad_AuthorsPerPerson=AuthorsPerPerson(AcademicIndiv(1:initAcademiaPop/2,end)); %get author/person for all individuals in academia in last generation
LastGenBotAcad_AuthorsPerPerson=AuthorsPerPerson(AcademicIndiv(initAcademiaPop/2+1:end,end)); %get author/person for all individuals in academia in last generation
[rankCapital, rankIDCapital]=sort(GenIndiv_capital_t(:),'descend','MissingPlacement','last'); %sort all individuals by capital
lastID=find(isnan(rankCapital),1)-1;
LastGenTopCap_AuthorsPerPerson=AuthorsPerPerson(GenIndiv_ID(rankIDCapital(1:initAcademiaPop)));
LastGenBotCap_AuthorsPerPerson=AuthorsPerPerson(GenIndiv_ID(rankIDCapital(lastID-initAcademiaPop+1:lastID)));

%record variables to be returned by function:
authorLikeRatio=mean(LastGenAcademics_AuthorsPerPerson)/mean(LastGenPop_AuthorsPerPerson); %LR for all academics
topAcadLikeRatio=mean(LastGenTopAcad_AuthorsPerPerson)/mean(LastGenPop_AuthorsPerPerson); %LR for academics ranked in top half
botAcadLikeRatio=mean(LastGenBotAcad_AuthorsPerPerson)/mean(LastGenPop_AuthorsPerPerson); %LR for academics ranked in bottom half
topCapLikeRatio=mean(LastGenTopCap_AuthorsPerPerson)/mean(LastGenPop_AuthorsPerPerson); %LR for pop with top 1% capital
botCapLikeRatio=mean(LastGenBotCap_AuthorsPerPerson)/mean(LastGenPop_AuthorsPerPerson); %LR for pop with bottom 1% capital
Port_top=GenTopMeritInAcad(end);
rankMeritSlope=rankMeritSlopes(end);
capMeritSlope=capMeritSlopes(end);
capMeritCor=capMeritCors(end);
Merit_pop=GenMerit(end);
Merit_top=GenMerit_TopMeritInPop(end);
Merit_acad=AcademicOutput(end);
Capital_pop=GenCapital(end);
Capital_top=GenCapital_TopMeritInPop(end);
Capital_acad=AcademicCapital(end);

mdl_merit_authorLR=fitlm(LastGenAcademics_AuthorsPerPerson/mean(LastGenPop_AuthorsPerPerson),GenIndiv_merit_t(rankID(1:initAcademiaPop))); %perform linear regression of merit on author LR within academia
LRMeritSlope=mdl_merit_authorLR.Coefficients.Estimate(2); %record regression slope
    

% text(0.15*.2,curYlim(2)*0.6,{['entrance LR in academia (top, bottom)=' num2str(authorLikeRatio,3) ' (' num2str(topAcadLikeRatio,3) ', ' num2str(botAcadLikeRatio,3) ')'];
%     ['entrance LR in (top, bottom) capital groups=(' num2str(topCapLikeRatio,3) ', '  num2str(botCapLikeRatio,3) ')'];
%     ['academic fairness=' num2str(Port_top,2)];
%     ['score-merit slope=' num2str(rankMeritSlope,2) ', capital-merit slope=' num2str(capMeritSlope,2)];
%     ['mean merit in: academia=' num2str(Merit_acad,3) ', top merit=' num2str(Merit_top,3) ', all=' num2str(Merit_pop,3)];
%     ['mean capital in: academia=' num2str(Capital_acad,3) ', top merit=' num2str(Capital_top,3) ', all=' num2str(Capital_pop,3)];
%     ['end # names=' num2str(sum(GenPop(:,end)>0)) ', end population=' num2str(sum(GenPop(:,end)))]
%     },'FontSize',fs);
% 
% GenMerit
% AcademicOutput
% GenCapital
% AcademicCapital