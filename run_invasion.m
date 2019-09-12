function SPP=run_invasion(ParamSet,SHAKE)

%SPP=run_invasion(ParamSet,SHAKE)
%
%This function runs an invasion analysis on a built-in parameter set.  There are 10 parameter sets possible, and you chose which one with the value ParamSet:
%1- The 'trade-off' set (here called tradeoff)
%2- The 'equal predation' set (here called justDisp).
%3- The 'equal dispersal' set (here called justPred).
%4- The 'yield difference' set (here called justYield).
%5- A set where species differ in dispersal, and predator are absent (the noPred set).
%6- The 'random parameters' set (here called rand).
%7- The 'random parameters 2' set (here called rand2).
%8- A third set of random parameters (called rand3).
%9- A set where species differ in their yields, and predators are absent (the yieldNoPred set).
%10- A set where species differ in sensitivity to predators and have essentially 100% dispersal (the farDisp set).
%
%Right now the parameters are hard-coded in.  You can change them by altering the code.
%
%Additionally, the variable SHAKE is set to 1 if you want to use the "shaken" condition which removes spatial structure, or 0 to run a normal situation.
%
%The code gives a junk output, but saves everything important in a file called ['20190603_invade_',name] or ['20190603_invade_shaken_',name].  The important values to care about are:
%res1- The actual invader growth rate. (a 120xSPP matrix, each value is for one species as invader during one replicate simulation)
%res2- This is the estimated growth rate, based on approximations. (a 120xSPP matrix, each value is for one species as invader during one replicate simulation)
%res3- This is the Delta mechanisms. (a 120xSPPx6 matrix, each value is for one species as invader during one replicate simulation; the third dimension is for the 6 different mechanisms, Delta Y, Delta P, Delta kappa_P, Delta kappa'_P, Delta kappa_C, and Delta kappa'_C)
%res5- This is what the actual growth rate would be without demographic stochasticity. (a 120xSPP matrix, each value is for one species as invader during one replicate simulation)
%res6- This is the pair correlation function of each species. (a 9xSPPxSPPx120 matrix.  The value res6(i,j,k,l) gives the pair correlation at distance i for species j when species k is the invader during replicate l)
%actual_stable- The mean invader growth rate across all species, measured in actual population change (here it is the natural log of lambda).
%actual_stable2- The mean invader growth rate across species, measured using expected population change (here it is the natural log of lambda).  
%guess_stable- The mean invader growth rate of the "shaken" DE approximation (here it is the natural log of lambda).
%stabilizing- A 1x6 vector.  Each value is the stabilizing effect of each mechanism (1-Delta Y, 2-Delta P, 3-Delta kappa_P, 4-Delta kappa'_P, 5-Delta kappa_C, 6-Delta kappa'_C).
%fit_diff- A 1xSPPx6 vector.  fit_diff(1,i,j) gives the mean fitness-difference effect of mechanism j on species i.
%minMax_fit_diff- A 1x6 vector.  Each value gives the range of fitness differences for each mechanism (i.e. the highest value minus the lowest value).
%boostGrowth- A 1x6 vector.  Each value is the fraction of species who have their invader growth rates boosted by that mechanism.  

rng('default');


ParamSet

%The number of species
SPP=8;

%The distribution of p and alpha values.
pmin=-.45;
pmax=-.15;
amin=6;
amax=4.5;

pred=[pmin:((pmax-pmin)/(SPP-1)):pmax]';
alpha=[amin:((amax-amin)/(SPP-1)):amax]';

%The intercept value for density-independent survival.
bint=pred*0-1;

%The amount of time before the invader is introduce.
TIME1=100;

%The amount of time the simulation runs after the invader is introduced.
TIME2=50;

%number of seeds of each species.
yield=ones(SPP,1)*7.8;

%chance that an adult dies each time step.
death=ones(SPP,1)*.4;

%This transforms alpha into the distance parameter for the simulations.
dDist=exp(alpha);

%Length and width of the community in number of adults.
LEN=200;

%The number of invaders that are put into the community.
ISTART=50;

%Here I alter the variables for the particular parameter set.
if(ParamSet==1)
    name='tradeoff'

elseif(ParamSet==2)
    name='justDisp'
    pred(:)=mean(pred);

elseif(ParamSet==3)
    name='justPred'
    dDist(:)=exp(mean(alpha))

elseif(ParamSet==4)
    name='justYield'
    yield=[7.5:(8.1-7.5)/(SPP-1):8.1]';
    pred(:)=mean(pred)*1.5;
    dDist(:)=exp(mean(alpha))

elseif(ParamSet==5)
    name='noPred'
    pred(:)=10^-6;
    dDist(:)=exp(alpha);

elseif(ParamSet==6)
    name='rand'
    yield=[7.5:(8.1-7.5)/(SPP-1):8.1]';
    pred(:)=pred([8     6     4     5     1     7     3     2]);
    dDist=exp(alpha([6     1     5     3     7     2     8     4]));

elseif(ParamSet==7)
    name='rand2'
    yield=[7.5:(8.1-7.5)/(SPP-1):8.1]';
    pred(:)=pred([6     3     7     8     5     1     2     4]);
    dDist=exp(alpha([8     3     6     7     5     1     2     4]));

elseif(ParamSet==8)
    name='rand3'
    yield=[7.5:(8.1-7.5)/(SPP-1):8.1]';
    pred(:)=pred([5     2     1     7     6     8     3     4]);
    dDist=exp(alpha([7     6     3     8     5     4     1     2]));

elseif(ParamSet==9)
    name='yieldNoPred'
    yield=[7.5:(8.1-7.5)/(SPP-1):8.1]';
    pred(:)=pred*0;
    dDist(:)=exp(mean(alpha));

elseif(ParamSet==10)
    name='farDisp'
    dDist(:)=exp(mean(alpha*10));
end

%I set it so that the potential NDD from distance- and density-responsive predators are the same.
predSeedling=pred./yield;

%Here I run a long simulation to determine who is actually persisting in the community.
[recordX] =...
    JC_invade1(pred,predSeedling,bint,...
    yield,dDist,death,2000,0,LEN,ISTART,SHAKE,0);

%Here I figure out which species persisted to time 2000.
nbar=recordX(2000,:)
alive=(nbar>0).*[1:SPP];
alive(alive==0)=[]

%Here I calculate the equilibrium frequencies of each species.
nbar(nbar==0)=[]

%Here I remove all of the species who did not persist.
SPP=length(alive)
pred=pred(alive)
predSeedling=predSeedling(alive)
bint=bint(alive)
yield=yield(alive);
death=death(alive);
dDist=dDist(alive);


%REPS is the number of replicate simulations I run.
REPS=120;

res1=zeros(REPS,SPP);
res2=zeros(REPS,SPP);
res3=zeros(REPS,SPP,6);
res5=res1;
res6=zeros(9,SPP,SPP,REPS);


%I let each invader grow for SKIP time steps before I begin tracking it.
SKIP=ceil(5/death(1));

%This is the actual time step where I start tracking the invader.
STARTHERE=TIME1+ceil(5/death(1));

TIME=TIME1+TIME2;

for INVADE=1:SPP

    %Here I remove the invader from the equilibrium community, and make every other species' frequency increase proportionately.
    nbar2=nbar;
    nbar2(INVADE)=0;
    nbar2=nbar2/sum(nbar2);
    
    nstart=nbar2*triu(ones(SPP));
    
    vals=ceil(LEN^2*nstart);
    
    %Input is the initial distribution of the species.  It will be shuffled later.
    input=ones(1,LEN.^2);
    
    for(ii=2:SPP)
        if(INVADE~=ii)
            input((1+vals(ii-1)):LEN^2)=ii;
        end
    end
        
    %Now I run the invasion analysis a bunch of times and record the results.
    parfor rep=1:REPS
        [rep INVADE]
        
        [record, record_guess,mechs, egrow,autoCor] =...
            JC_invade1(pred,predSeedling,bint,...
            yield,dDist,death,TIME1,TIME2,shuffle(input),ISTART,SHAKE,INVADE);
        
        ri=record(:,INVADE);
        ri(ri==0)=[];
        
        res1(rep,INVADE)=(log(ri(length(ri)))-log(ri(SKIP)))/(length(ri)-SKIP);
        res5(rep,INVADE)=mean(egrow(INVADE,(TIME1+SKIP):(TIME1+length(ri)-2)),2);
        
        
        ri2=log(record_guess(:,INVADE));
        ri2(isinf(ri2))=[];
        
        res2(rep,INVADE)=(ri2(length(ri2))-ri2(SKIP))/(length(ri2)-SKIP);
        
        here=(record(1:(TIME-1),INVADE)>0);
        here(TIME1:TIME1+SKIP)=0;
        RES=[1:SPP];
        RES(INVADE)=[];
        
        for ix=1:6
            if(ix==1)
                zz=1;
            else
                zz=-1;
            end
            res3(rep,INVADE,ix)=zz*...
                mean(mechs(INVADE,here,ix))-mean(mean(mechs(RES,here,ix)));
        end
        
        res6(:,:,INVADE,rep)=autoCor;
        
    end
    mean(res1)/.2
    mean(res2)/.2
    mean(res3)
end

actual_stable=mean(res1)
actual_stable2=mean(res2)
guess_stable=mean(res5)
stabilizing=reshape(mean(mean(res3)),1,6)
fit_diff=mean(res3)-mean(mean(res3))
minMax_fit_diff=reshape((max(mean(res3))-min(mean(res3))),1,6)
boostGrowth=reshape(mean(((mean(res3))<0)),1,6)



if(SHAKE)
    save(['20190603_invade_shaken_',name])
else
    save(['20190603_invade_',name])
end


plot(recordX)