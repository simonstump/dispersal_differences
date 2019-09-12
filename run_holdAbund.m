function SPP=run_holdAbund(ParamSet)

%function SPP=run_holdAbund(ParamSet)
%
%This function runs the community dynamics, holding the focal species at a set value (currently between 0.01 and 0.22).  There are 11 parameter sets possible, and you chose which one with the value ParamSet:
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
%11- This is the data used for Fig. 2a.
%
%Right now the parameters are hard-coded in.  You can change them by altering the code.
%
%Additionally, the variable SHAKE is set to 1 if you want to use the "shaken" condition which removes spatial structure, or 0 to run a normal situation.
%
%The code gives a junk output, but saves everything important in a file called ['20190612_multiHold_',name].  The important values to care about are:
%res1- This is the pair correlation function of each species. (a 9xSPPxSPPx120 matrix.  The value res6(i,j,k,l) gives the pair correlation at distance i for species j when species k is the invader during replicate l)
%res2- This is the estimated growth rate, based on approximations. (a 120xSPP matrix, each value is for one species as invader during one replicate simulation)
%res3- This is the Delta mechanisms. (a 50xSPPx6 matrix, each value is for one species as invader during one replicate simulation; the third dimension is for the 6 different mechanisms, Delta Y, Delta P, Delta kappa_P, Delta kappa'_P, Delta kappa_C, and Delta kappa'_C)

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
elseif(ParamSet==11)
    name='forFig'
    pmin=-.65;
    pmax=-.05;
    amin=7;
    amax=4;
    
    pred=[pmin:((pmax-pmin)/(SPP-1)):pmax]';
    alpha=[amin:((amax-amin)/(SPP-1)):amax]';
    
end

%I set it so that the potential NDD from distance- and density-responsive predators are the same.
predSeedling=pred./yield;

%Here I run a long simulation to determine who is actually persisting in the community.
[recordX] =...
    JC_invade1(pred,predSeedling,bint,...
    yield,dDist,death,2000,0,LEN,0,0,0);

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
REPS=50;


%I let each invader grow for SKIP time steps before I begin tracking it.
SKIP=ceil(5/death(1));

%This is the actual time step where I start tracking the invader.
STARTHERE=TIME1+ceil(5/death(1));

TIME=TIME1+TIME2;

%nval is the frequencies that I hold the focal species at.
nval=[.01:.01:.22];

res1=zeros(9,length(nval),SPP,REPS);
res2=zeros(length(nval),SPP,REPS);
res3=zeros(6,length(nval),SPP,SPP,REPS);


RUNS=length(nval)



for rep=1:REPS
    for INVADE=1:SPP
        [INVADE, rep]
        %parfor i=1:RUNS
            for i=1:RUNS
            
            [ripK,mechs,egrow] = ...
                JC_holdAbund2(pred,predSeedling,bint,...
                yield,dDist,death,TIME1,TIME2,LEN,round(LEN^2*nval(i)),0,INVADE,0);
            res1(:,i,INVADE,rep)=ripK(:,INVADE);
            res2(i,INVADE,rep)=mean(egrow(INVADE,:));
            
            
            for iii=1:SPP
                res3(:,i,INVADE,iii,rep)=mean(mechs(iii,TIME1:(TIME-1),:),2);
            end
        
        end
    end
    save(['20190612_multiHold_',name])
end


save(['20190612_multiHold_',name])
