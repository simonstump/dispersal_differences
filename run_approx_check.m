function SPP=run_approx_check(ParamSet)

%SPP=run_approx_check(ParamSet,SHAKE)
%
%This function runs a "shaken" simulation, along with the DE approximation.  It is used for checking how good the approximation is.  There are 10 parameter sets possible, and you chose which one with the value ParamSet:
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
%The code gives a junk output, but saves everything important in a file called ['20190614_approx_check_',name].  The important values to care about are:
%record- This is a TIMExSPP matrix.  It tracks the frequency of each species.  To make a graph of the population dynamics, you can simply type plot(nbarrec).
%rg- This is a TIMExSPP matrix.  While the main simulation is running, I also run a DE version of the shaken approximation in the background (whether or not SHAKE=1; this does not affect the simulation).  rg is the frequency of each species in the DE model. 
%mechs- This is a SPPx(TIME-1)x6 matrix.  The value mechs(i,j,x) gives the value for calculating Delta mechanism x for species i at time j.  To calculate the actual Delta terms, take the value for the invader and subtract the average value of the resident.  mechs(i,j,1) gives ln(YI) (for Delta Y).  mech(i,j,2) gives N_j*mean(p) (for Delta P).  mechs(i,j,3) gives (1-n)*(p)*(sum(d*dP)) (for Delta kappa_P.  mechs(i,j,4) gives the covariance used to make Delta kappa'_P.  mechs(i,j,5) gives the approximation for E[cov(d,ln{C})], eqn A55 (Delta kappa_C).  mechs(i,j,6) gives the covariance used for Delta kappa'_C.
%egrow- This is a SPPx(TIME-1) matrix.  It gives the expected fitness of each species at each time step (i.e. lambda-tilde).  
%ripK- This is a 9xSPP matrix.  The value ripK(i,j) gives the pair-correlation value of species i at distance j.
%covP- This is a SPPx(TIME-1) matrix.  It gives the actual covariance between the number of seeds of sp. j arriving at a site and the predation risk at that site for each species at each time step.
%covC- This is a SPPx(TIME-1) matrix.  It gives the actual covariance between the number of seeds of sp. j arriving at a site and the total number of seedlings at that site for each species at each time step.
%meanP- This is a (TIME-1)xSPP matrix.  This is the average risk of predation across all sites for each species.
%meanC- This is a (TIME-1)x1 vector.  It is the average number of seedlings (post-predation) at each site.


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

%The amount of time before the invader is introduce. (if not an invasion analysis, the simulation is run for TIME1+TIME2)
TIME1=100;

%The amount of time the simulation runs after the invader is introduced. (if not an invasion analysis, the simulation is run for TIME1+TIME2)
TIME2=200;

%number of seeds of each species.
yield=ones(SPP,1)*7.8;

%chance that an adult dies each time step.
death=ones(SPP,1)*.4;

%This transforms alpha into the distance parameter for the simulations.
dDist=exp(alpha);

%Length and width of the community in number of adults.
LEN=300;

%The number of invaders that are put into the community (not currently used)
ISTART=50;
SHAKE=1;

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

predSeedling=pred./yield;

%set INVADE=1 if you want to run the approximation checker as an invasion analysis.
INVADE=0;

[record, rg, mechs, egrow,ripK,covP,covC,meanP,meanC] =...
    JC_invade1(pred,predSeedling,bint,...
    yield,dDist,death,TIME1,TIME2,LEN,0,SHAKE,INVADE);

close all
save(['20190614_approx_check_',name])

