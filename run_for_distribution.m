%This code generates the data for fig. 2b and 2c (i.e. the distribution of individuals when their density is held at a particular value).  It saves the data in a file called '20190614holdDist'.  The distributions are saved as x1 (2b left), x2 (2b right), x3 (2c left), and x4 (2c right).  


clc, clear

rng('default');

%Number of species
SPP=8;

%The distribution of p and alpha values.
pmin=-.65;
pmax=-.05;
amin=7;
amax=4;

pred=[pmin:((pmax-pmin)/(SPP-1)):pmax]';
alpha=[amin:((amax-amin)/(SPP-1)):amax]';

%The intercept value for density-independent survival.
bint=pred*0-1;

%TIME1+TIME2 is the length of the simulation.
TIME1=150;
TIME2=150;

%number of seeds of each species.
yield=ones(SPP,1)*7.8;

%chance that an adult dies each time step.
death=ones(SPP,1)*.4;

%This transforms alpha into the distance parameter for the simulations.
dDist=exp(alpha);

%Length and width of the community in number of adults.
LEN=70;

%Leave as 0 to keep spatial structure.
SHAKE=0;

%I set it so that the potential NDD from distance- and density-responsive predators are the same.
predSeedling=pred./yield;

[recordX] =...
    JC_invade1(pred,predSeedling,bint,...
    yield,dDist,death,2000,0,LEN,50,0,0);


nbar=recordX(2000,:)
alive=(nbar>0).*[1:SPP];
alive(alive==0)=[]
nbar(nbar==0)=[]
SPP=length(alive)
pred=pred(alive)
predSeedling=predSeedling(alive)
bint=bint(alive)
yield=yield(alive);
death=death(alive);
dDist=dDist(alive);


%nval are the frequencies of focal species
nval=[0.03, .15];


i=1;

%INVADE1 is the species in fig. 2b
INVADE1=1

%INVADE2 is the species in fig. 2c
INVADE2=length(pred)

[K1, a, b, x1] = JC_holdAbund2(pred,predSeedling,bint,...
    yield,dDist,death,TIME1,TIME2,LEN,round(LEN^2*nval(i)),SHAKE,INVADE1,[300]);

i=2;

[K2, a, b, x2] = JC_holdAbund2(pred,predSeedling,bint,...
    yield,dDist,death,TIME1,TIME2,LEN,round(LEN^2*nval(i)),SHAKE,INVADE1,[300]);
i=1;
INVADE=length(pred)

[K3, a, b, x3] = JC_holdAbund2(pred,predSeedling,bint,...
    yield,dDist,death,TIME1,TIME2,LEN,round(LEN^2*nval(i)),SHAKE,INVADE2,[300]);
i=2;
INVADE=length(pred)

[K4 a, b, x4] = JC_holdAbund2(pred,predSeedling,bint,...
    yield,dDist,death,TIME1,TIME2,LEN,round(LEN^2*nval(i)),SHAKE,INVADE2,[300]);


save('20190614holdDist')