%This is the code used to generate fig. A6.  It sets up 15 values for alpha (the seed dispersal parameter) and p_j^A (the potential CNDD for distance-responsive predators).  It then creates random communities by sampling the p_j^A values (without replacement) to make communities of 15 species (the alpha parameters are listed from highest to lowest).  It lets those species compete for 1500 time steps, and then sees which species have persisted.  It repeats this 150 time.
%
%It saves the critical data in a file called 'test_emergentTradeoff'.  The important variables are:
% permPred- This is a 150x15 matrix. Each row contains the numbers 1 through 15.  predPerm(i,j) indicates that in replicate community i, species number j has that value of p_j^A.  Thus, the vector of p_j^A and p_j^S for simulation i are pred(permPred(i,:)) and predSeedling(permPred(i,:)).
% isAlive- This is a 150x15 matrix.  The value isAlive(i,j) will be 1 if species j persisted during simulation i, and 0 otherwise.
% startCor- This is a 150x1 vector.  Value i is the correlation between alpha and p_j^A in the initial 15 species community.
% endCor- This is a 150x1 vector.  Value i is the correlation between alpha and p_j^A amongst the species who persisted to the end of the simulation.

clc, clear

rng('default');

COMS=150;  %This is the number of communities that are constructed.

SPP=15;   %This is the number of species in the initial community.

%The predator susceptibility values are uniformly distributed between pmin and pmax.
pmin=-.5;
pmax=-.1;

%The alpha values used to calculate dispersal are uniformly distributed between amin and amax.
amin=7;
amax=5;

pred=[pmin:((pmax-pmin)/(SPP-1)):pmax]';
alpha=[amin:((amax-amin)/(SPP-1)):amax]';

%Here we set the impact of density-responsive predators; to get their true value, this must be multiplied by the yield.
predSeedling=pred/10;

%bint sets the baseline survival (on the log-odds scale).
bint=pred*0-1;

%TIME1 is how long the community is simulated for.
TIME1=1500;
TIME2=0;  %leave this as 0.

%This is the number of seeds each species produces.
yield=ones(SPP,1)*7.8;

%This is the adult death rate.
death=ones(SPP,1)*.4;

%This is used to calculate dispersal distances.
dDist=exp(alpha);

%We assumed a tree is a 10m square
TREESIZE=10;

%This is the length and width of the community in number of trees.
LEN=120;

%permPred is used to randomly assign the values for p_j (the impact of predators).  permPred(i,j) is p_j value for species j during simulation i.
permPred=zeros(COMS,SPP);

for i=1:COMS
    permPred(i,:)=randperm(SPP);
end

%This is what would happen if there is a perfect trade-off.
[record, record_guess, mechs, egrow,ripK,covP,covC,meanP,meanC] =...
    JC_invade1(pred,predSeedling,bint,...
    yield,dDist,death,TIME1,TIME2,LEN,0,0,0);
    endN(i,:)=(record(TIME1+TIME2,:))


baselineN=(record(TIME1+TIME2,:))

%return

%endN(i,j) is the frequency of species j at the end of simulation i.
endN=permPred*0;


%Here I run all of the simulations.
parfor i=1:COMS
i
[record, record_guess, mechs, egrow,ripK,covP,covC,meanP,meanC] =...
    JC_invade1(pred(permPred(i,:)),predSeedling(permPred(i,:)),bint,...
    yield,dDist,death,TIME1,TIME2,LEN,0,0,0);
    endN(i,:)=(record(TIME1+TIME2,:));
end

%isAlive(i,j) is 1 if species j has non-zero density at the end of simulation i, and 0 otherwise.
isAlive=(endN>0);

%startCor is the correlation between p_j and alpha in each initial community.
startCor=zeros(COMS,1);

%end is the correlation between p_j and alpha amongst persisting species at the end of each simulation.
endCor=zeros(COMS,1);

for i=1:COMS
    temp=pred(permPred(i,:));
    startCor(i)=corr(alpha, temp);
    endCor(i)=corr(alpha(isAlive(i,:)),temp(isAlive(i,:)));
end
startCor
endCor

plot(startCor,endCor,'.',[-1 1],[-1 1],'k--')

save('test_emergentTradeoff')
