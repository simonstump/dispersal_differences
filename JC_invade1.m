function [nbarrec, rg, mechs, egrow,ripK,covP,covC,meanP,meanC] = JC_invade1(pred,predSeedling,bint,...
   yield,dDist,death,TIME1,TIME2,LEN,ISTART,SHAKE,INVADE)

%[nbarrec, rg, mechs, egrow,ripK,covP,covC,meanP,meanC] = JC_invade1(pred,predSeedling,bint,...
%   yield,dDist,death,TIME1,TIME2,LEN,ISTART,SHAKE,INVADE)
%
%This function runs the main simulation for the paper.  It can be run as an invasion analysis, simulate dynamics in a non-equilibrium setting.  The function takes the following input variables:
%
%pred- This is a SPPx1 vector.  Each value is the amount that an adult of each species would lower conspecific seedling survival if it was 0m away.  It is essentially potential NDD times .752.  It affects the log-likelihood chance of survival, and thus corresponds to p-tilde, rather than p.
%predSeedling- This is a SPPx1 vector.  Each value is the amount that one seed decreases the survival of all conspecific seeds at the same site.  The potential NDD is thus this number times yield (that corresponds to p-tilde, rather than p).
%bint- This is a SPPx1 vector.  It is the log-likelihood chance that a seed will survive if it is not exposed to any predators.  This corresponds to I-tilde, rather than I.
%yield- This is a SPPx1 vector.  It is the number of seeds each adult produces.
%dDist- This is a SPPx1 vector.  It is the parameter used for calculating seed dispersal.  It should be input as exp(alpha), where alpha is the alpha parameter in a 2Dt distribution.
%death- This is a SPPx1 vector, which gives the chance an adult tree will die each time step.
%TIME1- This is a variable.  It is the time step in which the invader is introduced.
%TIME2- This is a variable.  It is the number of time steps the simulation runs after the invader is introduced.  
%LEN- This is either a variable, or a 1xN vector.  If it is a variable, then it corresponds to the length of the environment (L in the main text).  In this case, the environment is set up as a LENxLEN grid, where each species other than the invader occurs at equal distribution, and is randomly distributed across the environment.  If it is a vector, then it is the initial distribution of the species.  In this case, LEN's length must be a square number.
%ISTART- This is the number of invaders that are introduced at time TIME1.
%SHAKE- Make this 1 to run the "shaken model" (i.e. where the community is randomly rearranged each time step), or 0 to allow spatial structure to form.
%INVADE- This is the species number of the invader.  If set to 0, then the simulation runs normally with no invader.
%
%The function outputs the following:
%nbarrec- This is a TIMExSPP matrix.  It tracks the frequency of each species.  To make a graph of the population dynamics, you can simply type plot(nbarrec).
%rg- This is a TIMExSPP matrix.  While the main simulation is running, I also run a DE version of the shaken approximation in the background (whether or not SHAKE=1; this does not affect the simulation).  rg is the frequency of each species in the DE model. 
%mechs- This is a SPPx(TIME-1)x6 matrix.  The value mechs(i,j,x) gives the value for calculating Delta mechanism x for species i at time j.  To calculate the actual Delta terms, take the value for the invader and subtract the average value of the resident.  mechs(i,j,1) gives ln(YI) (for Delta Y).  mech(i,j,2) gives N_j*mean(p) (for Delta P).  mechs(i,j,3) gives (1-n)*(p)*(sum(d*dP)) (for Delta kappa_P.  mechs(i,j,4) gives the covariance used to make Delta kappa'_P.  mechs(i,j,5) gives the approximation for E[cov(d,ln{C})], eqn A55 (Delta kappa_C).  mechs(i,j,6) gives the covariance used for Delta kappa'_C.
%egrow- This is a SPPx(TIME-1) matrix.  It gives the expected fitness of each species at each time step (i.e. lambda-tilde).  
%ripK- This is a 9xSPP matrix.  The value ripK(i,j) gives the pair-correlation value of species i at distance j.
%covP- This is a SPPx(TIME-1) matrix.  It gives the actual covariance between the number of seeds of sp. j arriving at a site and the predation risk at that site for each species at each time step.
%covC- This is a SPPx(TIME-1) matrix.  It gives the actual covariance between the number of seeds of sp. j arriving at a site and the total number of seedlings at that site for each species at each time step.
%meanP- This is a (TIME-1)xSPP matrix.  This is the average risk of predation across all sites for each species.
%meanC- This is a (TIME-1)x1 vector.  It is the average number of seedlings (post-predation) at each site.
%
%The function also contains code to generate a movie of a simulation.  To turn this on, you need to alter the code.  There is a place early on where it says MOVIE=0; if you set this to MOVIE=1, it will generate a movie.




MOVIE=0;


TIME=TIME1+TIME2;  %length of the simulation
NUMFORDISP=5;   %leave this as 5; if is how far we track dispersal before we switch to uniform dispersal

SPP=length(pred);   %the number of species in the simulation.
TREESIZE=10;     %this is l, the length of a site.


%DISPMAT gives the fraction of seeds that disperse 0 through NUMFORDISP
%sites away, plus the fraction that are long-distance dispersed.
DISP_MAT=zeros(SPP,NUMFORDISP+2);

for i=1:SPP
    disp(i,:) = dispersal_kern(dDist(i),TREESIZE,NUMFORDISP);
end


%Here we set up the initial distribution of species.  If the variable
%initial is just a number, then we generate a random environment.  If
%initial is a matrix, then we use it as the initial setup.
if(length(LEN)==1)
    
    AREA=LEN^2;
    
    if(INVADE==0)
        N=ceil([1:AREA]*SPP/AREA);%randi(SPP,1,AREA);
    else
        N=ceil([1:AREA]*(SPP-1)/AREA);
        N(N>=INVADE)=N(N>=INVADE)+1;
    end
    N=shuffle(N);
else
    N=LEN;
    AREA=length(N);
    LEN=AREA^.5;
end

%Here we set up a bunch of matrices for recording data.

record = zeros(TIME,AREA);   %record(i,j) gives the species number of the adult in site j at time i.
record(1,:)=N;		%This records the initial distribution of adults.
record_guess=zeros(TIME,SPP);  %record_guess(i,j) gives the frequency of species j at time i, running the shaken approximation.  
mechs=zeros(SPP,TIME-1,6);   %This records the Delta mechanisms.
meanP=zeros(TIME-1,SPP);    %This records the mean predation risk for each species.
meanC=zeros(TIME-1,1);    %This records the mean number of seedlings competing for each site.
covP=zeros(SPP,TIME-1);   %This is the covariance between predation risk and number of seeds at a sitefor each species.
covC=zeros(SPP,TIME-1);    %This is the covariance between seedling competition and the number of seeds of a species at a site for each species.
egrow=zeros(SPP,TIME-1);    %This is the expected growth rate of each species.


%predScale is the dispersal kernel for distance-responsive predators.  
predScale=exp(-.2*[.25, (1:NUMFORDISP)]*TREESIZE);

%This is a matrix of the intercept values for seedling survival for each species.
intMat=(repmat(bint',AREA,1));

%here I make matricies that list which which sites are a given distance
%away.  D0 is all of the sites a distance 0 away, so they are simply a list
%of all sites in order.  D1 will be a 4xAREA matrix, where D1(:,1) are all
%of the neighbors of site 1, D1(:,25) are the neighbors of site 25, etc.
%If WARP=0, then I assume that things fall off the edge; if WAPR=1, then I
%assume the forest is on a hypertorroid and things warp around the edge.
%In the non-warping case, neigbors which are off the edge are set to
%(AREA+1).  
D0=[1:AREA];

D1=make_neighbor(1,LEN,1);
D2=make_neighbor(2,LEN,1);
D3=make_neighbor(3,LEN,1);
D4=make_neighbor(4,LEN,1);
D5=make_neighbor(5,LEN,1);

N2=reshape(N,LEN,LEN);  %This is N (distribution of spp), but in a square.

%scale(z) is the number of sites a distance z-1 away.
scale=[1 4 8 12 16 20];

%pmod is the derivative of seedling survival at the intercept.  It is used to convert the value of pred and predSeedling input into the program into the potential NDD values used for calculating each mechanism. 
pmod=1./(exp(bint)+1);

%These are the potential NDD as defined in the main text.
distMax=pred.*pmod/.752;
densMax=predSeedling.*pmod.*yield;

%p combines the impact of both potential NDD's.
p=(distMax+densMax);
p(p==0)=10^(-8);

%dP is the weighted mean of the predator dispersal kernels (as described in the Appendix).
dP=(distMax*predScale+densMax.*disp(:,1:6))./p;

%ddP is the product of the seed and predator dispersal kernel, summed over all distances, for each species.
ddP=sum(disp(:,1:6).*dP.*scale,2);

%d2 is the seed dispersal kernel squared, summed over all distances, for each species.
d2=sum(disp(:,1:6).^2.*scale,2);

%d2dP is the product of the seed dispersal kernel squared and the predator dispersal kernel, summed over all distances, for each species.
d2dP=sum(disp(:,1:6).^2.*dP.*scale,2);

%ddBar is the product of each species' dispersal kernel and the across-species mean dispersal kernel, summed over all distances.
ddBar=sum(disp(:,1:6).*mean(disp(:,1:6)).*scale,2);

%logy is the log of the yield times the log of the density-independent survival. 
logy=log(yield.*exp(bint)./(1+exp(bint)));

%logYIstar is the mean of the product of yield and density-independent survival for all species.
logYIstar=log(mean(exp(logy)));

%meany is the mean number of seeds produced by each species.
meany=mean(exp(logy));

%esurvive is the density-independent survival of each species.
esurvive=exp(bint)./(1+exp(bint));

%dtemp is the dispersal kernel, ignoring the seeds that are dispersed uniformly around the environment.
dtemp=disp(:,1:6);



%%%%%%%%%%%%


%If you are making a movie, this sets it up.
if(MOVIE)
    
    if(exist('movname')==0)
        movname='moviename';
    end
    DOTSIZE=ceil(600/LEN);
    %movname='movieNeut'
    M=VideoWriter(movname);
    M.FrameRate = 7;
    M.Quality = 25;
    open(M);
    axis tight
    set(gca,'nextplot','replacechildren');
    framenum=1; %tracks which frame we are on.
    forX=repmat([1:LEN],1,LEN);
    forY=zeros(1,LEN)+1;
    for i=2:LEN
        forY=[forY, ones(1,LEN)*i];
    end
    colors=manyCol(SPP);
    
    if(INVADE>0)
        colors(1,:)=colors(INVADE,:);
        colors(INVADE,:)=0;
    end
end

%nbar is the frequency of each species.
for i=1:SPP
    nbar(i)=mean(N==i);
end
record_guess(1,:)=nbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tt=2:TIME
    
    %If you are running the shaken model, shuffle the population now.  If this is the first time step, the program reminds you that you are running a shaken model.  
    if(SHAKE)
        N=shuffle(N);
        if(tt==2)
             'shake'
        end
    end
    
    %If you are making a movie, this records the distribution of adults.  
    if(MOVIE)
        
        tt
        clf
        hold on
        for i=1:SPP
            plot(forX(N==i),forY(N==i),'.',...
                'Color',colors(i,:),'MarkerSize',DOTSIZE)
        end
        framenum=getframe;
    end
    
    
    
    
    %First I calculate the mean density of each species across the
    %landscape.
    for i=1:SPP
        nbar(i)=mean(N==i);
    end
    

    %%%%%%%
    %Here I calculate the seeds that go to each site
    
    %next, I calculate seed_fall, which is an AREAxSPP matrix of the mean
    %number of seeds of each species that are expected to fall in each
    %site.
    for i=1:SPP
            
        %Basically, sum(N(D1)==i)) gives a AREAx1 vector, where each
        %entry is the number of sp. i individuals that are one site
        %away from a given site (i.., if site 1 has 2 neighbors of sp.
        %1, then sum(N(D1)==1)) will be 2).  I multiply this by yield
        %and the appropriate disp value to calculate the expected
        %number of seeds that will disperse in from 1 site away.
        %
        %I also add in yield*nbar*LDD to represent the number of seeds that
        %are being long-distance dispersed from around the community.
        
        neighborI=[N(D0)==i;sum(N(D1)==i); sum((N(D2)==i));...
            sum((N(D3)==i)); sum((N(D4)==i)); sum((N(D5)==i))];
        
        %t=disp(i,1:(NUMFORDISP+1));
        
        seed_fall(:,i) = yield(i)*(disp(i,1:(NUMFORDISP+1))*neighborI+...
            disp(i,(NUMFORDISP+2))*nbar(i));
        %return
        
        pred_effect(:,i)=pred(i)*predScale*neighborI;
        pred_effectSeedling(:,i)=seed_fall(:,i)*predSeedling(i);
        
    end
    

    %seed_fall is then converted from an average number of falling seeds to
    %an actual number of falling seeds by chosing values from a Poisson
    %distribution.
    %seed_fall=poissrnd(seed_fall);
    
    %here I calculate the chance that a seed actually survives.
    logit_survive=pred_effect+pred_effectSeedling+intMat;
    survive=exp(logit_survive)./(1+exp(logit_survive));

    %seedling is the number of seedlings who survive.
    seedling=seed_fall.*survive;

    
    
    %C is a vector of number of seedlings at each site.
    C=seedling*ones(SPP,1);  %effect of compettiion
     
    
    %%%%%%%
    %Next I choose a seed at random.  If the adult dies, then that seed
    %becomes the new adult.  If the adult dies, then this is ignored.
    
    this_rand = rand(AREA,1).*C;
    
    here2 = SPP+1-(this_rand<seedling(:,1));
    
    for i=2:SPP
        seedling(:,i) = seedling(:,i) + seedling(:,i-1);
        here2 = here2-(this_rand<seedling(:,i));
    end
    
    %lives is 1 if an adult survives, and 0 if it dies.
    lives = (rand([AREA, 1])>death(N))' ;  
    
    N=here2'.*(1-lives) + N.*lives;  %whatever number is left, that is who is there.
    
    
    record(tt,:)=N;

    
    %%%%%%%%%%%%%%%%%%%
    %record mechanisms

    DeltaY=logy;
    mechs(:,tt-1,1)=DeltaY;

    DeltaP=(nbar'.*(-p));
    mechs(:,tt-1,2)=DeltaP;
    
    DeltaKp=((1-nbar').*(-p).*ddP);
    mechs(:,tt-1,3)=DeltaKp;

    predTemp=log(survive./esurvive');
    DeltaKpX=-(mean(seed_fall.*predTemp./mean(seed_fall))'-mean(predTemp)'...
        +DeltaKp);
    mechs(:,tt-1,4)=DeltaKpX;

    zeta=-(-p).*(dtemp.*predScale-2.*nbar'.*dtemp.*predScale+...
        nbar'.*dtemp+nbar'.*predScale)+dtemp.*(1+logy-log(mean(exp(logy))));
    forCov=zeta- sum(nbar'.*zeta);
    forCov2=sum(disp(:,1:6).*scale.*forCov,2);
    DeltaKc=forCov2-nbar*forCov2;
    mechs(:,tt-1,5)=DeltaKc;
    
    DeltaKcX=mean(log(C).*seed_fall./mean(seed_fall))'-mean(log(C))-DeltaKc;
    mechs(:,tt-1,6)=DeltaKcX;
    

    %%%%%%%%%%
    %Here I record how the population changes, using the DE of the shaken model.
    
    nb=record_guess(tt-1,:);

    DeltaY=logy;

    DeltaP=(nb'.*(-p));
    
    DeltaKp=((1-nb').*(-p).*ddP);

    zeta=-(-p).*(dtemp.*predScale-2.*nb'.*dtemp.*predScale+...
        nb'.*dtemp+nb'.*predScale)+dtemp.*(1+logy-log(mean(exp(logy))));
    forCov=zeta- sum(nb'.*zeta);
    forCov2=sum(disp(:,1:6).*scale.*forCov,2);
    DeltaKc=forCov2-nb*forCov2;
    
    tilde=DeltaY-DeltaP-DeltaKp-DeltaKc;

    record_guess(tt,:)=nb'.*(1+death.*(tilde-mean(tilde.*nb')*SPP));
    
    
    %%%%%%%%%%
    %Here I record important quantities
    
    seedling=seed_fall.*survive;
    C=seedling*ones(SPP,1);
    
    covC(:,tt-1)=mean(log(C).*seed_fall./mean(seed_fall))-mean(log(C));
    covP(:,tt-1)=mean(seed_fall.*predTemp./mean(seed_fall))-mean(predTemp);
    meanC(tt-1)=mean(log(C));
    egrow(:,tt-1)=1-death'+death'.*mean(seedling./C)./nbar;
    meanP(tt-1,:)=(1-mean(survive)./(exp(bint')./(1+exp(bint'))));
    
    
    
    %%%%%%%%%%%%%%%
    %here I add invaders
    
    
    if((INVADE>0)*(tt==TIME1))
        %I select ISTART sites at random, and place them there.
        for i=1:ISTART
            here=randi(AREA);
            while(N(here)==INVADE)
                here=randi(AREA);
            end
            N(here)=INVADE;
        end
        
        %I also update this for the DE community.
        record_guess(TIME1,:)=record_guess(TIME1,:)*(1-ISTART/AREA);
        record_guess(TIME1,INVADE)=ISTART/AREA;
        
        
        %If I'm recording a movie, I make 5 slides where all of the residents fade out a bit, and the invaders are displace in big dots.
        if(MOVIE)
            
            for i=1:5
                clf
                hold on
                for i=1:SPP
                    plot(forX(N==i),forY(N==i),'.',...
                        'Color',(1-.3*(1-colors(i,:))),'MarkerSize',DOTSIZE)
                end
                plot(forX(N==INVADE),forY(N==INVADE),'k.',...
                    'MarkerSize',2*DOTSIZE)
                
                framenum=getframe;
                %here I save the current graph to my movie file, M.
                writeVideo(M,framenum);
            end
            colors=(1-.3*(1-colors));
            colors(INVADE,:)=0;
        end
        record(TIME1,:)=N;
    end

end
    
    
%Here I record the mean frequency of each species over time.
nbarrec=mean(record==1,2);

for i=2:SPP
    nbarrec(:,i)=mean(record==i,2);
end

%Here I record the mean frequency according to the DE expectation.  I set a lower bound of 1 individual on the population.
rg=record_guess;
rg(rg<(1/AREA))=0;


TIMEMAT=([1:TIME]-1)*death(1);

if(MOVIE)
    
    colors=manyCol(SPP);
    
    if(INVADE>0)
        colors(1,:)=colors(INVADE,:);
        colors(INVADE,:)=0;
    end
    
    clf
    semilogy(TIMEMAT,nbarrec(:,1),'-','Color',colors(1,:))
    hold on
    semilogy(TIMEMAT,rg(:,1),'--','Color',colors(1,:))
    for (i=2:SPP)
        semilogy(TIMEMAT,nbarrec(:,i),'-','Color',colors(i,:))
        semilogy(TIMEMAT,rg(:,i),'--','Color',colors(i,:))
    end
    hold off
    
    framenum=getframe;
    %here I save the current graph to my movie file, M.
    writeVideo(M,framenum);
    
    close(M)
end


%Here I run ripleyK to calculate the pair correlation values in ripK 
ripleyK;

