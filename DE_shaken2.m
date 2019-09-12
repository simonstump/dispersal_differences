%This function runs the "shaken" DE approximation.  It is designed to be run after running JC_invade1, and uses its data.  

DRAW=1;



TIME=TIME1+TIME2;   %length of the simulation.
NUMFORDISP=5;   %leave this as 5; if is how far we track dispersal before we switch to uniform dispersal

SPP=length(pred);   %the number of species in the simulation.
TREESIZE=10;     %this is l, the length of a site.


%DISPMAT gives the fraction of seeds that disperse 0 through NUMFORDISP
%sites away, plus the fraction that are long-distance dispersed.
DISP_MAT=zeros(SPP,NUMFORDISP+2);

for i=1:SPP
    %disp(i,:) = dispersal_kern(dDist(1),TREESIZE,NUMFORDISP);
    disp(i,:) = dispersal_kern(dDist(i),TREESIZE,NUMFORDISP);
end

AREA=LEN^2;

%scale(z) is the number of sites a distance z-1 away.
scale=[1 4 8 12 16 20];

%predScale is the dispersal kernel for distance-responsive predators.  
predScale=exp(-.2*[.25, (1:NUMFORDISP)]*TREESIZE);
predScale=predScale/(sum(predScale.*scale));

%This is a matrix of the intercept values for seedling survival for each species.
intMat=(repmat(bint',AREA,1));

%pmod is the derivative of seedling survival at the intercept.  It is used to convert the value of pred and predSeedling input into the program into the potential NDD values used for calculating each mechanism. 
pmod=1./(exp(bint)+1);%

%These are the potential NDD as defined in the main text.
distMax=pred.*pmod/.752;
densMax=predSeedling.*pmod.*yield;

%p combines the impact of both potential NDD's.
p=(distMax+densMax);
p(p==0)=10^(-8);
dP=(distMax*predScale+densMax.*disp(:,1:6))./p;


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

%dtemp is the dispersal kernel, ignoring the seeds that are dispersed uniformly around the environment.
dtemp=disp(:,1:6);



%%%%%%%%%%%%
%Here are the matrixes that record information for the simulation.

%the frequency of each species at each time
record_guess(1,:)=record(1,:);

%My estimate for the covariance between the number of seeds at a site and the probability of predation there.
covP_guess=covP*0;

%My estimate for the covariance between the number of seeds at a site and the amount of seedling competition at that site.
covC_guess=covC*0;

%My estimate for the average amount of predation across sites.
meanP_guess=covP_guess;

%My estimate for the growth rate of each species.
egrow_guess=zeros(SPP,TIME1+TIME2-1);

%My estimate for the mean level of seedling competition across sites.
meanC_guess=zeros(TIME1+TIME2-1,1);

DelP=zeros(TIME1+TIME2-1,SPP);
DelP_guess=DelP;
DelKP=DelP;
DelKP_guess=DelP;
DelKC_guess=DelP;
DelKC=DelP;


pp=-p;

for tt=2:TIME
    %tt
    
        
    nbar=record_guess(tt-1,:);
    
    %Here I calculate the Delta-Tilde values (described in the Appendix).  
    DeltaY=logy;

    DeltaP=(nbar'.*(-p));
    
    DeltaKp=((1-nbar').*(-p).*ddP);

    zeta=-(-p).*(dtemp.*predScale-2.*nbar'.*dtemp.*predScale+...
        nbar'.*dtemp+nbar'.*predScale)+dtemp.*(1+logy-log(mean(exp(logy))));
    forCov=zeta- sum(nbar'.*zeta);
    forCov2=sum(disp(:,1:6).*scale.*forCov,2);
    DeltaKc=forCov2-nbar*forCov2;

    %tilde is the sum of the Delta terms.
    tilde=DeltaY-DeltaP-DeltaKp-DeltaKc;

    %We expect each species' per capita growth rate to be (lambda-tilde)*death+1.  The lambda-tilde terms are the amount that each species benefits or is harmed by each of the Delta terms.  As described in the appendix, we lambda-tilde will equal each species' tilde minus the mean tilde across species (weighted by abundance).
    record_guess(tt,:)=nbar'.*(1+death.*(tilde-mean(tilde.*nbar')*SPP));
    
    
    %%%%%%%%%%%%%%%%%%%%%
    %estimate terms

    %In addition to running the simulation, we estimate if our approximations for each covariance is correct.  For this, we use the actual abundances from the simulation, rather than our estimate based on the DE.  If we used the DE estimates, then errors could be due either to a bad approximation or because the frequencies were wrong.  By using the frequencies from the simulation, we are able to just look at the quality of the approximation.
    
    nbar=record(tt-1,:);
    zeta=-(-p).*(dtemp.*predScale-2.*nbar'.*dtemp.*predScale+...
        nbar'.*dtemp+nbar'.*predScale)+dtemp.*(1+logy-log(mean(exp(logy))));
    forCov=zeta- sum(nbar'.*zeta);
    forCov2=sum(disp(:,1:6).*scale.*forCov,2);
    DeltaKc=forCov2-nbar*forCov2;

    covC_guess(:,tt-1)=DeltaKc;

    covP_guess(:,tt-1)=-((1-nbar').*(-p).*ddP);


    %%%%%%%%%%%%%%
    %If you are simulating an invasion analysis, put the invader in now.
    
    if((INVADE>0)*(tt==TIME1))

        record_guess(TIME1,:)=record_guess(TIME1,:)*(1-ISTART/AREA);
        record_guess(TIME1,INVADE)=ISTART/AREA;
        
        record_guess2(TIME1,:)=record_guess2(TIME1,:)*(1-ISTART/AREA);
        record_guess2(TIME1,INVADE)=ISTART/AREA;
    end

end
    
    

rg=record_guess;
rg(rg<(1/AREA))=0;

if(DRAW)
    bob=manyCol(8);
    figure();
    plot(covC(1,:),covC_guess(1,:),'k.')
    
    hold on
    for i=2:8
        plot(covC(i,:),covC_guess(i,:),'.','Color',bob(i,:))
    end
    plot([-.1 .1],[-.1 .1],'k--')
    hold off
    legend([num2str([1:8]')])
    INVADE


    figure();
    plot(covP(1,:),covP_guess(1,:),'k.')
    
    hold on
    for i=2:8
        plot(covP(i,:),covP_guess(i,:),'.','Color',bob(i,:))
    end
    plot([-.1 .1],[-.1 .1],'k--')
    hold off
    legend([num2str([1:8]')])
    
end