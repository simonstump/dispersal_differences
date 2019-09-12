function out = dispersal_kern(dDist,SIZE,NUMFORDISP)

%function out = dispersal_kern(dDist,SIZE,NUMFORDISP)
%
%This function generates a 1x(NUMFORDISP+2) vector.  Numbers 1 through NUMFORDISP+1 are the probability that a seed will disperse to a given site a distance 0 through NUMFORDISP away.  The last entry is the probability that it will disperse to any site that is farther than NUMFORDISP sites away.
%
%This function assumes a 2Dt distribution.  dDist is exp(alpha), where alpha is the a dispersal parameter.  SIZE is the length of a site, and is 10 for most simulations.


% %%%%%%%for testing
% out(1)=.2;
% out(NUMFORDISP+2)=1-out(1);
% return


dx=.01;
x=[dx:dx:SIZE*(NUMFORDISP+1)];

%2D-t kernel, multiplied by the area of a ring
pdf=1./(pi.*dDist.*(1+x.^2./dDist).^2).*2.*pi.*x;

out(1)=sum(pdf.*(x<=.5*SIZE))*dx;

%%%%%%%for testing
%out(1)=out(1)*2;
%out(NUMFORDISP+2)=1-out(1);
%return

for(i=1:NUMFORDISP)
    out(i+1)=sum(pdf.*(x>(SIZE*(.5+i-1))).*(x<=(SIZE*(.5+i))))*dx;
end
%out

%long distance dispersal
out(length(out)+1)=1-sum(out);

%before this, out is the number of seeds per site X away.  Here we change
%it to the number of seeds dispersing to each site X away.s
numSpot=4*[1:NUMFORDISP];
out=out./[1, numSpot, 1];