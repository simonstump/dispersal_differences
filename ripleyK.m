%This calculates the pair correlation function for all species, for distances 1 through 9.    It puts it into a 9xSPP (#of species) matrix called ripK.

if(~exist('D1'))
D1=make_neighbor(1,LEN,1);
D2=make_neighbor(2,LEN,1);
D3=make_neighbor(3,LEN,1);
D4=make_neighbor(4,LEN,1);
D5=make_neighbor(5,LEN,1);
end
if(~exist('D9'))
D6=make_neighbor(6,LEN,1);
D7=make_neighbor(7,LEN,1);
D8=make_neighbor(8,LEN,1);
D9=make_neighbor(9,LEN,1);
end


for i=1:SPP
    %neigh is the number of neighbors a distance x from each adult of that species.
    neigh=[sum(N(D1)==i)/4;
        sum(N(D2)==i)/8;...
        sum(N(D3)==i)/12;...
        sum(N(D4)==i)/16;...
        sum(N(D5)==i)/20;...
        sum(N(D6)==i)/24;...
        sum(N(D7)==i)/28;...
        sum(N(D8)==i)/32;...
        sum(N(D9)==i)/36];
    
    %neigh2 is the expectation under uniformity.
    neigh2=neigh(:,(N==i));
    
    ripK(:,i)=mean(neigh2,2)/mean(N==i);
    
end

