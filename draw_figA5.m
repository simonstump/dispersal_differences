%This code draws fig. A6.  To generate the data, run the code for_emergent_tradeoff.m.  

clc, clear

load('test_emergentTradeoff')


figA=figure();
plot(startCor,endCor,'k.','MarkerSize',15)
hold on
plot([-1 1],[-1 1],'k--')


xlabel('Initial correlation between $p_j^A$ and $\alpha_j$','interpreter','latex')
ylabel('Final correlation between $p_j^A$ and $\alpha_j$','interpreter','latex')
title('(a) Change in correlation between $p_j^A$ and $\alpha_j$','interpreter','latex')


set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA5a.pdf')

total=0;
for i=1:COMS
    
    XX=isAlive(i,:);
    
    temp=pred(permPred(i,:));
    
    vals=temp(XX);
    
    same=(vals==(sort(vals,'descend')));
    %same2=(val
    
    if(min(same)==1)
        i
        total=total+1;
    end
    
end
total

'mean start cor'
mean(startCor)

'std start cor'
std(startCor)

'mean end cor'
mean(endCor)

'corr is reduced'
mean(endCor<startCor)

'end corr is negative'
mean(endCor<0)

%%%%%%%%%%%%%%%%%%%

i=34;
%i=129

figA=figure(); 

XX=isAlive(i,:);

temp=-pred(permPred(i,:));

plot(temp(XX),alpha(XX),'ko',temp(~XX),alpha(~XX),'kx',...
    'MarkerSize',8,'LineWidth',3)

startCor(i)
endCor(i)

%a1=polyfit(temp,alpha,1);
%a2=polyfit(temp(XX),alpha(XX),1);
%
%x=[min(temp),max(temp)];
%hold on
%plot(x,x*a1(1)+a1(2),'k--',x,x*a2(1)+a2(2),'k-','LineWidth',2)

xlabel('Potential NDD, $p_j^A$','interpreter','latex')
ylabel('Dispersal parameter, $\alpha_j$','interpreter','latex')
title('(b) Change in $p_j^A$ and $\alpha_j$ via community assembly','interpreter','latex')

l=legend(['Persisting';'Excluded  ']);
set(l,'Interpreter','latex','Location','SouthWest');



set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA5b.pdf')