%This draws figures A4 and A5.  To generate the data for it, run the following code:
% run_holdAbund(1)
% run_holdAbund(2)
% run_holdAbund(3)

close all
clc, clear

load('20190612_multiHold_tradeoff')

z=3;

temp=mean(res1(z,:,:,1:rep),4)-1;
temp2=reshape(temp,length(nval),SPP);


X=(log(nval));

pred=nval.*[1:SPP]'*0;
comp=nval.*[1:SPP]'*0;
L=length(nval);
for i=1:SPP
    x=mean(res3(4,:,i,i,1:rep),5);
    pred(i,:)=reshape(x,1,L);
    x=mean(res3(6,:,i,i,1:rep),5);
    comp(i,:)=reshape(x,1,L);    
end

figA=figure();
semilogx(temp2,comp','k-','LineWidth',3)

title('(a) Dispersal-susceptibility trade-off','interpreter','latex')
xlabel('Clustering, ln$\{g(3)-1\}$','interpreter','latex')
ylabel('$\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,C_x(t))$','interpreter','latex')

axis([.1 10 -.01 .06])
set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA4a.pdf')

figA=figure();
semilogx(temp2,pred','k-','LineWidth',3)

title('(a) Dispersal-susceptibility trade-off','interpreter','latex')
xlabel('Clustering, ln$\{g(3)-1\}$','interpreter','latex')
ylabel('$\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,\mbox{(predation)}_x)$','interpreter','latex')

axis([.1 10 -.01 .06])
set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA3a.pdf')



%%%%%%%%%%%%%%%%%


load('20190612_multiHold_justDisp')

z=3;

temp=mean(res1(z,:,:,1:rep),4)-1;
temp2=reshape(temp,length(nval),SPP);


X=(log(nval));

pred=nval.*[1:SPP]'*0;
comp=nval.*[1:SPP]'*0;
L=length(nval);
for i=1:SPP
    x=mean(res3(4,:,i,i,1:rep),5);
    pred(i,:)=reshape(x,1,L);
    x=mean(res3(6,:,i,i,1:rep),5);
    comp(i,:)=reshape(x,1,L);    
end

figA=figure();
semilogx(temp2,comp','k-','LineWidth',3)

title('(b) Differences in dispersal','interpreter','latex')
xlabel('Clustering, ln$\{g(3)-1\}$','interpreter','latex')
ylabel('$\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,C_x(t))$','interpreter','latex')

axis([.1 10 -.01 .06])
set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA4b.pdf')

figA=figure();
semilogx(temp2,pred','k-','LineWidth',3)

title('(b) Differences in dispersal','interpreter','latex')
xlabel('Clustering, ln$\{g(3)-1\}$','interpreter','latex')
ylabel('$\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,\mbox{(predation)}_x)$','interpreter','latex')

axis([.1 10 -.01 .06])
set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA3b.pdf')



%%%%%%%%%%%%%%%%%


load('20190612_multiHold_justPred')

z=3;

temp=mean(res1(z,:,:,1:rep),4)-1;
temp2=reshape(temp,length(nval),SPP);


X=(log(nval));

pred=nval.*[1:SPP]'*0;
comp=nval.*[1:SPP]'*0;
L=length(nval);
for i=1:SPP
    x=mean(res3(4,:,i,i,1:rep),5);
    pred(i,:)=reshape(x,1,L);
    x=mean(res3(6,:,i,i,1:rep),5);
    comp(i,:)=reshape(x,1,L);    
end

figA=figure();
semilogx(temp2,comp','k-','LineWidth',3)

title('(c) Differences in predation','interpreter','latex')
xlabel('Clustering, ln$\{g(3)-1\}$','interpreter','latex')
ylabel('$\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,C_x(t))$','interpreter','latex')

axis([.1 10 -.01 .06])
set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA4c.pdf')

figA=figure();
semilogx(temp2,pred','k-','LineWidth',3)

title('(c) Differences in predation','interpreter','latex')
xlabel('Clustering, ln$\{g(3)-1\}$','interpreter','latex')
ylabel('$\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,\mbox{(predation)}_x)$','interpreter','latex')

axis([.1 10 -.01 .06])
set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA3c.pdf')



%%%%%%%%%%%%%%%%%


