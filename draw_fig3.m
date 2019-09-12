%This draws fig. A3.  To generate the data for it, run the following code:
%run_invasion(1,0)
%run_invasion(2,0)
%run_invasion(3,0)
%run_invasion(4,0)
%run_invasion(6,0)
%run_invasion(7,0)

close all
clc, clear
load('20190603_invade_tradeoff')

figA=figure();
plot((mean(res2)*death(1))+1,log(mean(res5))*death(1)+1,'ko','MarkerSize',8,'LineWidth',3)
hold on
plot([-.02 2],[-.02 2],'k--','LineWidth',2)


title('(d) Dispersal-susceptibility trade-off','interpreter','latex')
ylabel('Growth rate with spatial structure, $\tilde{\lambda}_i$','interpreter','latex')
xlabel('Growth rate without spatial structure, $\tilde{\lambda}_i$','interpreter','latex')

set(gca,'fontsize', 14);

axis([.995, 1.025, .995, 1.025])

'for the trade-off'
[mean(res2); yield'; pred';dDist']

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'fig3d.pdf')

zz=mean(res3(:,:,4))+ mean(res3(:,:,6));

discrepancy(1)=mean(log(mean(res5))-mean(res2));

isItDeltas(1)=mean(log(mean(res5))-mean(res2)+zz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('20190603_invade_justPred')

figA=figure();
plot((mean(res2)*death(1))+1,log(mean(res5))*death(1)+1,'ko','MarkerSize',8,'LineWidth',3)
hold on
plot([-.02 2],[-.02 2],'k--','LineWidth',2)


title('(c) Differences in susceptibility','interpreter','latex')
ylabel('Growth rate with spatial structure, $\tilde{\lambda}_i$','interpreter','latex')
xlabel('Growth rate without spatial structure, $\tilde{\lambda}_i$','interpreter','latex')

set(gca,'fontsize', 14);

axis([.995, 1.025, .995, 1.025])

'for the predators only'
[mean(res2); yield'; pred';dDist']


set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'fig3c.pdf')

zz=mean(res3(:,:,4))+ mean(res3(:,:,6));

discrepancy(2)=mean(log(mean(res5))-mean(res2));

isItDeltas(2)=mean(log(mean(res5))-mean(res2)+zz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('20190603_invade_justDisp')

figA=figure();
plot((mean(res2)*death(1))+1,log(mean(res5))*death(1)+1,'ko','MarkerSize',8,'LineWidth',3)
hold on
plot([-.02 2],[-.02 2],'k--','LineWidth',2)


title('(b) Differences in dispersal','interpreter','latex')
ylabel('Growth rate with spatial structure, $\tilde{\lambda}_i$','interpreter','latex')
xlabel('Growth rate without spatial structure, $\tilde{\lambda}_i$','interpreter','latex')

set(gca,'fontsize', 14);

axis([.995, 1.025, .995, 1.025])

'for the dispersal differences'
[mean(res2); yield'; pred';dDist']


set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'fig3b.pdf')

zz=mean(res3(:,:,4))+ mean(res3(:,:,6));

discrepancy(3)=mean(log(mean(res5))-mean(res2));

isItDeltas(3)=mean(log(mean(res5))-mean(res2)+zz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('20190603_invade_justYield')

figA=figure();
plot((mean(res2)*death(1))+1,log(mean(res5))*death(1)+1,'ko','MarkerSize',8,'LineWidth',3)
hold on
plot([-.02 2],[-.02 2],'k--','LineWidth',2)


title('(a) Baseline with yield differences','interpreter','latex')
ylabel('Growth rate with spatial structure, $\tilde{\lambda}_i$','interpreter','latex')
xlabel('Growth rate without spatial structure, $\tilde{\lambda}_i$','interpreter','latex')

set(gca,'fontsize', 14);

axis([.995, 1.025, .995, 1.025])

'for yield differences'
[mean(res2); yield'; pred';dDist']


set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'fig3a.pdf')

zz=mean(res3(:,:,4))+ mean(res3(:,:,6));

discrepancy(4)=mean(log(mean(res5))-mean(res2));

isItDeltas(4)=mean(log(mean(res5))-mean(res2)+zz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('20190603_invade_rand')

figA=figure();
plot((mean(res2)*death(1))+1,log(mean(res5))*death(1)+1,'ko','MarkerSize',8,'LineWidth',3)
hold on
plot([-.02 2],[-.02 2],'k--','LineWidth',2)


title('(e) Different dispersal \& susceptibility 1','interpreter','latex')
ylabel('Actual invader growth rate, $\tilde{\lambda}_i$','interpreter','latex')
xlabel('Estimated invader growth rate, $\tilde{\lambda}_i$','interpreter','latex')

set(gca,'fontsize', 14);

axis([.995, 1.025, .995, 1.025])

'for the random 1'
[mean(res2); yield'; pred';dDist']


set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'fig3e.pdf')

zz=mean(res3(:,:,4))+ mean(res3(:,:,6));

discrepancy(5)=mean(log(mean(res5))-mean(res2));

isItDeltas(5)=mean(log(mean(res5))-mean(res2)+zz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('20190603_invade_rand2')

figA=figure();
plot((mean(res2)*death(1))+1,log(mean(res5))*death(1)+1,'ko','MarkerSize',8,'LineWidth',3)
hold on
plot([-.02 2],[-.02 2],'k--','LineWidth',2)


title('(f) Different dispersal \& susceptibility 2','interpreter','latex')
ylabel('Actual invader growth rate, $\tilde{\lambda}_i$','interpreter','latex')
xlabel('Estimated invader growth rate, $\tilde{\lambda}_i$','interpreter','latex')

set(gca,'fontsize', 14);

axis([.995, 1.025, .995, 1.025])

'for random 2'
[mean(res2); yield'; pred';dDist']


set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'fig3f.pdf')

zz=mean(res3(:,:,4))+ mean(res3(:,:,6));

discrepancy(6)=mean(log(mean(res5))-mean(res2));

isItDeltas(6)=mean(log(mean(res5))-mean(res2)+zz);



discrepancy
mean(discrepancy)

isItDeltas
mean(isItDeltas)