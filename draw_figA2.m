%This draws figure A2.  To generate the data for it, run:
%run_approx_check(2)
%run_approx_check(3)
%run_approx_check(6)

clc, clear


load('20190614_approx_check_justDisp')

DE_shaken2
close all

bob=manyCol(8);

figA=figure();
plot([1:TIME],record(:,1),'k.','LineWidth',3)
hold on
plot([1:TIME],record_guess(:,1),'k--','LineWidth',2)

hold on
for i=2:8
    plot([1:TIME],record(:,i),'.','Color',bob(i,:),'LineWidth',3)
    plot([1:TIME],record_guess(:,i),'--','LineWidth',2,'Color',bob(i,:))

end
hold off

axis([1,TIME,0,.25])

title('(a) Differences in dispersal','interpreter','latex')

xlabel('Time, $t$','interpreter','latex')
ylabel('Frequency, $\overline{N_j(t)}$','interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA2a.pdf')

figA=figure();
plot(covC(1,:),covC_guess(1,:),'k.')

hold on
for i=2:8
    plot(covC(i,:),covC_guess(i,:),'.','Color',bob(i,:))
end
plot([-.1 .1],[-.1 .1],'k--')
hold off


xlabel('Actual, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,C_x(t))$','interpreter','latex')
ylabel('Estimated, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,C_x(t))$','interpreter','latex')
title('(c) Differences in dispersal','interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA2c.pdf')



figA=figure();

plot(covP(1,:),covP_guess(1,:),'k.')

hold on
for i=2:8
    plot(covP(i,:),covP_guess(i,:),'.','Color',bob(i,:))
end
plot([-.1 .1],[-.1 .1],'k--')
hold off

title('(b) Differences in dispersal','interpreter','latex')
xlabel('Actual, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,\mbox{(predation)}_x)$','interpreter','latex')
ylabel('Estimated, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,\mbox{(predation)}_x)$','interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA2b.pdf')



%%%%%%%%%%%%%%%%%%%%%


clc, clear


load('20190614_approx_check_justPred')

DE_shaken2
close all

bob=manyCol(8);

figA=figure();
plot([1:TIME],record(:,1),'k.','LineWidth',3)
hold on
plot([1:TIME],record_guess(:,1),'k--','LineWidth',2)

hold on
for i=2:8
    plot([1:TIME],record(:,i),'.','Color',bob(i,:),'LineWidth',3)
    plot([1:TIME],record_guess(:,i),'--','LineWidth',2,'Color',bob(i,:))

end
hold off

axis([1,TIME,0,.3])

title('(d) Differences in predation','interpreter','latex')
xlabel('Time, $t$','interpreter','latex')
ylabel('Frequency, $\overline{N_j(t)}$','interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA2d.pdf')

figA=figure();
plot(covC(1,:),covC_guess(1,:),'k.')

hold on
for i=2:8
    plot(covC(i,:),covC_guess(i,:),'.','Color',bob(i,:))
end
plot([-.1 .1],[-.1 .1],'k--')
hold off

title('(f) Differences in predation','interpreter','latex')
xlabel('Actual, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,C_x(t))$','interpreter','latex')
ylabel('Estimated, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,C_x(t))$','interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA2f.pdf')



figA=figure();

plot(covP(1,:),covP_guess(1,:),'k.')

hold on
for i=2:8
    plot(covP(i,:),covP_guess(i,:),'.','Color',bob(i,:))
end
plot([-.1 .1],[-.1 .1],'k--')
hold off

title('(e) Differences in predation','interpreter','latex')
xlabel('Actual, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,\mbox{(predation)}_x)$','interpreter','latex')
ylabel('Estimated, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,\mbox{(predation)}_x)$','interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA2e.pdf')



%%%%%%%%%%%%%%%%%%%


clc, clear


load('20190614_approx_check_rand')

DE_shaken2
close all

bob=manyCol(8);

figA=figure();
plot([1:TIME],record(:,1),'k.','LineWidth',3)
hold on
plot([1:TIME],record_guess(:,1),'k--','LineWidth',2)

hold on
for i=2:8
    plot([1:TIME],record(:,i),'.','Color',bob(i,:),'LineWidth',3)
    plot([1:TIME],record_guess(:,i),'--','LineWidth',2,'Color',bob(i,:))

end
hold off

axis([1,TIME,0,.4])
title('(g) Different dispersal \& susceptibility 1','interpreter','latex')
xlabel('Time, $t$','interpreter','latex')
ylabel('Frequency, $\overline{N_j(t)}$','interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA2g.pdf')

figA=figure();
plot(covC(1,:),covC_guess(1,:),'k.')

hold on
for i=2:8
    plot(covC(i,:),covC_guess(i,:),'.','Color',bob(i,:))
end
plot([-.1 .1],[-.1 .1],'k--')
hold off

title('(i) Different dispersal \& susceptibility 1','interpreter','latex')
xlabel('Actual, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,C_x(t))$','interpreter','latex')
ylabel('Estimated, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,C_x(t))$','interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA2i.pdf')



figA=figure();

plot(covP(1,:),covP_guess(1,:),'k.')

hold on
for i=2:8
    plot(covP(i,:),covP_guess(i,:),'.','Color',bob(i,:))
end
plot([-.1 .1],[-.1 .1],'k--')
hold off

title('(h) Different dispersal \& susceptibility 1','interpreter','latex')
xlabel('Actual, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,\mbox{(predation)}_x)$','interpreter','latex')
ylabel('Estimated, $\mbox{\textbf{cov}}(\mbox{(\# seeds)}_x,\mbox{(predation)}_x)$','interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'figA2h.pdf')

