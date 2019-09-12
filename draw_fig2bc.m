%This draws Fig. 2b and 2c.  To do it, first run the code run_for_distribution.m.  

clc, clear, close all

load('20190614holdDist')


this=reshape(x1(:,1),LEN,LEN);
INVADE
nval(1)

[a,b]=size(this);

theX=repmat([1:a],b,1);
theY=repmat([1:b]',1,a);

theX2=reshape(theX,a*b,1);
theY2=reshape(theY,a*b,1);

%[x2 y2]

num1=(this==1);
figA=figure();
plot(theX2(num1),theY2(num1),'k.','MarkerSize',12)
pbaspect([1 1 1])
axis off

N=reshape(num1,LEN^2,1);
ripleyK
ripK(3,1)

title(['High dispersal, $\overline{N_j(t)}=$',num2str(nval(1))],'interpreter','latex')

set(gca,'fontsize', 12);



set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'fig2b1.pdf')


'-------------------'

this=reshape(x2(:,1),LEN,LEN);
INVADE
nval(1)

[a,b]=size(this);

theX=repmat([1:a],b,1);
theY=repmat([1:b]',1,a);

theX2=reshape(theX,a*b,1);
theY2=reshape(theY,a*b,1);

%[x2 y2]

num1=(this==1);
figA=figure();
plot(theX2(num1),theY2(num1),'k.','MarkerSize',12)
pbaspect([1 1 1])
axis off

N=reshape(num1,LEN^2,1);
ripleyK
ripK(3,1)

title(['High dispersal, $\overline{N_j(t)}=$',num2str(nval(2))],'interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'fig2b2.pdf')


'-------------------'
this=reshape(x3(:,1),LEN,LEN);
INVADE
nval(1)

[a,b]=size(this);

theX=repmat([1:a],b,1);
theY=repmat([1:b]',1,a);

theX2=reshape(theX,a*b,1);
theY2=reshape(theY,a*b,1);

%[x2 y2]

num1=(this==INVADE2);
figA=figure();
plot(theX2(num1),theY2(num1),'k.','MarkerSize',12)
pbaspect([1 1 1])
axis off

N=reshape(num1,LEN^2,1);
ripleyK
ripK(3,1)

title(['Low dispersal, $\overline{N_j(t)}=$',num2str(nval(1))],'interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'fig2c1.pdf')


'-------------------'
this=reshape(x4(:,1),LEN,LEN);
INVADE
nval(1)

[a,b]=size(this);

theX=repmat([1:a],b,1);
theY=repmat([1:b]',1,a);

theX2=reshape(theX,a*b,1);
theY2=reshape(theY,a*b,1);

%[x2 y2]

num1=(this==INVADE2);
figA=figure();
plot(theX2(num1),theY2(num1),'k.','MarkerSize',12)
pbaspect([1 1 1])
axis off

N=reshape(num1,LEN^2,1);
ripleyK
ripK(3,1)

title(['Low dispersal, $\overline{N_j(t)}=$',num2str(nval(2))],'interpreter','latex')

set(gca,'fontsize', 12);

set(figA,'Units','Inches');
pos = get(figA,'Position');
set(figA,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
saveas(figA,'fig2c2.pdf')


'-------------------'
z=3

[K1(z,INVADE1) K2(z,INVADE1) K3(z,INVADE2) K4(z,INVADE2)]