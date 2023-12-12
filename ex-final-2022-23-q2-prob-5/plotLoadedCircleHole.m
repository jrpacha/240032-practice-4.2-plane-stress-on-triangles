clearvars
close all
forceLoad = [1000.0;0];
vertexs=[-1,-1;1,-1;1,1;-1,1;-1,-1];
plot(vertexs(:,1),vertexs(:,2),'k-','lineWidth',2)
axis([-1,2,-1.1,1.1])
set(gca,'xtick',[-1,-0.5,0,0.5,1]);
set(gca,'ytick',[-1,-0.5,0,0.5,1]);
set(gca,'box','off');
hold on
theta=linspace(0,2*pi,361);
radius=0.5;
xx=radius*cos(theta);
yy=radius*sin(theta);
plot(xx,yy,'k-','lineWidth',2)
ndiv=20;
scale=1;
signe=1;
plotEdgeConstantBC(vertexs(2,:),vertexs(3,:),forceLoad,ndiv,scale,signe, ...
    'Color','black','LineWidth',1.5)
text(1.51,0,'$\tau = 10^{3}\,\mathrm{N/mm}$','interpreter','LaTeX','fontSize',14)
hold off
saveas(gcf,'loadedCircle.png');
