clearvars
close all

clearvars
close all

ndiv = 15;
t = 20;
scale=t/50;
signe=1;

delta = -0.4;

triangle = [0, 0;  ...
            10, 0; ...
            5, 10; ...
            0, 0];

edgeLength = norm(triangle(1,:)-triangle(2,:));
ndivEdgeFixed = 68;
h = 1/ndivEdgeFixed;

plot(triangle(:,1),triangle(:,2),'-','LineWidth',2,'Color','black')
hold on 

axis equal
axis off

for i=0:h:1
    p = (1-i)*triangle(1,:)+i*triangle(2,:);
    line = [p; p+delta];
    plot(line(:,1),line(:,2),'-k')
end

plotEdgeConstantBC(triangle(2,:),triangle(3,:),t,ndiv,scale,signe,...
    'Color','black','LineWidth',1.5)
text(0.5,0.4,'$1$','FontSize',14,'Interpreter','latex','Color','red')
text(9.2,0.4,'$2$','FontSize',14,'Interpreter','latex','Color','red')
text(4.85,9.0,'$3$','FontSize',14,'Interpreter','latex','Color','red')
text(9.0,6.0,'$\tau_{0}$','FontSize',18,'Interpreter','latex',...
    'Color','black')
hold off
saveas(gcf,'figure1.png')