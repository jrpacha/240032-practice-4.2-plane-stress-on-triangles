clearvars
close all

nodes=[0,0;
       120, 0;
       120, 160;
       0, 160];
elem=[1 2 3;
      3 4 1];

[numNod,ndim]=size(nodes);
numElem=size(elem,1);

%plot Elements
numbering=1; %=0;
plotElements(nodes, elem, numbering);

%Real constants
E=3.e7;    %N/mm^2
nu=0.25;
th=3.6e-2; %Thickness (in mm^2)
c11=E/(1-nu^2);
c22=c11;
c12=nu*c11;
c33=0.5*E/(1+nu);

C=[c11, c12, 0; c12, c22, 0; 0, 0, c33];
K = zeros(ndim*numNod);
Q = zeros(ndim*numNod,1);
F = zeros(ndim*numNod,1);

for e = 1:numElem
    v1=nodes(elem(e,1),:);
    v2=nodes(elem(e,2),:);
    v3=nodes(elem(e,3),:);
    beta=[v2(2)-v3(2),v3(2)-v1(2),v1(2)-v2(2)];
    gamma=-[v2(1)-v3(1),v3(1)-v1(1),v1(1)-v2(1)];
    Area=0.5*det([v1 1; v2 1; v3 1]);
    B=[beta(1), 0, beta(2), 0, beta(3), 0; 
       0, gamma(1), 0, gamma(2), 0, gamma(3);
       gamma(1), beta(1), gamma(2), beta(2), gamma(3), beta(3)]/(2*Area);
    fprintf('\nMatrix B of element %d:\n',e)
    fprintf('%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n',B')
    Ke = Area*th*B'*C*B;
    fprintf('\nMatrix Ke of element %d:\n',e)
    fprintf('%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e\n',Ke')
    row=[2*elem(e,1)-1; 2*elem(e,1); ...
         2*elem(e,2)-1; 2*elem(e,2); ...
         2*elem(e,3)-1; 2*elem(e,3)];
    col=row;
    K(row,col)=K(row,col)+Ke;
end

%Apply BC
%Constant traction on the right (Natural BC)
t0=1.0e3; %in N/mm, so we shall not multiply by th. in Qe
nodLoads=[2,3]; 
L23=norm(nodes(3,:)-nodes(2,:));
Qe=0.5*t0*L23*[1;0;1;0];%<-- We do not multiply by thickness 
row=[2*nodLoads(1)-1;2*nodLoads(1);2*nodLoads(2)-1;2*nodLoads(2)];
Q(row)=Q(row)+Qe;
%Fix displacements (essential BC)
u=zeros(2*numNod,1);
LeftNod=[1,4];
fixedNod=[2*LeftNod-1,2*LeftNod];
u(fixedNod)=0.0;

%Reduced System
freeNod=setdiff(1:2*numNod,fixedNod);
%Fm=F(freeNod)-K(freeNod,fixedNod)*u(fixedNod); %Note: here u(fixedNod)=0,
                                                %so this is not necessary
                                                %in this case
%Fm=Fm+Q(freeNod);
Qm=Q(freeNod);
Km=K(freeNod,freeNod);

%Solve the reduced system
um=Km\Qm;
u(freeNod)=um;

%Post Process
%Compute the Stress
stress=zeros(numElem,3);
VonMisses=zeros(2,1);
for e=1:numElem
    v1=nodes(elem(e,1),:);
    v2=nodes(elem(e,2),:);
    v3=nodes(elem(e,3),:);
    beta=[v2(2)-v3(2),v3(2)-v1(2),v1(2)-v2(2)];
    gamma=-[v2(1)-v3(1),v3(1)-v1(1),v1(1)-v2(1)];
    Area=0.5*det([v1 1; v2 1; v3 1]);
    B=[beta(1), 0, beta(2), 0, beta(3), 0; 
       0 gamma(1), 0 gamma(2), 0 gamma(3);
       gamma(1), beta(1), gamma(2), beta(2), gamma(3), beta(3)]/(2*Area);   
    row=[2*elem(e,1)-1; 2*elem(e,1); ...
         2*elem(e,2)-1; 2*elem(e,2); ...
         2*elem(e,3)-1; 2*elem(e,3)];
    sigma=C*B*u(row);
    stress(e,:)=sigma';
    VonMisses(e,1)=sqrt(sigma(1)^2+sigma(2)^2-sigma(1)*sigma(2)+3*sigma(3)^2);
end

%Output
%Displacements:
fprintf('\n%35s\n\n','Displacements')
fprintf('%7s%8s%12s%12s%12s\n',...
    'Num.Nod.','X','Y','U','V')
fprintf('%4d%16.4e%12.4e%12.4e%12.4e\n', ...
    [(1:numNod)',nodes(:,1),nodes(:,2),u(1:2:end),u(2:2:end)]')
%Stress:
fprintf('\n%31s\n\n','Stress')
fprintf('%7s%10s%12s%12s%11s\n','Elem.','SXX','SYY','SXY','VM')
fprintf('%4d%16.4e%12.4e%12.4e%12.4e\n',[(1:numElem)',stress,VonMisses]')

%Graphical output
title='Desp.X';
colorScale='jet';
valueToShow=u(1:2:end);
plotContourSolution(nodes,elem,valueToShow,title,colorScale);
%Strain
esc=100;
plotPlaneNodElemDespl(nodes, elem, u, esc)
%Von Misses Stress
%plotStressVM(nodes,elem, VonMisses)