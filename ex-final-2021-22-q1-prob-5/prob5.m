% 240032 ExFinal Q1 2022-23
% Problema 5

clearvars
close all

E = 9.0e8;               % N/m^2
nu = 0.52;               % Poisson Ratio
forceLoad = [0;-1075.0]; % N/m
elemStressVM = 421;

%Define the plane elasticity problem: 
%modelProblem=1; %plane stress
modelProblem=2; %plane strain


%eval('AirFoilmesh01');
%save 'AirFoilmesh01.mat' nodes elem -mat
load AirFoilmesh01.mat;

numNodes = size(nodes,1);
numElem = size(elem,1);
ndim = size(nodes,2);

%figure()
numbering=0;
plotElements(nodes, elem, numbering);
hold on
%Boundary nodes
indNodBd = boundaryNodes(nodes, elem);
indNodLeft = find(nodes(:,1) < -1.99);
indNodRight = find(nodes(:,1) > 11.99);
indNodBottom = find(nodes(:,2) < -1.99);
indNodTop = find(nodes(:,2) > 2.99);

%indNodRBT = unique([indNodTop',indNodBottom',indNodRight']);
indNodExternalBd = unique([indNodTop;indNodBottom;indNodRight;indNodLeft]);
indNodInternalBd = setdiff(indNodBd,indNodExternalBd);

plot(nodes(indNodTop,1),nodes(indNodTop,2),'o','MarkerFaceColor',...
    'red','MarkerSize',10)
plot(nodes(indNodBottom,1),nodes(indNodBottom,2),'o','MarkerFaceColor',...
    'green','MarkerSize',10)
plot(nodes(indNodInternalBd,1),nodes(indNodInternalBd,2),...
    'o','MarkerFaceColor','blue','MarkerSize',10)
hold off

switch modelProblem
    case 1
        c11=E/(1-nu^2);
        c22=c11;
        c12=nu*c11;
        c21=c12;
        c33=E/(2*(1+nu));
        problemType = 'Plane stress problem';
    case 2
        th=1.0;
        c11=E*(1-nu)/((1+nu)*(1-2*nu));
        c22=c11;
        c12=c11*nu/(1-nu);
        c21=c12;
        c33=E/(2*(1+nu));
        problemType = 'Plane strain problem';
    otherwise
        error('modelProblem should be 1 (stress) or 2 (strain)');
end
C=[c11, c12, 0; c21, c22, 0; 0, 0, c33];

clc
fprintf('\tPROBLEM 5 (%s)\n',problemType)
%PART (A)
numNodInternalBd = length(indNodInternalBd);
avXintBd = sum(nodes(indNodInternalBd,1))/numNodInternalBd;
fprintf(['Part (a)\n',...
         '*** Number of nodes in the airfoil shaped boundary of\n',...
         '*** the domain: %d\n', ...
         '*** Hint. The mean of the x-component of the nodes in\n',...
         '*** this boundary is: %.6e\n\n'],numNodInternalBd,avXintBd)



K=zeros(ndim*numNodes);  
F=zeros(ndim*numNodes,1);
Q=zeros(ndim*numNodes,1);

for e=1:numElem
    Ke=planeElastTriangStiffMatrix(nodes,elem,e,C,th);
    %
    % Assemble the stiffness matrices
    %
    row=[ndim*elem(e,1)-1; ndim*elem(e,1); ...
         ndim*elem(e,2)-1; ndim*elem(e,2); ...
         ndim*elem(e,3)-1; ndim*elem(e,3)];
    col=row;
    K(row,col)=K(row,col)+Ke;
end
%Boundary conditions
%Natural BC: constant traction on the right boundary
nodLoads=indNodBottom'; %nodes where the load is applied (transposed)
Q=applyLoadsTriang(nodes,elem,nodLoads,Q,forceLoad);
%Essential BC: left boundary fixed
u=zeros(ndim*numNodes,1);
fixedNod=[ndim*indNodInternalBd'-1,ndim*indNodInternalBd',ndim*indNodTop'];
u(fixedNod)=0.0;

%Plot BC
u1=[min(nodes(indNodBottom,1));nodes(indNodBottom(1),2)];
u2=[max(nodes(indNodBottom,1));nodes(indNodBottom(1),2)];
ndiv = 50;
scale = 0.35;
plotEdgeConstantBC(u1,u2,-forceLoad(2),ndiv,scale)

%Reduced System
freeNod=setdiff(1:ndim*numNodes,fixedNod);
%Fm=F(freeNod)-K(freeNod,fixedNod)*u(fixedNod); %Note: here u(fixedNod)=0,
                                                %so this is not necessary
                                                %in this case
%Fm=Fm+Q(freeNod);
Qm=Q(freeNod);
Km=K(freeNod,freeNod);

%Solve the reduced system
um=Km\Qm;
u(freeNod)=um;

%PART (B)
fprintf(['\nPart (b)\n',...
         '*** The maximum of the absolute value of the vertical\n',...
         '*** displacement of the nodes is max |UY| = %.4e\n',...
         '*** Hint: The maximum of the absolute value of the\n',...
         'horizontal displacement of the nodes is: max |UX| = %.4e\n\n'],...
         max(abs(u(2:2:end))),max(abs(u(1:2:end)))) 

%PART (C)
displ = zeros(numElem,1);
for e = 1:numElem
    nods = elem(e,:);
    UX = u(ndim*nods-1);
    UY = u(ndim*nods);
    displ(e) = norm([UX(1);UY(1)]) + ...
        norm([UX(2);UY(2)]) + ...
        norm([UX(3);UY(3)]);
end
[maxDispl,elemMaxDispl]=max(displ);
fprintf(['Part (c)\n',...
         '*** Elem of max. displacement: %d\n\n'],elemMaxDispl)

%PART (D)
%Post Process
%Compute the Stress & Strain
stress=zeros(numElem,3);
strain=zeros(numElem,3);
vonMisses=zeros(numElem,1);
for e=1:numElem
    v1=nodes(elem(e,1),:);
    v2=nodes(elem(e,2),:);
    v3=nodes(elem(e,3),:);
    beta=[v2(2)-v3(2),v3(2)-v1(2),v1(2)-v2(2)];
    gamma=-[v2(1)-v3(1),v3(1)-v1(1),v1(1)-v2(1)];
    Area=0.5*det([v1 1; v2 1; v3 1]);
    B=[beta(1), 0, beta(2), 0, beta(3), 0; 
       0, gamma(1), 0 gamma(2), 0 gamma(3);
       gamma(1), beta(1), gamma(2), beta(2), gamma(3), beta(3)]/(2*Area);   
    row=[2*elem(e,1)-1; 2*elem(e,1); ...
         2*elem(e,2)-1; 2*elem(e,2); ...
         2*elem(e,3)-1; 2*elem(e,3)];
    ue=u(row,:);
    strain(e,:)=B*ue;   %This is constant for each element
    stress(e,:)=C*B*ue; %This as constant for each element
    sxx=stress(e,1);
    syy=stress(e,2);
    sxy=stress(e,3);
    vonMisses(e)=sqrt(sxx^2+syy^2-sxx*syy+3*sxy^2);
end
fprintf(['Part (d)\n',...
        '*** The Von Misses stress of the element %d is\n',...
        '*** stressVM = %.4e\n\n'],elemStressVM,vonMisses(elemStressVM))
%Plot the domain and its deformation
esc=2.5e4;
plotPlaneNodElemDespl(nodes, elem, u, esc)
%Plot VM stress
plotStressVM(nodes,elem,vonMisses)