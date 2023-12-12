clearvars
close all

%Material properties
th=0.5;                %thickness (in mm)
%forceLoad=[1.0e3;0.0];%[Fx; Fy] traction on the right boundary (in N/mm)
forceLoad = 1900.00;    % N/mm
E=1.0e+6;              %Young's modulus (in N/mm^2)
nu=0.25;               %Poisson's ratio (adimensional)
nodHint=200;


% Part (a)

%Define the plane elasticity problem: 
%modelProblem=1; %plane stress
%modelProblem=2; %plane strain 

modelProblem=1; %Plane stress problem
          
%eval('CircleHolemesh01');  %load the mesh 
%load CircleHolemesh01.mat  %We read the nodes and the connectivity
                            %matrix from a .mat
eval('meshTriangWHole');

[numNod,ndim]=size(nodes);
numElem=size(elem,1);

numbering=0; %=1 shows nodes and element numbering
plotElementsOld(nodes,elem,numbering);

%Find Boundary points
indNodBottom=find(nodes(:,2)<-0.99);
%indRight=find(nodes(:,1)>0.99);

hold on
plot(nodes(indNodBottom,1),nodes(indNodBottom,2),'o','color','black',...
    'markerFaceColor','green','markerSize',7)
%plot(nodes(indRight,1),nodes(indRight,2),'o','color','black',...
%    'markerFaceColor','blue','markerSize',7)
hold off

switch modelProblem
    case 1
        c11=E/(1-nu^2);
        c22=c11;
        c12=nu*c11;
        c21=c12;
        c33=E/(2*(1+nu));
        fprintf('Plane stress problem\n')
    case 2
        th=1.0;
        c11=E*(1-nu)/((1+nu)*(1-2*nu));
        c22=c11;
        c12=c11*nu/(1-nu);
        c21=c12;
        c33=E/(2*(1+nu));
        fprintf('Plane strain problem\n')
    otherwise
        error('modelProblem should be 1 (stress) or 2 (strain)');
end
C=[c11, c12, 0; c21, c22, 0; 0, 0, c33];

%Computation of the stiffness matrix
K=zeros(ndim*numNod);
F=zeros(ndim*numNod,1);
Q=zeros(ndim*numNod,1);
for e=1:numElem
    Ke=planeElastTriangStiffMatrix(nodes,elem,e,C,th);
    %
    % Assemble the stiffness matrices
    %
    row=[2*elem(e,1)-1; 2*elem(e,1); ...
         2*elem(e,2)-1; 2*elem(e,2); ...
         2*elem(e,3)-1; 2*elem(e,3)];
    col=row;
    K(row,col)=K(row,col)+Ke;
end
%Boundary conditions
%Natural BC: constant traction on the right boundary
nodLoads=indNodBottom'; %nodes where the load is applied (transposed)
Q=applyLoadsTriang(nodes,elem,nodLoads,Q,[0;-forceLoad]);
%Essential BC: left boundary fixed
u=zeros(ndim*numNod,1);

fixedNod = [1;3];
fixedNod=[ndim*fixedNod-1;ndim*fixedNod];
u(fixedNod)=0.0;

%Reduced System
freeNod=setdiff(1:ndim*numNod,fixedNod);
%Fm=F(freeNod)-K(freeNod,fixedNod)*u(fixedNod); %Note: here u(fixedNod)=0,
                                                %so this is not necessary
                                                %in this case
%Fm=Fm+Q(freeNod);
Qm=Q(freeNod);
Km=K(freeNod,freeNod);

%Solve the reduced system
um=Km\Qm;
u(freeNod)=um;

displacements = [u(1:2:end),u(2:2:end)];

fprintf('Problem 5(a) Maximum of the horizontal displacements of the bottom boundary nodes: %.4e\n',...
    max(displacements(indNodBottom,1)))
fprintf('             Hint. The horizontal displacement of node %d is: %.4e\n',...
    nodHint, displacements(nodHint, 1))

% %Graphical output
esc=2.0;
plotPlaneNodElemDespl(nodes, elem, u, esc);
% valueToShow=u(1:2:end);
% title='Desp.X';
% colorScale='jet';
% plotContourSolution(nodes,elem,valueToShow,title,colorScale);
% plotStressVM(nodes,elem,vonMisses);