clearvars
close all

%Material properties
th=0.5;                 %thickness (in mm)
%forceLoad=[1.0e3;0.0]; %[Fx; Fy] traction on the right boundary (in N/mm)
forceLoad = 1900.00;    % N/mm
E=1.0e6;                %Young's modulus (in N/mm^2)
nu=0.25;                %Poisson's ratio (adimensional)
nodHint = 200;

% Part (b)
%Define the plane elasticity problem: 
%modelProblem=1; %plane stress
%modelProblem=2; %plane strain (2)

modelProblem=1; %Plane stress
          
%eval('CircleHolemesh01');  %load the mesh 
%load CircleHolemesh01.mat  %We read the nodes and the connectivity
                            %matrix from a .mat
eval('meshTriangWHole');

[numNod,ndim]=size(nodes);
numElem=size(elem,1);

numbering=0; %=1 shows nodes and element numbering
plotElementsOld(nodes,elem,numbering);

%Find Boundary points
indNodCirc = find(sqrt(nodes(:,1).^2 + (nodes(:,2)+0.4).^2) < 0.201);
indRightEdge = find(nodes(:,2)+1.732*nodes(:,1) - 0.732 > -0.01);
indNodBottom = find(nodes(:,2) < -0.99);

hold on
plot(nodes(indNodCirc,1),nodes(indNodCirc,2),'o','color','black',...
    'markerFaceColor','green','markerSize',7)
plot(nodes(indRightEdge,1),nodes(indRightEdge,2),'o','color','black',...
    'markerFaceColor','blue','markerSize',7)
plot(nodes(indNodBottom,1),nodes(indNodBottom,2),'o','color','black',...
    'markerFaceColor','yellow','markerSize',7)
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
    row=[ndim*elem(e,1)-1; ndim*elem(e,1); ...
         ndim*elem(e,2)-1; ndim*elem(e,2); ...
         ndim*elem(e,3)-1; ndim*elem(e,3)];
    col=row;
    K(row,col)=K(row,col)+Ke;
end
%Boundary conditions
%

%Natural BC: constant traction on the right boundary
xi = [0; 0.732; 0] - [1; -1; 0];
eta = cross(xi, [0; 0; 1]);
normalVector = [eta(1); eta(2)];
normalVector = normalVector/norm(normalVector);
tractionForce = forceLoad * normalVector;
ndiv = 20;
scale = 0.8;
signe = 1.0;
plotEdgeConstantBC([1,-1],[0,0.732],forceLoad,ndiv,scale,signe)

%nodLoads=indNodBottom';
%Q=applyLoadsTriang(nodes,elem,nodLoads,Q,[0;-forceLoad]);

nodLoads=indRightEdge'; %nodes where the load is applied (transposed)
Q=applyLoadsTriang(nodes,elem,nodLoads,Q,tractionForce);

%Essential BC: left boundary fixed
u=zeros(ndim*numNod,1);

fixedNod = indNodCirc';
fixedNod = [ndim*fixedNod-1, ndim*fixedNod];
freeNod = setdiff(1:ndim*numNod,fixedNod);
u(fixedNod)=0.0;

%
%Reduced System
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


%fprintf('(b) Maximum of the horizontal displacements of the bottom boundary nodes: %.4e\n',...
%    max(displacements(indNodBottom,1)))
%
% Error: 
%
% The text should say "the vertical displacements", NOT "the horizontal
% displacementet"
%
fprintf('Problem 5(b) Maximum of the vertical displacements of the bottom boundary nodes: %.4e\n',...
     max(displacements(indNodBottom,2)))
fprintf('         Hint. The vertical displacement of node %d is: %.4e\n',...
     nodHint, displacements(nodHint, 2))

%Graphical output
%figure()
esc=2.0;
plotPlaneNodElemDespl(nodes, elem, u, esc);
% valueToShow=u(1:2:end);
% title='Desp.X';
% colorScale='jet';
% plotContourSolution(nodes,elem,valueToShow,title,colorScale);
% plotStressVM(nodes,elem,vonMisses);