clearvars
close all

%Material properties
E=3.0e+7; %Young Modulus (N/cm^2)
nu=0.3; %Poisson's ratio (adimensional)
th=1.0;  % thickness (cm)
t=20; %traction force N/cm^2 (in the normal direction of the 2nd edge).
nodelProblem=1;

%Geometry
v1=[0,0];
v2=[10,0];
v3=[5,10];
nodes=[v1;v2;v3];
elem=[1,2,3];
numNod=size(nodes,1);
ndim=size(nodes,2);
numElem=size(elem,1);
plotElements(nodes,elem,1);

modelProblem=1;
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

F=zeros(ndim*numNod,1);
%Q=zeros(ndim*numNod,1); %here pre-allocation of Q is not necessary
%Compute the Stiffness matrix
K=planeElastTriangStiffMatrix(nodes,elem,1,C,th);

%Boundary conditions
edgeVector=(v3-v2)/norm(v3-v2);
normal=[-edgeVector(2),edgeVector(1)];
f=t*normal; %constant traction applied on normal vector direction
%
%plot traction vector (origin at the edge)
%
ndiv=10; %number of arrows in the plot
vini=v2;
vfin=v3;
scale=t/100; %scale arrows according to t value
signe=1; %force orientation
plotEdgeConstantBC(vini,vfin,t,ndiv,scale,signe) %
%
%Natural B.C.: constant traction on the second edge of the triangle
%
Q=-0.5*th*norm(v3-v2)*[0;0;f(1);f(2);f(1);f(2)];
%
%Essential B.C.
fixedNod=[];
%nod 1;
nNod=1; 
fixedNod=[fixedNod,ndim*nNod-1]; %u1_x=0
fixedNod=[fixedNod,ndim*nNod];   %u1_y=0
%nod 2;
nNod=2; 
fixedNod=[fixedNod,ndim*nNod-1]; %u2_x=0
fixedNod=[fixedNod,ndim*nNod];   %u2_y=0
%
u=zeros(ndim*numNod,1);
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

%Post Process
%Compute the Stress & Strain
beta=[v2(2)-v3(2),v3(2)-v1(2),v1(2)-v2(2)];
gamma=-[v2(1)-v3(1),v3(1)-v1(1),v1(1)-v2(1)];
Area=0.5*det([v1 1; v2 1; v3 1]);
B=[beta(1), 0, beta(2), 0, beta(3), 0; 
   0 gamma(1), 0 gamma(2), 0 gamma(3);
   gamma(1), beta(1), gamma(2), beta(2), gamma(3), beta(3)]/(2*Area); 
strain=B*u;
stress=C*strain;
sxx=stress(1);
syy=stress(2);
sxy=stress(3);
vonMisses=sqrt(sxx^2+syy^2-sxx*syy+3*sxy^2);

%Output (fancy output: don't try this at the exams!)
%Displacements only the last 10 nodes:
displacements=[u(1:2:end),u(2:2:end)];
tableDispl=[(1:numNod)',nodes(:,1),nodes(:,2),displacements];
fprintf('\n%35s\n\n','Displacements')
fprintf('%7s%8s%12s%12s%11s\n',...
    'Num.Nod.','X','Y','U','V')
fprintf('%4d%16.4e%12.4e%12.4e%12.4e\n',tableDispl');    
%Stress:
tableStress=[(1:numElem)',stress',vonMisses]; 
fprintf('\n%31s\n\n','Stress')
fprintf('%5s%12s%12s%12s%15s\n','Elem.','SXX','SYY','SXY','vonMisses')
fprintf('%4d%16.4e%12.4e%12.4e%12.4e\n',tableStress')

%Graphical output: displacements
esc=5.0e4;
plotPlaneNodElemDespl(nodes, elem, u, esc);