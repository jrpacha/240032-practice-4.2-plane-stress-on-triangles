function plotEdgeConstantBC(v1,v2,t,ndiv,scale,signe,varargin)
%plot traction vector (origin at the edge)
long=norm(v2-v1);
edgeVector=(v2-v1)/long;
normal=[edgeVector(2),-edgeVector(1)];
f=signe*t*normal;
%ndiv=10;
h=1/ndiv;
j=0;
p=[];
for i=0:h:1
    j=j+1;
    p(j,:)=(1-i)*v1+i*v2;
end
X=p(:,1);
Y=p(:,2);
U=f(1)*ones(ndiv+1,1);
V=f(2)*ones(ndiv+1,1);
hold on;
%scale=0.2*long/t;
quiver(X,Y,U,V,scale,varargin{:});
axis equal;
hold off;
%end plot
