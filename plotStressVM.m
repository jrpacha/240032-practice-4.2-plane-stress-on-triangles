function plotStressVM(nodes,elem, VonMisses)
% Plot VonMisses Stress  for each Triangle
numElem=size(elem,1);
figure()
xmax=max(nodes(:,1));
xmin=min(nodes(:,1));
ymax=max(nodes(:,2));
ymin=min(nodes(:,2));
shift=2*max(0.02*max(abs([xmax, xmin,ymin,ymax])),0.01);

for e=1:numElem
    v1=nodes(elem(e,1),:);
    v2=nodes(elem(e,2),:);
    v3=nodes(elem(e,3),:);
    vertices=[v1;v2;v3];
    X=vertices(:,1);
    Y=vertices(:,2);
    val=VonMisses(e);
    values=[val;val;val]; %equal for the three vetices beacuse is constant.
    plot(X,Y,'k')
    fill(X,Y,values);
    hold on;
end
   axis([xmin-2*shift, xmax+2*shift, ymin-2*shift, ymax+2*shift]);
   axis equal;
   title('Von Misses Stress') ;         
colorbar;
