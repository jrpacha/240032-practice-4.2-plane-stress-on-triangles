function plotPlaneNodElemDespl(nodes, elem, u, esc)
[num_nod, dim]=size(nodes);
num_elem=size(elem,1);
nodu=zeros(num_nod,dim);

nodu(:,1)=nodes(:,1)+esc*u(1:2:num_nod*dim,1);
nodu(:,2)=nodes(:,2)+esc*u(2:2:num_nod*dim,1);
ymax=max(nodu(:,2));
ymin=min(nodu(:,2));
xmax=max(nodu(:,1));
xmin=min(nodu(:,1));
figure()
plotElementsOld(nodes, elem, 0);
hold on
plotElementsOld(nodu, elem, 0);
