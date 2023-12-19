function Ke=planeElastTriangStiffMatrix(nodes,elem,e,C,th)
v1=nodes(elem(e,1),:);
v2=nodes(elem(e,2),:);
v3=nodes(elem(e,3),:);
vertices=[v1;v2;v3];
beta=[v2(2)-v3(2),v3(2)-v1(2),v1(2)-v2(2)];
gamma=-[v2(1)-v3(1),v3(1)-v1(1),v1(1)-v2(1)];
%alpha=[v2(1)*v3(2)-v2(2)*v3(1),v3(1)*v1(2)-v3(2)*v1(1),v1(1)*v2(2)-v1(2)*v2(1)];
Area=1/2*det([ones(3,1),vertices]);
B=[beta(1),  0, beta(2), 0,     beta(3), 0;
   0 ,  gamma(1),0 ,   gamma(2),    0,  gamma(3); 
   gamma(1), beta(1),gamma(2), beta(2),gamma(3), beta(3)]/(2*Area);
Ke=(B'*C*B)*th*Area;
