function Ke=planeElastTriangStiffMatrix(nodes,elem,e,C,th)
%function PLANEELASTTRIANGSTIFFMATRIX 
%Computes the stiff matrix of a triangrular element in the isotropic plane 
%elasticity problem 
% INPUT
% nodes: Matrix with the vertexs coordinates of the element, as usual.
%  elem: Connectivity matrix of the mesh.
%     e: Number of the current element, according to the connectivity.
%     C: Elasticity matrix (plane isotropic case).
%    th: Thicknes of the element. Recall that, for the model problem 2
%        (plane strain), the thickness is not significative and we take
%        th=1
% OUTPUT
%    Ke: Current elemt's stiff matrix.

v1=nodes(elem(e,1),:);
v2=nodes(elem(e,2),:);
v3=nodes(elem(e,3),:);
beta=[v2(2)-v3(2),v3(2)-v1(2),v1(2)-v2(2)];
gamma=-[v2(1)-v3(1),v3(1)-v1(1),v1(1)-v2(1)];
Area=0.5*det([v1 1; v2 1; v3 1]);
B=[beta(1), 0, beta(2), 0, beta(3), 0; 
   0, gamma(1), 0, gamma(2), 0, gamma(3);
   gamma(1), beta(1), gamma(2), beta(2), gamma(3), beta(3)]/(2*Area);    
Ke = Area*th*B'*C*B;
end %end of function PLANEELASTTRIANGSTIFFMATRIX

