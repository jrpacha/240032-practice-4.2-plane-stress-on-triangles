function Ke = planeElastTriangStiffMatrixJR(nodes,elem,e,C,h)
%
% For testing purposes only. Do not use this!
%
%plaNeElastTriangStiffMatrixJR
% INPUT
% nodes: matrix with the nodes' position (as usual)
%  elem: connectivity matrix (as usual)
%     e: element number
%     C: stress-strain matrix
%     h: thicknes (for strain problems h = 1)
%
% OUTPUT
%    Ke: element stifness matrix
%

nods = elem(e,:);
x = nodes(nods,1);
y = nodes(nods,2);

beta = [y(2)-y(3), y(3)-y(1), y(1)-y(2)];
gamma = [x(3)-x(2), x(1)-x(3), x(2)-x(1)];

Area = 0.5*det([ones(3,1),x,y]);

B = [
     beta(1),  0,        beta(2),  0,        beta(3),  0; ...
     0,        gamma(1), 0,        gamma(2), 0,        gamma(3); ...
     gamma(1), beta(1),  gamma(2), beta(2),  gamma(3), beta(3) ...
     ];
coef=0.25*h/Area; 
Ke = coef*B'*C*B;

end % end of function planeElastTriangStiffMatrixJR

