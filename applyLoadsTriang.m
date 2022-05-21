function Q=applyLoadsTriang(nodes,elem,nodLoads,Q, forceLoad)
%------------------------------------------------------------------------
% (c) Numerical Factory 2018
%
% Apply a convection BC on a boundary of a triangulated domain. 
%
% Input:
% nodLoads -> List of nodes on the load boundary. They can be
%          non-connected but it needs: 
%          i) Each connected component contains at least two nodes 
%          ii) No corners are accepted (two edges of the same triangle). 
%          If an unique triangle is a corner of the load bondary, we have two call twice 
%          this function (one time with each edge) and sume contributions
%
%          *           *               
%          | \         |               
%          |  \    =   |   +       
%          *---*       *     *---*  
%
% Q    -> Global Secondary variables Vector
% nodes-> all nodes (matrix of size nunNodx2)
% elem -> all elements (matrix of size numElemx3)
%
%------------------------------------------------------------------------
% Output:
% Q -> Modified Q with loads
%
%------------------------------------------------------------------------
numElem=size(elem,1); 
[numNod,ndim]=size(nodes);
numCov=length(nodLoads);
if numCov==1 
	error('applyLoadTriang: Not unic node allow'); 
end
for k=1:numElem
  aux=[0,0,0]; %initial values (no node is found)
  for inod=1:3 %loop for the three nodes
      r=find(nodLoads==elem(k,inod)); %find if one node of element k is in the convection node list 
      if(~isempty(r)) 
          aux(inod)=1; %put 1 on the local position found
      end
  end
    %      
    % Now the aux vector is one of the following possibilities
    % aux=[0,0,0]  -> this element does not contain any convection node
    % aux=[1,0,0], [0,1,0], [0,0,1] -> contains only one convection node 
    % aux=[1,1,0], [1,0,1], [0,1,1] -> contains two convection nodes
    % To classify which is the edge in local numbering, we consider aux as
    % a binary number expressed as power of 2: aux(1)*2^0+aux(2)*2^1+aux(3)*2^2
    % only 3 (first edge), 5 (third edge) or 6 (second edge)
    % are significant.
  number = aux(1)+2*aux(2)+4*aux(3); 
  switch (number) %identify the appropriate edge
      case 3
          ij=[1,2];
      case 5
          ij=[3,1];
      case 6
          ij=[2,3];
      case 7
          error('applyLoadTriang: Corners not allowed !!!!\n');           
      otherwise, ij=[0,0];
  end
  if ( ij(1) > 0) %it's an existing edge 
    n1=elem(k,ij(1)); 
    n2=elem(k,ij(2));
    h=norm(nodes(n1,:)-nodes(n2,:));
    posit=[n1*ndim-1, n1*ndim];
    Q(posit)=Q(posit)+0.5*h*forceLoad; %same (x,y) force on both nodes
    posit=[n2*ndim-1, n2*ndim];   
    Q(posit)=Q(posit)+0.5*h*forceLoad; %same (x,y) force on both nodes
  end
end