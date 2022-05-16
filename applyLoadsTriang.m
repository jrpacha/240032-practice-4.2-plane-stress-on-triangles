function Q=applyLoadsTriang(nodes,elem,nodLoads,Q,forceLoad)
%------------------------------------------------------------------------
% (c) Numerical Factory 2019
%------------------------------------------------------------------------
% Apply a load BC on a boundary of a quadrilateral meshed domain. 
% The behabiour is very similar to the applyConvecQuad funtion used to
% apply convection BC.
%-------------------------------------------------------------------------
numElem=size(elem,1); 
[numNod,ndim]=size(nodes);
numCov=size(nodLoads,2);
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
    L=norm(nodes(n1,:)-nodes(n2,:));
    posit=[n1*ndim-1, n1*ndim];
    Q(posit)=Q(posit)+0.5*L*forceLoad; %same (x,y) force on both nodes
    posit=[n2*ndim-1, n2*ndim];   
    Q(posit)=Q(posit)+0.5*L*forceLoad; %same (x,y) force on both nodes
  end
end
end