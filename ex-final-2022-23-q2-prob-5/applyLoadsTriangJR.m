function Q=applyLoadsTriangJR(nodes,elem,nodLoads,Q,forceLoad)
%
% Just for testing purposes. Do not use this!
%
%applyLoadsTriangJR
%INPUT
%    nodes: matrix with nodes' position
%     elem: connectivity
%        h: thickness (h=1) for strain problems
% nodLoads: row vector with the indices of the loaded nodes 
%        Q: Loads
%forceLoad: column vector with the loads (Fx,Fy)
%
%OUTPUT
%         Q: vector Q at the output is the incoming Q plus with tractions 
%            (or compressions) added.
%
% Remark: as in applyConvecTriang.m, no corners are allowed

nodLoads=nodLoads(:)';

if size(nodLoads,2) < 2
    error('applyLoadTriangJR: a unique node not allowed')
end

ndim=size(nodes,2);

for e=1:size(elem,1)
    nodElem=elem(e,:);
    r=intersect(nodElem,nodLoads);
    if ~isempty(r)
        nnodes=size(r,2);
        if nnodes > 1
            if nnodes < 3 
                n1=r(1,1);
                n2=r(1,2);
                L=norm(nodes(n1,:)-nodes(n2,:));
                rows=[ndim*n1-1;ndim*n1];
                Q(rows) = Q(rows) + 0.5*L*forceLoad;
                rows=[ndim*n2-1;ndim*n2];
                Q(rows) = Q(rows) + 0.5*L*forceLoad;
            else
                error('applyLoadTriangJR: corners not allowed!!!');
            end
        end                
    end
end
end %end of function applyLoadsTriangJR

