function Td = connectDOFs(n_el,n_nod,n_i,Tnod)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_el     Total number of elements
%                  n_nod    Number of nodes per element
%                  n_i      Number of DOFs per node
%   - Tn    Nodal connectivities table [n_el x n_nod]
%            Tn(e,a) - Nodal number associated to node a of element e
%--------------------------------------------------------------------------
% It must provide as output:
%   - Td    DOFs connectivities table [n_el x n_el_dof]
%            Td(e,i) - DOF i associated to element e
%--------------------------------------------------------------------------
% Hint: Use the relation between the DOFs numbering and nodal numbering.

Td=zeros(n_el,n_nod*n_i);
for i=1:n_el
    Td(i,1)=3*Tnod(i,1)-2;
    Td(i,2)=3*Tnod(i,1)-1;
    Td(i,3)=3*Tnod(i,1);
    Td(i,4)=3*Tnod(i,2)-2;
    Td(i,5)=3*Tnod(i,2)-1;
    Td(i,6)=3*Tnod(i,2);
end

end