function Fext = computeF(n_i,n_dof,Tmat,mat,x,Tnod,n_el,g,M)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  n_i         Number of DOFs per node
%                  n_dof       Total number of DOFs
%   - Fdata  External nodal forces [Nforces x 3]
%            Fdata(k,1) - Node at which the force is applied
%            Fdata(k,2) - DOF (direction) at which the force acts
%            Fdata(k,3) - Force magnitude in the corresponding DOF
%--------------------------------------------------------------------------
% It must provide as output:
%   - Fext  Global force vector [n_dof x 1]
%            Fext(I) - Total external force acting on DOF I
%--------------------------------------------------------------------------
% Hint: Use the relation between the DOFs numbering and nodal numbering to
% determine at which DOF in the global system each force is applied.

Fext=zeros(n_dof,1); 
%{
for i=1:size(Fdata,1)
    if Fdata(i,2) == 3 
        Fext(Fdata(i,1)*3,1)=Fdata(i,3);
    end
    if Fdata(i,2) == 2
       Fext(Fdata(i,1)*3-1,1)=Fdata(i,3);
    end
    if Fdata(i,2) == 1
       Fext(Fdata(i,1)*3-2,1)=Fdata(i,3);
    end
end
%}

EM=0; % masas de las barras (absoluto vertical descendiente)
Emom7=0; % momentos de las barras respecto nodo 7.
xnod7=x(7,1);
znod7=x(7,3);

for e=1:n_el % element weight
    x1e=x(Tnod(e,1),1);
    y1e=x(Tnod(e,1),2);
    z1e=x(Tnod(e,1),3);
    x2e=x(Tnod(e,2),1);
    y2e=x(Tnod(e,2),2);
    z2e=x(Tnod(e,2),3);
    le=sqrt((x2e-x1e)^2+(y2e-y1e)^2+(z2e-z1e)^2);
    m=le*mat(Tmat(e),2)*mat(Tmat(e),3); % m=le*E*A
    
    EM=EM+m;
   
    Emom7=Emom7+(x2e-xnod7)*(m*g/2)+(x1e-xnod7)*(m*g/2); %(y+)
    
    Fext(3*Tnod(e,1),1)=Fext(3*Tnod(e,1),1)-m*g/2; 
    Fext(3*Tnod(e,2),1)=Fext(3*Tnod(e,2),1)-m*g/2;
end
    
    
syms L D T; %(L+ en z, T+ en x, D+ en x)
eq1= L-M*g-EM*g==0; % 1. sumatorio fuerzas verticales
eq2= T+D==0;        % 2. sumatorio fuerzas horizontales

EmL7=0; %sumatorio momento distribución L
EmD7=0; %sumatorio momento distribución D
for n=3:7 %upper surface nodes
    EmL7=EmL7+(x(n,1)-xnod7)*(-L/5);  %(y+)
    EmD7=EmD7+(x(n,3)-znod7)*(D/5);   %(y+)
end

EmM7=0; %sumatorio momento distribución M
EmT7=0; %sumatorio momento distribución T
for n=1:2 %lower bar nodes
    EmM7=EmM7+(x(n,1)-xnod7)*(M*g/2); %(y+)
    EmT7=EmT7+(x(n,3)-znod7)*(T/2);  %(y+)
end

eq3=Emom7+EmL7+EmD7+EmM7+EmT7==0; % 3. sumatorio momentos (y+) respecto nodo 7.

S=vpasolve(eq1,eq2,eq3,L,D,T);

% complete Fext:
for n=3:7 %upper surface nodes
    Fext(3*n,1)=Fext(3*n,1)+vpa(S.L)/5;
    Fext(3*n-2,1)=Fext(3*n-2,1)+vpa(S.D)/5;
end

for n=1:2 %lower bar nodes
    Fext(3*n,1)=Fext(3*n,1)-M*g/2;
    Fext(3*n-2,1)=Fext(3*n-2,1)+vpa(S.T)/2;
end



end