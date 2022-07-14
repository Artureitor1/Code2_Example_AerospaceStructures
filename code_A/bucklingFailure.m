function [FB] = bucklingFailure(mat,Tmat,x,Tnod,n_el, sig, D1, d1, D2)
    %  Tn(e,a) = global nodal number associated to node a of element e
    %  x(a,j) = coordinate of node a in the dimension j
    %  Total number of elements
    %  Tmat(e) = Row in mat corresponding to the material associated to element e 
    
    % Material data
    %  mat(m,1) = Young modulus of material m
    %  mat(m,2) = Section area of material m
        
   
    FB = zeros(n_el,1); %Bars that have failed
    for e = 1:n_el

        x1e=x(Tnod(e,1),1);
        y1e=x(Tnod(e,1),2);
        z1e=x(Tnod(e,1),3);
        x2e=x(Tnod(e,2),1);
        y2e=x(Tnod(e,2),2);
        z2e=x(Tnod(e,2),3);
        le=sqrt((x2e-x1e)^2+(y2e-y1e)^2+(z2e-z1e)^2);

        E = mat(Tmat(e),1);
        A = mat(Tmat(e),2);
        I = mat(Tmat(e),4);
        
        sigPan = ((pi^2)*E*I)/((le^2)*A); % Formula of buckling failure tension
        
        if sig(e)<0
            if abs(sig(e)) > sigPan
                FB(e) = 1; %If the bars fails because buckling effect, 1 will appear
            end
        end
    end
end