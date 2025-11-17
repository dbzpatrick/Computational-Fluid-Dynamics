function [A] = findAmatrix(Nx, Ny, dx2, gamma2)

A = zeros((Nx-1)*(Ny-1), (Nx-1)*(Ny-1));

% Interior points
for j = 2:Ny-2
    for i = 2:Nx-2
        po = i+(j-1)*(Nx-1);
        A(po, po) = -(2 * gamma2 + 2)/dx2; %i, j
        A(po, po+1) = 1/dx2; %i+1, j
        A(po, po-1) = 1/dx2; %i-1, j
        A(po, po+(Nx-1)) = gamma2/dx2; %i, j+1
        A(po, po-(Nx-1)) = gamma2/dx2; %i, j-1
    end 
end 

% South
j = 1;
for i = 2:Nx-2
    po = i+(j-1)*(Nx-1);
    A(po, po) = -(2 * gamma2 + 2)/dx2; %i, j
    A(po, po+1) = 1/dx2; %i+1, j
    A(po, po-1) = 1/dx2; %i-1, j
    A(po, po+(Nx-1)) = gamma2/dx2; %i, j+1
end 

% North
j = Ny-1;
for i = 2:Nx-2
    po = i+(j-1)*(Nx-1);
    A(po, po) = -(2 * gamma2 + 2)/dx2; %i, j
    A(po, po+1) = 1/dx2; %i+1, j
    A(po, po-1) = 1/dx2; %i-1, j
    A(po, po-(Nx-1)) = gamma2/dx2; %i, j-1
end 

% West
i = 1; 
for j = 2:Ny-2
    po = i+(j-1)*(Nx-1);
    A(po, po) = -(2 * gamma2 + 2)/dx2; %i, j
    A(po, po+1) = 1/dx2; %i+1, j
    A(po, po+(Nx-1)) = gamma2/dx2; %i, j+1
    A(po, po-(Nx-1)) = gamma2/dx2; %i, j-1
end 

% East
i = Nx-1; 
for j = 2:Ny-2
    po = i+(j-1)*(Nx-1);
    A(po, po) = -(2 * gamma2 + 2)/dx2; %i, j
    A(po, po-1) = 1/dx2; %i-1, j
    A(po, po+(Nx-1)) = gamma2/dx2; %i, j+1
    A(po, po-(Nx-1)) = gamma2/dx2; %i, j-1
end 

% South West
i = 1; j = 1;
po = i+(j-1)*(Nx-1);
A(po, po) = -(2 * gamma2 + 2)/dx2; %i, j
A(po, po+1) = 1/dx2; %i+1, j
A(po, po+(Nx-1)) = gamma2/dx2; %i, j+1
 
% South East
i = Nx-1; j = 1;
po = i+(j-1)*(Nx-1);
A(po, po) = -(2 * gamma2 + 2)/dx2; %i, j
A(po, po-1) = 1/dx2; %i-1, j
A(po, po+(Nx-1)) = gamma2/dx2; %i, j+1

% North West
i = 1; j = Ny-1;
po = i+(j-1)*(Nx-1);
A(po, po) = -(2 * gamma2 + 2)/dx2; %i, j
A(po, po+1) = 1/dx2; %i+1, j
A(po, po-(Nx-1)) = gamma2/dx2; %i, j-1

% North East
i = Nx-1; j = Ny-1;
po = i+(j-1)*(Nx-1);
A(po, po) = -(2 * gamma2 + 2)/dx2; 
A(po, po-1) = 1/dx2; 
A(po, po-(Nx-1)) = gamma2/dx2; 

end