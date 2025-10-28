function [A, b] = findNewAb(Nx, Ny)
Lx = 1; Ly = 1;
dx = Lx/Nx;
dy = Ly/Ny;
gamma = (dx/dy)^2;
A = zeros((Nx-1)*(Ny-1), (Nx-1)*(Ny-1));
b = zeros((Nx-1)*(Ny-1), 1);
dx2 = dx*dx;
C= 25;
%Creating interior grids without boundary condition
for j = 2:Ny-2 %y-grid
    for i = 2:Nx-2 %x-grid
        po = i+(j-1)*(Nx-1); %Determine the position on the grid
        A(po, po) = -(2*gamma/9 + 1/8); %i, j
        A(po, po-1) = 1/16; %i-1, j
        A(po, po+1) = 1/16; %i+1, j
        A(po, po+(Nx-1)) = gamma/9; %i, j+1
        A(po, po-(Nx-1)) = gamma/9; %i, j-1
        bx = dx*i; %x value at particular point
        by = dy*j; %y value at particular point
        b(po) = dx2 * C * (-13/9) * pi^2 * cos(2*pi*(2*bx+by));
    end
end
%Edge and corner grids with the boundary condition
%South edge boundary
j = 1;
for i = 2:Nx-2
    po = i + (j-1)*(Nx-1);
    A(po, po) = -(2*gamma/9 + 1/8); %i, j
    A(po, po-1) = 1/16; %i-1, j
    A(po, po+1) = 1/16; %i+1, j
    A(po, po+(Nx-1)) = gamma/9; %i, j+1
    %A(po, po-(Nx-1)) = gamma2/9 %i, j-1 (boundary)
    bx = dx*i; %x value at particular point
    by = dy*j; %y value at particular point
    b(po) = dx2 * C * (-13/9) * pi^2 * cos(2*pi*(2*bx+by)) - (gamma/9)* (C * cos(4*pi*bx));
end
%North edge boundary
j = Ny-1;
for i = 2:Nx-2
    po = i + (j-1)*(Nx-1);
    A(po, po) = -(2*gamma/9 + 1/8); %i,j
    A(po, po-1) = 1/16; %i-1, j
    A(po, po+1) = 1/16; %i+1, j
    %A(po, po+(Nx-1)) = gamma2/9; %i,j+1 (boundary)
    A(po, po-(Nx-1)) = gamma/9; %i,j-1
    bx = dx*i; %x value at particular point
    by = dy*j; %y value at particular point
    b(po) = dx2 * C * (-13/9) * pi^2 * cos(2*pi*(2*bx+by)) - (gamma/9) * (C * cos(4*pi*bx));
end
%West edge boundary
i = 1;
for j = 2:Ny-2
    po = i+(j-1)*(Nx-1);
    A(po, po) = -(2*gamma/9 + 1/8); %i,j
    %A(po, po-1) = 1/16; %i-1, j
    A(po, po+1) = 1/16; %i+1, j
    A(po, po+(Nx-1)) = gamma/9; %i, j+1
    A(po, po-(Nx-1)) = gamma/9; %i, j-1
    bx = dx*i; %x value at particular point
    by = dy*j; %y value at particular point
    b(po) = dx2 * C * (-13/9) * pi^2 * cos(2*pi*(2*bx+by)) - 1/16* (C * cos(2*pi*by));
end
%East edge boundary
i = Nx-1;
for j = 2:Ny-2
    po = i+(j-1)*(Nx-1);
    A(po, po) = -(2*gamma/9 + 1/8); %i, j
    A(po, po-1) = 1/16; %i-1, j
    %A(po, po+1) = 1/16; %i+1, j
    A(po, po+(Nx-1)) = gamma/9; %i, j+1
    A(po, po-(Nx-1)) = gamma/9; %i, j-1
    bx = dx*i; %x value at particular point
    by = dy*j; %y value at particular point
    b(po) = dx2 * C * (-13/9) * pi^2 * cos(2*pi*(2*bx+by)) - 1/16* (C * cos(2*pi*by));
end
%Southwest boundary
i = 1; j = 1;
po = i+(j-1)*(Nx-1);
A(po, po) = -(2*gamma/9 + 1/8); %i,j
%A(po, po-1) = 1/16; %i-1, j
A(po, po+1) = 1/16; %i+1,j
A(po, po+(Nx-1)) = gamma/9; %i, j+1
%A(po, po-(Nx-1)) = gamma2/9; %i, j-1
bx = dx*i; %x value at particular point
by = dy*j; %y value at particular point
b(po) = dx2 * C * (-13/9) * pi^2 * cos(2*pi*(2*bx+by)) -  (gamma/9)* (C * cos(4*pi*bx)) - 1/16* (C * cos(2*pi*by));
%Southeast boundary
i = Nx-1; j = 1;
po = i+(j-1)*(Nx-1);
A(po, po) = -(2*gamma/9 + 1/8); %i, j
A(po, po-1) = 1/16; %i-1, j
%A(po, po+1) = 1/16; %i+1, j
A(po, po+(Nx-1)) = gamma/9; %i, j+1
%A(po, po-(Nx-1)) = gamma2/9; %i, j-1
bx = dx*i; %x value at particular point
by = dy*j; %y value at particular point
b(po) = dx2 * C * (-13/9) * pi^2 * cos(2*pi*(2*bx+by)) - (gamma/9)* (C * cos(4*pi*bx)) - 1/16* (C * cos(2*pi*by));
%Northwest boundary
i = 1; j = Ny-1;
po = i+(j-1)*(Nx-1);
A(po, po) = -(2*gamma/9 + 1/8); %i, j
%A(po ,po-1) = 1/16; %i-1, j
A(po, po+1) = 1/16; %i+1, j
%A(po, po+(Nx-1)) = gamma2/9; %i, j+1
A(po, po-(Nx-1)) = gamma/9; %i, j-1
bx = dx*i; %x value at particular point
by = dy*j; %y value at particular point
b(po) = dx2 * C * (-13/9) * pi^2 * cos(2*pi*(2*bx+by)) - (gamma/9)* (C * cos(4*pi*bx)) - 1/16*(C * cos(2*pi*by));
%Northeast boundary
i = Nx-1; j = Ny-1;
po = i+(j-1)*(Nx-1);
A(po, po) = -(2*gamma/9 + 1/8); %i, j
A(po, po-1) = 1/16; %i-1, j
%A(po, po+1) = 1/16; %i+1, j
%A(po, po+(Nx-1)) = gamma2/9; %i, j+1
A(po, po-(Nx-1)) = gamma/9; %i, j-1
bx = dx*i; %x value at particular point
by = dy*j; %y value at particular point
b(po) = dx2 * C * (-13/9) * pi^2 * cos(2*pi*(2*bx+by)) - (gamma/9)* (C * cos(4*pi*bx)) - 1/16* (C * cos(2*pi*by));
end