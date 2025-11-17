function A = findMatrixA(Nx,Ny,Lx,Ly) 

% initialize
I = Nx-1;
J = Ny-1;
IJ = I*J;
A = zeros(IJ);
dx = Lx/Nx; dy = Ly/Ny;
gamma = (4/3)*(dx/dy);
D = -2*(gamma^2+1);
% fill in D
i = IJ;
for j = IJ:-1:1
    A(j,i) = D;
    i = i-1;
end
% fill in gamma
i =IJ-Nx+1;
for j = IJ:-1:1
    A(j,i) = gamma^2;
    i = i-1;
    if i*j == 0
        break
    end
end
i =IJ;
for j = IJ-Nx+1:-1:1
    A(j,i) = gamma^2;
    i = i-1;
    if i*j == 0
        break
    end
end
% fill in 1
for i = 1:IJ-1
    if i/I == floor(i/I)
        A(i,i+1) =0;
        A(i+1,i) =0;
    else 
        A(i,i+1) =1;
        A(i+1,i) =1;
    end
end

end
