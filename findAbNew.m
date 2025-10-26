function [A,b] = findAbNew(Nx,Ny)
    % C value
    C = mean([26,025,024]);
    Lx = 1; Ly = 1;
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
    % find b
    % boundary values
    B = zeros(Ny+1,Nx+1);
    % bottom
    for i = 1:Ny+1
        x = (i-1)*dx;
        y = (i-1)*dy;
        B(end,i) = C*cos(4*pi*y);
    end
    % top
    for i = 1:Ny+1
        x = (i-1)*dx;
        y = (i-1)*dy;
        B(1,i) = C*cos(4*pi*y);
    end
    % left
    for i = 1:Ny+1
        x = (i-1)*dx;
        y = (i-1)*dy;
        B(i,1) = C*cos(2*pi*x);
    end
    % right
    for i = 1:Ny+1
        x = (i-1)*dx;
        y = (i-1)*dy;
        B(i,end) = C*cos(2*pi*x);
    end
    
    % g values and finalize b
    idx = 0;
    for j = Ny:-1:2
        for i =2:Nx
            idx = idx+1;
            b(idx,1) = B(j-1,i)+B(j,i-1)+B(j,i)+B(j+1,i)+B(j,i+1);
            x = (i-1)*dx;
            y = (Ny-j+1)*dy;
            g = -C*(13/9)*pi^2*cos(2*pi*(2*x+y));
            b(idx,1) = 16*g*dx^2 -b(idx,1);
        end
    end
end
