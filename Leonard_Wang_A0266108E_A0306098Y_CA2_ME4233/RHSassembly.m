function [RHS] = RHSassembly(Nx, Ny, stmfunc, vort, Re, dx, dy, t)

%The boundary conditions at the boundaries
U_south = zeros(Nx-1, 1);
U_north = ones(Nx-1, 1);
U_west = zeros(Ny-1, 1);
U_east = zeros(Ny-1, 1);
dx2 = dx*dx; dy2 = dy*dy;

nu = 1/Re;

%interior
for j = 2:Ny-2
    for i = 2:Nx-2
        fac1 = -(stmfunc(i,j+1) -  stmfunc(i,j-1))/(2*dy);
        fac2 = (stmfunc(i+1,j) - stmfunc(i-1,j))/(2*dx);

        RHS(i,j) = fac1 * ( ( vort(i+1,j) - vort(i-1,j) )/ (2*dx) ) + ...
                   fac2 * ( ( vort(i,j+1) - vort(i,j-1) )/ (2*dy) ) + ...
                   nu * ( ( vort(i+1,j) - 2* vort(i,j) + vort(i-1,j) )/(dx2) + ...
                          ( vort(i,j+1) - 2* vort(i,j) + vort(i, j-1) )/(dy2) );

    end 
end 

%South
j = 1;
for i = 2:Nx-2
    vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy )/(0.5 * dy2); %Imposing the south vorticity bc
    fac1 = -(stmfunc(i,j+1) - 0             )/(2*dy);
    fac2 =  (stmfunc(i+1,j) - stmfunc(i-1,j))/(2*dx);

    RHS(i,j) = fac1 * ( ( vort(i+1,j) - vort(i-1,j) )/ (2*dx) ) + ...
               fac2 * ( ( vort(i,j+1) - vortsouthbc )/ (2*dy) ) + ...
               nu   * ( ( vort(i+1,j) - 2*vort(i,j) + vort(i-1,j) )/ (dx2) + ...
                      ( vort(i,j+1) - 2*vort(i,j) + vortsouthbc )/ (dy2) );

end

%North
j = Ny-1;
for i = 2:Nx-2
    vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*sin(2*pi*5*t)*dy )/(0.5 * dy2); %Imposing the north vorticity bc 
    fac1 = -(0              - stmfunc(i,j-1))/(2*dy);
    fac2 =  (stmfunc(i+1,j) - stmfunc(i-1,j))/(2*dx);

    RHS(i,j) = fac1 * ( ( vort(i+1,j) - vort(i-1,j) )/ (2*dx)) + ...
               fac2 * ( ( vortnorthbc - vort(i,j-1) )/ (2*dy)) + ...
               nu   * ( ( vort(i+1,j) - 2*vort(i,j) + vort(i-1,j) )/ (dx2) + ...
                        ( vortnorthbc - 2*vort(i,j) + vort(i,j-1) )/ (dy2)); 

end 

%West
i = 1;
for j = 2:Ny-2
    vortwestbc = ( 0 - stmfunc(i,j) - U_west(j)*dx )/(0.5 * dx2); %Imposing the west vorticity bc
    fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/(2*dy); 
    fac2 =  (stmfunc(i+1,j) - 0             )/(2*dx);

    RHS(i,j) = fac1 * ( ( vort(i+1,j) - vortwestbc )/ (2*dx)) + ...
               fac2 * ( ( vort(i,j+1) - vort(i,j-1) )/ (2*dy)) + ...
               nu   * ( ( vort(i+1,j) - 2*vort(i,j) + vortwestbc )/ (dx2) + ...
                        ( vort(i,j+1) - 2*vort(i,j) + vort(i,j-1) )/ (dy2));

end 

%East
i = Nx-1;
for j = 2:Ny-2
    vorteastbc = ( 0 - stmfunc(i,j) + U_east(j)*dx )/(0.5 * dx2); %Imposing the east vorticity bc
    fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/(2*dy);
    fac2 =  (0              - stmfunc(i-1,j))/(2*dx);

    RHS(i,j) = fac1 * ( ( vorteastbc - vort(i-1,j) ) )/ (2*dx) + ...
               fac2 * ( ( vort(i,j+1) - vort(i,j-1) ) )/ (2*dy) + ...
               nu   * ( ( vorteastbc - 2*vort(i,j) + vort(i-1,j) )/ (dx2) + ...
                        ( vort(i,j+1) - 2*vort(i,j) + vort(i,j-1) )/ (dy2));
end 

%South West
i = 1; j=1;
vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/ (0.5 * dy2); %Imposing the south vorticity bc
vortwestbc  = ( 0 - stmfunc(i,j) - U_west(j)*dx )/ (0.5 * dx2); %Imposing the west vorticity bc
    fac1 = -(stmfunc(i,j+1) - 0             )/(2*dy);
    fac2 =  (stmfunc(i+1,j) - 0             )/(2*dx);

    RHS(i,j) = fac1 * ( (vort(i+1,j) - vortwestbc )/ (2*dx)) + ...
               fac2 * ( (vort(i,j+1) - vortsouthbc )/ (2*dy)) + ...
               nu   * ( (vort(i+1,j) - 2*vort(i,j) + vortwestbc )/ (dx2) + ...
                        (vort(i,j+1) - 2*vort(i,j) + vortsouthbc )/ (dy2));


%South East
i = Nx-1; j=1;
vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/ (0.5 * dy2); %Imposing the south vorticity bc
vorteastbc  = ( 0 - stmfunc(i,j) + U_east(j)*dx )/ (0.5 * dx2); %Imposing the east vorticity bc
    fac1 = -(stmfunc(i,j+1) - 0             )/(2*dy);
    fac2 =  (0              - stmfunc(i-1,j))/(2*dx);

    RHS(i,j) = fac1 * ( ( vorteastbc - vort(i-1,j) )/ (2*dx)) + ...
               fac2 * ( ( vort(i,j+1) - vortsouthbc )/ (2*dy)) + ...
               nu   * ( ( vorteastbc - 2*vort(i,j) + vort(i-1,j) )/ (dx2) + ...
                        ( vort(i,j+1) - 2*vort(i,j) + vortsouthbc)/ (dy2));


%North West
i = 1; j = Ny-1;
vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*sin(2*pi*5*t)*dy)/ (0.5 * dy2); %Imposing the north vorticity bc
vortwestbc  = ( 0 - stmfunc(i,j) - U_west(j)*dx )/ (0.5 * dx2); %Imposing the west vorticity bc
    fac1 = -(0               - stmfunc(i,j-1))/(2*dy);
    fac2 =  (stmfunc(i+1,j) - 0              )/(2*dx);

    RHS(i,j) = fac1 * ( ( vort(i+1,j) - vortwestbc )/ (2*dx)) + ...
               fac2 * ( ( vortnorthbc - vort(i,j-1) )/ (2*dy)) + ...
               nu   * ( ( vort(i+1,j) - 2*vort(i,j) + vortwestbc )/ (dx2) + ...
                        ( vortnorthbc - 2*vort(i,j) + vort(i,j-1) )/ (dy2));

   
%North East
i = Nx-1; j = Ny-1;
vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*sin(2*pi*5*t)*dy)/ (0.5 * dy2); %Imposing the north vorticity bc
vorteastbc  = ( 0 - stmfunc(i,j) + U_east(j)*dx )/ (0.5 * dx2); %Imposing the east vorticity bc
    fac1 = -(0               - stmfunc(i,j-1))/(2*dy);
    fac2 =  (0               - stmfunc(i-1,j))/(2*dx);

    RHS(i,j) = fac1 * ( ( vorteastbc - vort(i-1,j) )/ (2*dx)) + ...
               fac2 * ( ( vortnorthbc - vort(i,j-1) )/ (2*dy)) + ...
               nu   * ( ( vorteastbc - 2*vort(i,j) + vort(i-1,j) )/ (dx2) + ...
                        ( vortnorthbc - 2*vort(i,j) + vort(i,j-1) )/ (dy2));



