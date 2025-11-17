function [LHS,RHS] = assembleLHSRHS(Nx, Ny, stmfunc,vort, dx, dy ,dt,Re, t)

U_south = 0;
U_north = sin(2*pi*5*t);
U_west = 0;
U_east = 0;

dx2 = dx*dx; dy2 = dy*dy;
gamma = dx/dt; gamma2 = dx2/dy2;
C = dx2/dt;

%At timestep n
vort = 4* C* vort(:);
stmfunc = stmfunc(:);

RHS = zeros((Nx-1)*(Ny-1),1);
LHS = zeros((Nx-1)*(Ny-1), (Nx-1)*(Ny-1)); %Vorticity coefficient at timestep n+1

%Interior points
for j = 2:Ny-2
    for i = 2:Nx-2
        pos = i + (j-1)*(Nx-1);

        LHS(pos,pos) = 4 * C + 8/Re*(1+gamma2); %Coefficient for i,j
        LHS(pos, pos+1) =    gamma * stmfunc(pos+(Nx-1)) - gamma * stmfunc(pos-(Nx-1)) - 4/Re ; %Coefficient for i+1,j
        LHS(pos, pos-1) = -( gamma * stmfunc(pos+(Nx-1)) - gamma * stmfunc(pos-(Nx-1)) + 4/Re); %Coefficient for i-1,j
        LHS(pos, pos+(Nx-1)) = -( gamma * stmfunc(pos+1) - gamma * stmfunc(pos-1) + 4*gamma2/Re); %Coefficient for i,j+1
        LHS(pos, pos-(Nx-1)) =    gamma * stmfunc(pos+1) - gamma * stmfunc(pos-1) - 4*gamma2/Re ; %Coefficient for i,j-1
        RHS(pos) = vort(pos);
    end
end 

%South
j = 1;
for i = 2:Nx-2
    pos = i + (j-1)*(Nx-1);
    vortsouthbc = (0 - stmfunc(pos) + U_south*dy)/(0.5*dy2);
    LHS(pos,pos) = 4 * C + 8/Re*(1+gamma2); %Coefficient for i,j
    LHS(pos, pos+1) =    gamma * stmfunc(pos+(Nx-1)) - gamma * 0              - 4/Re ; %Coefficient for i+1,j
    LHS(pos, pos-1) = -( gamma * stmfunc(pos+(Nx-1)) - gamma * 0              + 4/Re); %Coefficient for i-1,j
    LHS(pos, pos+(Nx-1)) = -( gamma * stmfunc(pos+1) - gamma * stmfunc(pos-1) + 4*gamma2/Re); %Coefficient for i,j+1
    %LHS(pos, pos-(Nx-1)) =    gamma * stmfunc(pos+1) - gamma * stmfunc(pos-1) - 4*gamma2/Re ; %Coefficient for i,j-1
    RHS(pos) = vort(pos) - ( (gamma * stmfunc(pos+1) - gamma * stmfunc(pos-1) - 4*gamma2/Re) * vortsouthbc );

end

%North
j = Ny-1;
for i = 2:Nx-2
    pos = i + (j-1)*(Nx-1);
    vortnorthbc = (0 - stmfunc(pos) - U_north*dy)/(0.5*dy2);
    LHS(pos,pos) = 4 * C + 8/Re*(1+gamma2); %Coefficient for i,j
    LHS(pos, pos+1) =    gamma * 0 -              gamma * stmfunc(pos-(Nx-1)) - 4/Re ; %Coefficient for i+1,j
    LHS(pos, pos-1) = -( gamma * 0 -              gamma * stmfunc(pos-(Nx-1)) + 4/Re); %Coefficient for i-1,j
    %LHS(pos, pos+(Nx-1)) = -( gamma * stmfunc(pos+1) - gamma * stmfunc(pos-1) + 4*gamma2/Re); %Coefficient for i,j+1
    LHS(pos, pos-(Nx-1)) =    gamma * stmfunc(pos+1) - gamma * stmfunc(pos-1) - 4*gamma2/Re ; %Coefficient for i,j-1
    RHS(pos) = vort(pos) + ( (gamma * stmfunc(pos+1) - gamma * stmfunc(pos-1) + 4*gamma2/Re) * vortnorthbc );
 
end 

%West
i = 1; 
for j = 2:Ny-2
    pos = i + (j-1)*(Nx-1);
    vortwestbc = (0 - stmfunc(pos) - U_west*dx)/(0.5*dx2);
    LHS(pos,pos) = 4 * C + 8/Re*(1+gamma2); %Coefficient for i,j
    LHS(pos, pos+1) =    gamma * stmfunc(pos+(Nx-1)) - gamma * stmfunc(pos-(Nx-1)) - 4/Re ; %Coefficient for i+1,j
    %LHS(pos, pos-1) = -( gamma * stmfunc(pos+(Nx-1)) - gamma * stmfunc(pos-(Nx-1)) + 4/Re); %Coefficient for i-1,j
    LHS(pos, pos+(Nx-1)) = -( gamma * stmfunc(pos+1) - gamma * 0 +             4*gamma2/Re); %Coefficient for i,j+1
    LHS(pos, pos-(Nx-1)) =    gamma * stmfunc(pos+1) - gamma * 0 -             4*gamma2/Re ; %Coefficient for i,j-1
    RHS(pos) = vort(pos) + ( (gamma * stmfunc(pos+(Nx-1)) - gamma * stmfunc(pos-(Nx-1)) +4/Re) * vortwestbc );

end 

%East 
i = Nx-1;
for j = 2:Ny-2
    pos = i + (j-1)*(Nx-1);
    vorteastbc = (0 - stmfunc(pos) + U_east*dx)/(0.5*dx2);
    LHS(pos,pos) = 4 * C + 8/Re*(1+gamma2); %Coefficient for i,j
    %LHS(pos, pos+1) =    gamma * stmfunc(pos+(Nx-1)) - gamma * stmfunc(pos-(Nx-1)) - 4/Re ; %Coefficient for i+1,j
    LHS(pos, pos-1) = -( gamma * stmfunc(pos+(Nx-1)) - gamma * stmfunc(pos-(Nx-1)) + 4/Re); %Coefficient for i-1,j
    LHS(pos, pos+(Nx-1)) = -( gamma * 0 -             gamma * stmfunc(pos-1) + 4*gamma2/Re); %Coefficient for i,j+1
    LHS(pos, pos-(Nx-1)) =    gamma * 0 -             gamma * stmfunc(pos-1) - 4*gamma2/Re ; %Coefficient for i,j-1
    RHS(pos) = vort(pos) + ( -(gamma * stmfunc(pos+(Nx-1)) - gamma * stmfunc(pos-(Nx-1)) - 4/Re) * vorteastbc );
end 

%South West
i = 1; j = 1;
pos = i + (j-1)*(Nx-1);
vortsouthbc = (0 - stmfunc(pos) + U_south*dy)/(0.5*dy2);
vortwestbc  = (0 - stmfunc(pos) - U_west*dx)/(0.5*dx2);
LHS(pos,pos) = 4 * C + 8/Re*(1+gamma2); %Coefficient for i,j
LHS(pos, pos+1) =    gamma * stmfunc(pos+(Nx-1)) - gamma * 0 - 4/Re ; %Coefficient for i+1,j
%LHS(pos, pos-1) = -( gamma * stmfunc(pos+(Nx-1)) - gamma * 0 + 4/Re); %Coefficient for i-1,j
LHS(pos, pos+(Nx-1)) = -( gamma * stmfunc(pos+1) - gamma * 0 + 4*gamma2/Re); %Coefficient for i,j+1
%LHS(pos, pos-(Nx-1)) =    gamma * stmfunc(pos+1) - gamma * 0 - 4*gamma2/Re ; %Coefficient for i,j-1
RHS(pos) = vort(pos) + ( (gamma * stmfunc(pos+(Nx-1)) - gamma * 0 + 4/Re) * vortwestbc ) ...
                     - ( (gamma * stmfunc(pos+1) - gamma * 0 - 4*gamma2/Re) * vortsouthbc );

%South East
i = Nx-1; j = 1;
pos = i + (j-1)*(Nx-1);
vortsouthbc  = (0 - stmfunc(pos) + U_south*dy)/(0.5*dy2);
vorteastbc   = (0 - stmfunc(pos) + U_east*dx)/(0.5*dx2);
LHS(pos,pos) = 4 * C + 8/Re*(1+gamma2); %Coefficient for i,j
%LHS(pos, pos+1) =    gamma * stmfunc(pos+(Nx-1)) - gamma * 0 - 4/Re ; %Coefficient for i+1,j
LHS(pos, pos-1) = -( gamma * stmfunc(pos+(Nx-1)) - gamma * 0 + 4/Re); %Coefficient for i-1,j
LHS(pos, pos+(Nx-1)) = -( gamma * 0 - gamma * stmfunc(pos-1) + 4*gamma2/Re); %Coefficient for i,j+1
%LHS(pos, pos-(Nx-1)) =    gamma * 0 - gamma * stmfunc(pos-1) - 4*gamma2/Re ; %Coefficient for i,j-1
RHS(pos) = vort(pos) - ( (gamma * stmfunc(pos+(Nx-1)) - gamma * 0 - 4/Re) * vorteastbc ) ...
                     - ( (gamma * 0 - gamma * stmfunc(pos-1) - 4*gamma2/Re) * vortsouthbc );

%North West
i = 1; j = Ny-1;
pos = i + (j-1)*(Nx-1);
vortnorthbc = (0 - stmfunc(pos) - U_north*dy)/(0.5*dy2);
vortwestbc  = (0 - stmfunc(pos) - U_west*dx)/(0.5*dx2);
LHS(pos,pos) = 4 * C + 8/Re*(1+gamma2); %Coefficient for i,j
LHS(pos, pos+1) =    gamma * 0 - gamma * stmfunc(pos-(Nx-1)) - 4/Re ; %Coefficient for i+1,j
%LHS(pos, pos-1) = -( gamma * 0 - gamma * stmfunc(pos-(Nx-1)) + 4/Re); %Coefficient for i-1,j
%LHS(pos, pos+(Nx-1)) = -( gamma * stmfunc(pos+1) - gamma * 0 + 4*gamma2/Re); %Coefficient for i,j+1
LHS(pos, pos-(Nx-1)) =    gamma * stmfunc(pos+1) - gamma * 0 - 4*gamma2/Re ; %Coefficient for i,j-1
RHS(pos) = vort(pos) + ( (gamma * 0 - gamma * stmfunc(pos-(Nx-1)) + 4/Re) * vortwestbc ) ...
                     + ( (gamma * stmfunc(pos+1) - gamma * 0 + 4*gamma2/Re) *vortnorthbc);

%North East
i = Nx-1; j = Ny-1;
pos = i + (j-1)*(Nx-1);
vortnorthbc = (0 - stmfunc(pos) - U_north*dy)/(0.5*dy2);
vorteastbc  = (0 - stmfunc(pos) + U_east*dx)/(0.5*dx2);
LHS(pos,pos) = 4 * C + 8/Re*(1+gamma2); %Coefficient for i,j
%LHS(pos, pos+1) =    gamma * 0 - gamma * stmfunc(pos-(Nx-1)) - 4/Re ; %Coefficient for i+1,j
LHS(pos, pos-1) = -( gamma * 0 - gamma * stmfunc(pos-(Nx-1)) + 4/Re); %Coefficient for i-1,j
%LHS(pos, pos+(Nx-1)) = -( gamma * 0 - gamma * stmfunc(pos-1) + 4*gamma2/Re); %Coefficient for i,j+1
LHS(pos, pos-(Nx-1)) =    gamma * 0 - gamma * stmfunc(pos-1) - 4*gamma2/Re ; %Coefficient for i,j-1
RHS(pos)  =  vort(pos) - ( (gamma * 0 - gamma * stmfunc(pos-(Nx-1)) - 4/Re) * vorteastbc ) ...
                       + ( (gamma * 0 - gamma * stmfunc(pos-1) + 4*gamma2/Re) * vortnorthbc );

