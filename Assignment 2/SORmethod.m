function [u0]=SORmethod(A,b,u0,resobj)

omega=1.5;

% SOR
L=tril(A,-1); %Lower triangular matrix
D=diag(diag(A)); %Diagonal matrix
U=A-L-D; %Upper triangular matrix

k=0;

while 1
    u1=(D+omega*L)\(omega*b-(omega*U+(omega-1)*D)*u0);  % SOR method
    res=norm(b-A*u1);
    if res<resobj
        break
    end   
    
    u0=u1;
    k=k+1;
end