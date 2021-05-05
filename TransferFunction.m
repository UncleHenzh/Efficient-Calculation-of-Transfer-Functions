function [e,S,n] = TransferFunction(A,B,C,D,dt)

% The algorithm presented here follows Appendix D of
% Seem, J. E., et al. "Transfer functions for efficient calculation of 
% multidimensional transient heat transfer." (1989): 5-12.

n=size(A);
n=n(1);
m=size(D);
m=m(2);
phi=expm(A*dt);
gamma1=inv(A)*(phi-eye(n))*B;
gamma2=inv(A)*(gamma1/dt-B);

R0=eye(n);
S=zeros(n+1,m);
S0=C*R0*gamma2+D;
S(1,:)=S0;

Rnew=eye(n);
for j=1:n-1
    Rold=Rnew;
    e(j)=-trace(phi*Rold)/j;
    Rnew=phi*Rold+e(j)*eye(n);
    S(j+1,:)=C*(Rold*(gamma1-gamma2)+Rnew*gamma2)+e(j)*D;
end

e(n)=-trace(phi*Rnew)/n;
Sn=C*Rnew*(gamma1-gamma2)+e(n)*D;
S(n+1,:)=Sn;
end

