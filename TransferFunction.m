function [e,S,n] = TransferFunction(A,B,C,D,dt)

% You running this script/function means you will not blame the author(s) 
% if this breaks your stuff. This script/function is provided AS IS without 
% warranty of any kind. Author(s) disclaim all implied warranties 
% including, without limitation, any implied warranties of merchantability 
% or of fitness for a particular purpose. The entire risk arising out of 
% the use or performance of the sample scripts and documentation remains 
% with you. In no event shall author(s) be held liable for any damages 
% whatsoever (including, without limitation, damages for loss of business 
% profits, business interruption, loss of business information, or other 
% pecuniary loss) arising out of the use of or inability to use the script
% or documentation. Neither this script/function, nor any part of it other 
% than those parts that are explicitly copied from others, may be 
% republished without author(s) express written permission. Author(s) 
% retain the right to alter this disclaimer at any time.

% The codes presented here is for the algorithm in Appendix D of
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

