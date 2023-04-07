% GevaSkinner Compute 2-state distribution function with motional narrowing
%    DIST = GevaSkinner(K0,K1,TAU) computes the probability density, 
%    DIST(X), as a function of X.
%    K0, the rate constant for the system to jump out of state 0;
%    K1, the rate constant for the system to jump out of state 1;
%    TAU, the average time scale for each measurement.
%
%    Reference
%       1. E. Geva and J. L. Skinner, Chem. Phys. Lett., 288, 225-229
%          (1998).
%
%
%    Revision History
%       HY-20061103, fixed speed problem and scaling error.
%       HY-20061020, start coding and first working code.
%
function dist = GevaSkinner(k0,k1,tau)
  k = k0 + k1;
  p = k0/k; 
  x = (0:0.025:1)';
  dx = x(2)-x(1);
  nx = length(x);
  dist = zeros(length(x),2);

  p0 = (1-p).*exp(-k0.*tau);
  p1 = p.*exp(-k1.*tau);
  dist(2:nx-1,2) = real(quadv(@(y)integrand(y,p,k.*tau,...
    x(2:nx-1),p0,p1),-25.0,25.0))/2/pi;
  dist(1,2) = p0/dx;
  dist(nx,2) = p1/dx;
  dist(:,1) = x;
end

function out = g2(y,p,ktau,p0,p1)

  phi = sqrt((ktau./2).^2-i.*(p-0.5).*ktau.*y-(y./2).^2);
  alpha = ktau./2-i.*(p-0.5).*y;
  tmp = exp(-(ktau+i.*y)./2).*(cosh(phi)+alpha./phi.*sinh(phi));
  out = tmp - p0 - p1.*exp(-i.*y);
end

function out = integrand(y,p,ktau,x,p0,p1)

  out = exp(i.*x.*y).*g2(y,p,ktau,p0,p1);
end
