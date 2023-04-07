% fit_GevaSkinner_model.m
%    fit_GevaSkinner_model Finds the parameters that will fit with the
%    experimentally deconvolved conformational distribution.
%    See PNAS, 2007, 104, 46, 18055 
%    
%    Example
%       parm = fit_GevaSkinner_model(k0,k1,tau,sigma,a,b,x0,y0)
%          k0, rate constant for escaping conformational state 0
%          k1, rate constant for escaping conformational state 1
%          tau, measurement time scale
%          sigma, the standard deviation of molecule-to-molecule variation
%          a, the position of state 1, in unit of angstrom
%          b, the position of state 0, in unit of angstrom
%          x0, experimental x axis
%          y0, experimental y axis
%
%    Version 20061030: sigma now is a vector to allow variations in each
%          time resolution.
%
function [parm,exitf] = fit_GevaSkinner_model_global(k0,k1,tau,sigma,a,b,x0,y0)
  
  parm0 = [k0; k1; a; b; sigma];
  
  [parm,val,exitf] = fminsearch(@(a0) model_sq(a0,x0,y0,tau), parm0,optimset('maxiter',1e5,'maxfunevals',1e5));
return

function sq = model_sq(a0,x0,y0,tau)
  persistent N_iter
  
  if isempty(N_iter)
    N_iter = 0;
  end
  N_iter = N_iter + 1;

  
  k0 = a0(1);
  k1 = a0(2);
  a = a0(3);
  b = a0(4);
  sigma = a0(5:length(tau)+4);


  fprintf(1,'Iteration %5i\n', N_iter);
  fprintf(1,'  k0    = %8.3f\n', k0);
  fprintf(1,'  k1    = %8.3f\n', k1);
  fprintf(1,'  a     = %8.3f\n', a);
  fprintf(1,'  b     = %8.3f\n\n', b);
  for k=1:length(tau),
    fprintf(1,'  sigma(%i) = %8.3f\n', k, sigma(k)*a);
  end
  cd
  s2 = 2.*sigma.*sigma;

  sq = 0;
  for k=1:length(tau),

    dist = GevaSkinner(k0,k1,tau(k)); 

    dx = dist(2,1)-dist(1,1);
    nl = floor(abs(-2.0/dx)); 
    nr = floor(abs(2.0/dx));  
    xl = dist(1,1) - (nl:-1:1)'.*dx;
    xr = dist(length(dist(:,1)),1) + (1:nr)'.*dx;
    x = [xl; dist(:,1); xr];
    
    g = exp(-x.^2./s2(k))./sqrt(2*pi)./sigma(k);
    gg = conv(dist(:,2),g);
    model = gg(1:length(x));
   
    xm = a.*x + b;
    ym = interp1(xm, model, x0(:,k)); % use interpolation
    ym = ym / (sum(ym)*(x0(2,k)-x0(1,k))); % renormalize
    
    sq = sq + mean((ym-y0(:,k)).^2);

  end
return


