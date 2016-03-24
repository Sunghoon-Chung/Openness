function sigma = trade_elasticity(aa, Y, P, omega, yy, parameters , step, impshare)

gamma = parameters(1);
theta = parameters(2);
N     = parameters(3);
S     = parameters(4);
FD    = parameters(5);
FX    = parameters(6);


omegaold = omega;
yyold    = yy; 

for it = 1 : 200
   
omega  = omegaold;
yy     = yyold;

ee = (1/gamma*(1-omega) + 1/theta*omega).^(-1); %% cournot
%ee  = (gamma*(1-omega)+theta*omega);             %% bertrand

pp  = ee./(ee-1)./aa;

profit = (pp - 1./aa).*yy - [FD*ones(N,S); FX*ones(N,S)];

if it < 50    
  
phi  = profit >= 0;   % don't iterate on phi further, discreteness gives jumps, no convergence

end

Pj              = mean(2*phi.*pp.^(1-gamma)).^(1/(1-gamma));
Pj(isinf(Pj))   = 10^16;

omega = 1/N*(pp./repmat(Pj,2*N,1)).^(1-gamma).*phi;

yy    = (pp./repmat(Pj,2*N,1)).^(-gamma).*(repmat(Pj,2*N,1)./P).^(-theta).*Y;

omega = real(omega);
yy    = real(yy); 


error = norm(omega - omegaold);

   %fprintf('%4i %6.2e  \n',[it, error]);

if error < 1e-7, break, end

omegaold     = 1/10*omega +  9/10*omegaold;
yyold        = 1/10*yy    +  9/10*yyold;

end

yy     = yy.*(phi>0);

pH     = pp(1:N,:); 
pF     = pp(N+1:2*N,:); 
yH     = yy(1:N,:); 
yF     = yy(N+1:2*N,:);

impsharenew   = sum(pF(:).*yF(:))/sum(pH(:).*yH(:) + pF(:).*yF(:)); % can't use the sum of omegas since P const


sigma = -log(impsharenew./(1-impsharenew)/(impshare/(1-impshare)))/log(step);
