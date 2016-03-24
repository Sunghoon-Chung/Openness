rng(0);   % reset random number generator

% Idiosyncratic firm productivity

[ux,uy] = gumbel_copula(N*S,0); %% copula function, independent since kendall's tau=0 

a  = (1-ux).^(-1/xi_x);         % $a = F^{-1}(u_x)$ is the inverse of Pareto CDF, $u_x$, with shape param. $xi_x$
as = (1-uy).^(-1/xi_x);         % a: productivity of domestic firms, as: productivity for foreign firms

clear ux uy; 

nn = geornd(zeta,S,1);          %% sample from geometric distribution with parameter zeta

phi = zeros(N, S);              %% matrix with producer number in each sector


for i = 1:S
    temp = zeros(N,1); 
    temp(1:nn(i)) = 1;          %% generate temp with 1 = # of producers, or 0
    temp = temp(randperm(N));   %% randperm: random permutation from 1 to N ==> reshuffling the position
    phi(:,i) = temp; 
end


a   = reshape(a, N, S);         %% mapping productivities to assigned numbers above (i.e., firms)
as  = reshape(as, N, S);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sector productivity 

mu = [0, 0];
zcovar = xi_z^2 .* [1 rho; rho 1];
ln_sectorz = mvnrnd(mu, zcovar, S);
sectorz = exp(ln_sectorz);
z  = sectorz(:,1);
zs = sectorz(:,2);

a  = bsxfun(@times, a, z');     %% bsxfun = element-by-element binary operation to arrays with singleton expansion enabled: a(NxS) x z'(1xS)
as = bsxfun(@times, as, zs'); 

a  = a.*phi  + eps.*(1-phi);    %% eps = 2.2204e-16, floating-point relative accuracy
as = as.*phi + eps.*(1-phi);


a  = sort(a, 1,'descend');      %% sort each column (row if 2) in descending order
as = sort(as,1,'descend');

% MAKE COUNTRIES PERFECTLY SYMMETRIC: Drop 1/2 observations, switch around

a  = [a, as];
as = [a(:, S+1:end), a(:,1:S)];


S  = S*2;

a  = sort(a, 1, 'descend');
as = sort(as, 1, 'descend');


% initial guess: economy with p = gamma/(gamma-1)*W/a

aa = [a;as/(1+tau)]; 

pp  = gamma/(gamma-1)./aa;                      % price value: 2*NxS (=558x2000) matrix

Pj    = (mean(2*pp.^(1-gamma))).^(1/(1-gamma));  % original

P     = (mean(Pj.^(1-theta))).^(1/(1-theta));    % Aggregate Price: original

phi = ones(2*N, S);

% repmat(Pj,2*N,1): replication matrix of Pj by 2*N times of row and 1 times of column
A    = mean(mean(2*phi.*aa.^(-1).*(pp./repmat(Pj,2*N,1)).^(-gamma).*repmat(Pj/P, 2*N,1).^(-theta))).^(-1); % calculation based on (24)
Y    = A*L;
Yold = Y;

% iterate until convergence on firm decision rules and exit decisions

omegaold  =  1/N*(pp./repmat(Pj,2*N,1)).^(1-gamma); 
yyold      = (pp./repmat(Pj,2*N,1)).^(-gamma).*(repmat(Pj,2*N,1)./P).^(-theta).*Y;


for it = 1 : 200
   
Y      = Yold;  
omega  = omegaold;
yy     = yyold;

%ee = (1/gamma*(1-omega) + 1/theta*omega).^(-1); % Cournot model
ee = (gamma*(1-omega)+theta*omega);           % Bertrand model

pp  = ee./(ee-1)./aa;

profit = (pp - 1./aa).*yy - [FD*ones(N,S); FX*ones(N,S)];

if it < 50    
  
phi  = profit >= 0;   % don't iterate on phi further, discreteness gives jumps, no convergence

end

Pj              = mean(2*phi.*pp.^(1-gamma)).^(1/(1-gamma));
Pj(isinf(Pj))   = 10^16;

P               = mean(Pj.^(1-theta)).^(1/(1-theta));

omega = 1/N*bsxfun(@times, pp.^(1-gamma).*phi, Pj.^(gamma-1));

yy = bsxfun(@times, pp.^(-gamma), Pj.^(gamma-theta).*P.^theta.*Y);

omega = real(omega);            %% real() = real part of complex number
yy    = real(yy); 

A  = mean(mean(2*phi.*aa.^(-1).*yy/Y)).^(-1);   

fcost = phi.*[FD*ones(N,S); FX*ones(N,S)];
Y     = A*(L-2*mean(fcost(:)));

error = max([norm(omega - omegaold); norm(Y-Yold)]);        %% norm() = Euclidean norm


%fprintf('%4i %6.2e \n',[it, error]);


if error < 1e-7, break, end
    
Yold         = 1/2*Y      +  1/2*Yold;
omegaold     = 1/10*omega +  9/10*omegaold;
yyold        = 1/10*yy    +  9/10*yyold;

end

yy  = yy.*(phi>0);
omega = omega.*(phi>0);

ll = yy./aa;

if norm(omega-omegaold)>1e-3
    disp('Warning: omega has not converged')
end

mu  = ee./(ee-1); 

% Naive Trade Elasticity (fix measures of producers)

omegaH = omega(1:N,:); 
omegaF = omega(N+1:2*N,:); 
pH     = pp(1:N,:); 
pF     = pp(N+1:2*N,:);    % this is (1+tau)*pF
yH     = yy(1:N,:); 
yF     = yy(N+1:2*N,:); 
phiH   = phi(1:N,:); 
phiF   = phi(N+1:2*N,:);
muH    = mu(1:N,:); 
muF    = mu(N+1:2*N,:); 
lH     = ll(1:N,:); 
lF     = ll(N+1:2*N,:);


% domestic and foreign markups
mudom = (mean(mean(muH.^(-1).*pH.*yH))/mean(mean(pH.*yH)))^(-1);
mufor = (mean(mean((1+tau).*muF.^(-1).*pF.*yF))/mean(mean((1+tau)*pF.*yF)))^(-1);
domshare = mean(mean(pH.*yH))/(P*Y);

if domshare>0.999,
    
    mucheck = mudom; %% autarky, mufor not defined
    
else

    mucheck = (mudom^(-1)*domshare + mufor^(-1)*(1-domshare))^(-1); %% should be same as aggregate markup

end
    
muagg   = P*A;


lambdaj = sum(omegaH);
sj      = (Pj./P).^(1-theta);
lambda  = mean(lambdaj.*sj);
aweight = mean(sj.*lambdaj/lambda.*(1-lambdaj)/(1-lambda));
naiveArm     = gamma*aweight + theta*(1-aweight);

if abs(lambda-1) < 1e-3   % don't divide by 0
  
    aweight = 1;
    naiveArm = gamma;
    
end

relmarkups  = mean(mean(phiF.*(pF./repmat(Pj,N,1)).^(1-gamma).*(repmat(Pj,N,1)/P).^(1-theta)))/...
              mean(mean(phiH.*(pH./repmat(Pj,N,1)).^(1-gamma).*(repmat(Pj,N,1)/P).^(1-theta)));

impshare    = sum(omegaF)';
import_pen  = impshare; 

% use the fact that everything is symmetric (a(:,j) = as(:,S/2+j), as(:,j)
% = a(:,S/2+j), to back out "export share" of industry j from 
% = import share of industry S/2 + j

pHs     = [pF(:,S/2+1:S), pF(:,1:S/2)];
yHs     = [yF(:,S/2+1:S), yF(:,1:S/2)];
phiHs   = [phiF(:,S/2+1:S), phiF(:,1:S/2)];



omegaHs = [omegaF(:,S/2+1:S), omegaF(:,1:S/2)];

expshare   = sum(omegaHs)';

ttrade     = import_pen + expshare;
intraindex  = mean(sj(ttrade>0)'.*abs(impshare(ttrade>0)-expshare(ttrade>0))./(impshare(ttrade>0) + expshare(ttrade>0)))./mean(sj(ttrade>0)');
     
impshare   = sum(pF(:).*yF(:))/sum(pH(:).*yH(:) + pF(:).*yF(:));

agg_fexporters = sum(omegaH(:)>0 & omegaF(:)>0)/sum(omegaH(:)>0);

% Correct Trade Elasticity (trade costs affect markups and measures of producers)

if comp_trade_elasticity
    
parameters   = [gamma, theta, N, S, FD, FX ];

step   = 1.005;  
taunew = (1+tau)*step - 1;
aanew  = [a; as/(1+taunew)]; 

sigma = trade_elasticity_bertrand(aanew, Y, P, omega, yy, parameters, step, impshare); 

%Arm = log(impsharenew./(1-impsharenew)/(impshare/(1-impshare)))/log(step);

else
    
sigma = 0;
    
end


clear a as
S = S/2;




