% domestic and import/foreign markups

mu_domestic = muH;
mu_import   = muF;

mudom = (mean(mean(muH.^(-1).*pH.*yH))/mean(mean(pH.*yH)))^(-1);
mufor = (mean(mean((1+tau).*muF.^(-1).*pF.*yF))/mean(mean((1+tau)*pF.*yF)))^(-1);
domshare = mean(mean(pH.*yH))/(P*Y);

if domshare>0.999,
    
    mucheck = mudom; %% autarky, mufor not defined
    
else

    mucheck = (mudom^(-1)*domshare + mufor^(-1)*(1-domshare))^(-1); %% should be same as aggregate markup

end
    
muagg   = P*A;

mu_domestic = mu_domestic(phiH>0);
mu_domestic = mu_domestic(~isnan(mu_domestic));

mu = mu_domestic;




display('* domestic markup distribution *');
fprintf('\n');  
fprintf('aggregate markup  = %7.3f \n', mudom);
fprintf('\n');
fprintf('mean markup       = %7.3f \n', 1./mean(1./mu));
fprintf('\n');
fprintf('markup p10        = %7.3f \n', prctile(mu,10));
fprintf('markup p25        = %7.3f \n', prctile(mu,25));
fprintf('markup p50        = %7.3f \n', prctile(mu,50));
fprintf('markup p75        = %7.3f \n', prctile(mu,75));
fprintf('markup p90        = %7.3f \n', prctile(mu,90));
fprintf('markup p95        = %7.3f \n', prctile(mu,95));
fprintf('markup p99        = %7.3f \n', prctile(mu,99));
fprintf('\n');
fprintf('s.d. log markup   = %7.3f \n', std(log(mu)));
fprintf('log p95/p50       = %7.3f \n', log(prctile(mu,95)/prctile(mu,50)));
fprintf('\n');

mu_import = mu_import(phiF>0);
mu_import = mu_import(~isnan(mu_import));

mu = mu_import;

display('* import markup distribution *');
fprintf('\n');  
fprintf('aggregate markup  = %7.3f \n', mufor);
fprintf('\n');
fprintf('mean markup       = %7.3f \n', 1./mean(1./mu));
fprintf('\n');
fprintf('markup p10        = %7.3f \n', prctile(mu,10));
fprintf('markup p25        = %7.3f \n', prctile(mu,25));
fprintf('markup p50        = %7.3f \n', prctile(mu,50));
fprintf('markup p75        = %7.3f \n', prctile(mu,75));
fprintf('markup p90        = %7.3f \n', prctile(mu,90));
fprintf('markup p95        = %7.3f \n', prctile(mu,95));
fprintf('markup p99        = %7.3f \n', prctile(mu,99));
fprintf('\n');
fprintf('s.d. log markup   = %7.3f \n', std(log(mu)));
fprintf('log p95/p50       = %7.3f \n', log(prctile(mu,95)/prctile(mu,50)));
fprintf('\n');