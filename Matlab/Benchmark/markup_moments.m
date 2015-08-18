

mu = mu(phi>0);
mu = mu(~isnan(mu));
    
muj = muj(~isnan(muj));



fprintf('\n');
display('Markup Moments');
fprintf('\n');  
fprintf('aggregate markup  = %7.3f \n', P*A);
fprintf('\n');
display('* unconditional markup distribution *');
fprintf('\n');
fprintf('mean markup       = %7.3f \n', mean(mu));
fprintf('\n');
fprintf('markup p50        = %7.3f \n', prctile(mu,50));
fprintf('markup p75        = %7.3f \n', prctile(mu,75));
fprintf('markup p90        = %7.3f \n', prctile(mu,90));
fprintf('markup p95        = %7.3f \n', prctile(mu,95));
fprintf('markup p99        = %7.3f \n', prctile(mu,99));
fprintf('\n');
fprintf('s.d. log markup   = %7.3f \n', std(log(mu)));
fprintf('log p95/p50       = %7.3f \n', log(prctile(mu,95)/prctile(mu,50)));
fprintf('\n');


fprintf('\n');
display('* sectoral markup distribution *');
fprintf('\n');
fprintf('sectoral mean     = %7.3f \n', mean(muj));
fprintf('sectoral s.d. log = %7.3f \n', std(log(muj)));
fprintf('\n');
fprintf('sectoral p50      = %7.3f \n', prctile(muj,50));
fprintf('sectoral p75      = %7.3f \n', prctile(muj,75));
fprintf('sectoral p90      = %7.3f \n', prctile(muj,90));
fprintf('sectoral p95      = %7.3f \n', prctile(muj,95));
fprintf('sectoral p99      = %7.3f \n', prctile(muj,99));
fprintf('\n');
fprintf('sectoral s.d. log = %7.3f \n', std(log(muj)));
fprintf('log p95/p50       = %7.3f \n', log(prctile(muj,95)/prctile(muj,50)));
fprintf('\n');




