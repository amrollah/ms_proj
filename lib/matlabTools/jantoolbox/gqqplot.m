function gqqplot (x,dist,marker,N)
%GQQPLOT Display a generic Q-Q plot of sample vs. distribution
%	GQQPLOT(X,DIST) makes an plot of the quantiles of the data 
%	set X versus the inverse of the cdf of a distribution specified in 
%	DIST, evaluated at probabilities equal to the quantiles of X. The 
%	parameters of the distribution are calculated from X.
%
%	GQQPLOT(X,DIST,MARKER) uses the marker string specified
%
%	The purpose of the QQ plot is to determine whether the sample in X 
%	is drawn from a given distribution. If it is, the plot will be linear. 
%	In the case of the Binomial distribution, an additional parameter is 
%	needed: N, number of trials, e.g. GQQPLOT(X,'binom',MARKER,N). 
%
%	The name of the distribution can be (case insensitive):
%
%	'norm'     or 'normal'
%	'gam'      or 'gamma'
%	'logn'     or 'Lognormal',
%	'exp'      or 'Exponential',
%	'wbl'      or 'Weibull'.
%	'beta'
%	'bin'      or 'Binomial'
%	'ev'       or 'Extreme Value',
%	'gev'      or 'Generalized Extreme Value',
%	'gp'       or 'Generalized Pareto',
%	'nbin'     or 'Negative Binomial',
%	'poiss'    or 'poisson'
%	'unif'     or 'uniform'
%	'rayl'     or 'rayleigh'
%
%	EXAMPLE:
%	x= gamrnd(2,0.5,1000,1);
%	gqqplot(x,'normal') %Bad fit
%	gqqplot(x,'gamma') %Good fit
%
%	Feb. 2008
%	Requires: Statistics Toolbox

if nargin<3, marker=[]; end
if nargin<4, N=1; end

if isempty(marker), marker='.'; end

x= sort(x(:))';
qd = [];
np= 1000;
if length(x)>np
  p = ((1:np)-0.5)/np;
else
  p = ((1:length(x))-0.5)/length(x);
  qd = x;
end

has_stats = license('test','statistics_toolbox');
if has_stats
  if isempty(qd), qd= quantile(x,p); end
  tit= 'QQ Plot of Sample Data versus ';
else
  if isempty(qd), qd = interp1(((1:length(x))-0.5)/length(x),x,p); end
  tit= 'QQ plot of sample data versus ';
end

switch lower(dist)
	case {'norm','normal'}
    if has_stats
		  [mu,sigma]= normfit(x);
		  y= icdf('normal',p,mu,sigma);
    else
      mu = mean(x);
      sigma = std(x);
      y = (-sqrt(2)*sigma).*erfcinv(2*p) + mu;
    end
		tit= [tit,'Normal'];

	case {'gam','gamma'}
		pgam= gamfit(x);
		y= icdf('gamma',p,pgam(1),pgam(2));
		tit= [tit,'Gamma'];

	case {'logn','lognormal'}
		plogn= lognfit(x);
		y= icdf('logn',p,plogn(1),plogn(2));
		tit= [tit, 'Lognormal'];

	case {'exp','exponential'}
    if has_stats
  		mu= expfit(x);
	  	y= icdf('exp',p,mu);
    else
      mu = mean(x);
      y =  -mu*log(1-p);
    end
		tit= [tit,'Exponential'];

	case {'wbl','weibull'}
		pwbl= wblfit(x);
		y= icdf('wbl',p,pwbl(1),pwbl(2));
		tit= [tit, 'Weibull'];

	case 'beta'
		betap= betafit(x);
		y= icdf('beta',p,betap(1),betap(2));
		tit= [tit, 'Beta'];

	case {'bino','binomial'}
		pbin= binofit(sum(x),N*length(x));
		y= binoinv(p,N,pbin)
		tit= [tit, 'Binomial'];

	case {'ev','extreme value'}
		pev= evfit(x);
		y= icdf('ev',p,pev(1),pev(2));
		tit= [tit, 'Extreme Value'];

	case {'gev','generalized extreme value'} 
		pgev= gevfit(x);
		y= icdf('gev',p,pgev(1),pgev(2),pgev(3));
		tit= [tit, 'Generalized Extreme Value'];

	case {'gp','generalized pareto'}
		warning ('Location parameter (p) must be 0')
		[xi,sigma]= gp_fit(x);
		y = gpwinv(p,[xi 1 sigma]);
		tit= [tit, 'Generalized Pareto'];
    
  case {'gpw'}
		pgpw = gpwfit(x);
		y = gpwinv(p,pgpw);
		tit = [tit 'GPW (' num2str(pgpw(1),'%.2f') ', ' num2str(pgpw(2),'%.2f') ', ' num2str(pgpw(3),'%.2f') ')'];

	case {'nbin','negative binomial'}
		pnbin= nbinfit(x);
		y= icdf('nbin',p,pnbin(1),pnbin(2));
		tit= [tit, 'Negative Binomial'];

	case {'poiss','poisson'}
		r= poissfit(x);
		y= icdf('poisson',p,r);
		tit= [tit, 'Poisson'];

	case {'unif','uniform'}
		y= unifinv(p,min(x),max(x));
		tit= [tit, 'Uniform'];

	case {'rayl','rayleigh'}
		mu= raylfit(x);
		y= icdf('rayleigh',p,mu);
		tit= [tit,'Rayleigh'];

	otherwise 
		error ('Unrecognized distribution name')
end

plot(y,qd,marker,[min(y) max(y)],[min(y) max(y)],'k--');
xlabel('Theoretical Quantiles'); 
ylabel('Quantiles of Input Sample');
title(tit);
