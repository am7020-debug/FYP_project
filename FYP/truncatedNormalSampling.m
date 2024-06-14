function output = truncatedNormalSampling(mu, sigma, N, varargin)


checkPositive = @(x) assert(isnumeric(x) && isscalar(x) ...
    && (x > 0),'Wrong input');

% Parse inputs
p = inputParser();
addParameter(p,'k', 1, checkPositive) ;  % Default: k = 1
parse(p,varargin{:});
k = p.Results.k ;

xmean = mu ;
xmin = xmean - k * sigma ;
if xmin < 0
    xmin = 0 ;
end
xmax = xmean + k * sigma ;
pd = makedist('Normal', 'mu', xmean, 'sigma', sigma) ;
t = truncate(pd, xmin, xmax) ;
output = random(t, N, 1);