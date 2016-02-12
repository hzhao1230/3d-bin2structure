function r = exprndBounded(mu, sizeOut, r1, r2)

minE = exp(-r1/mu); 
maxE = exp(-r2/mu);

randBounded = minE + (maxE-minE).*rand(sizeOut);
r = -mu .* log(randBounded);