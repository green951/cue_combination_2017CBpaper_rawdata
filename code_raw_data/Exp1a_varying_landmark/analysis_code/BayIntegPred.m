% predictions of bayesian integration model for von mises distribution

function [mu_comb, leta_comb,std_comb]=BayIntegPred(mu1,mu2,leta1,leta2)

mu_comb=atan2( (leta1*sin(mu1)+leta2*sin(mu2)),(leta1*cos(mu1)+leta2*cos(mu2)));

leta_comb=sqrt( leta1^2+leta2^2+2*leta1*leta2*cos(mu1-mu2));

std_comb=sqrt( 1-besseli(1,leta_comb) /besseli(0,leta_comb));

end
