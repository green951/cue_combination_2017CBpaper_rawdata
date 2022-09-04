function output=alternation_model_normal_2sides(p1_left,mu1_left,mu2_left,std1_left,std2_left,...
    p1_right,mu1_right,mu2_right,std1_right,std2_right, n1, n2)

% altstd_left=alternation_model_normal(p1_left,mu1_left,mu2_left,std1_left,std2_left);
% altstd_right=alternation_model_normal(p1_right,mu1_right,mu2_right,std1_right,std2_right);

if length(mu1_left)>1  % for 2 dimensional
mu_left=sqrt( (mu1_left(1)-mu2_left(1))^2 + (mu1_left(2)-mu2_left(2))^2);
mu_right=sqrt( (mu1_right(1)-mu2_right(1))^2 + (mu1_right(2)-mu2_right(2))^2);

else   % for 1 dimensional
 mu_left=abs(mu1_left-mu2_left);  
  mu_right=abs(mu1_right-mu2_right);
end

altstd_left=sqrt((1-p1_left)*(mu_left^2+std2_left^2) + p1_left *(0^2+std1_left^2)-((1-p1_left)*mu_left+p1_left *0)^2);
altstd_right=sqrt((1-p1_right)*(mu_right^2+std2_right^2) + p1_left *(0^2+std1_right^2)-((1-p1_right)*mu_right+p1_right *0)^2);

std1=altstd_left;
std2=altstd_right;
output= sqrt((std1^2*n1+std2^2*n2)/(n1+n2));
end




