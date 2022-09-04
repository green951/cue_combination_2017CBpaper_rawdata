function output=alternation_model_normal_2sidesTogether(p1,mu1,mu2,std1,std2)


if length(mu1)>1  % for 2 dimensional
mu=sqrt( (mu1(1)-mu2(1))^2 + (mu1(2)-mu2(2))^2);

else
 mu=abs(mu1-mu2);  
  
end

altstd=sqrt((1-p1)*(mu^2+std2^2) + p1 *(0^2+std1^2)-((1-p1)*mu+p1 *0)^2);


output=altstd;

end




