function output=alternation_model_normal(p1,mu1,mu2,std1,std2)

Ntrials=500;
data=zeros(Ntrials,3);
if length(mu1)>1  % for 2 dimensional
mu=sqrt( (mu1(1)-mu2(1))^2 + (mu1(2)-mu2(2))^2);
else
 mu=abs(mu1-mu2);  
end
for i=1:Ntrials
% distribParams1(1) = mu1;
%         distribParams1(2) = leta1;
        data(i,1)=randn(1,1)*std1+0;
        %data(i,1) = randraw('vonmises', distribParams1, 1);
        
% distribParams2(1) = mu2;
%         distribParams2(2) = leta2;
        data(i,2)=randn(1,1)*std2+mu;
       % data(i,2) = randraw('vonmises', distribParams2, 1);
        
        if rand<p1
            data(i,3)=data(i,1);
        else
            data(i,3)=data(i,2);
        end
             
        
        
        
end
output=std(data(:,3));
%output=circ_std(data(:,3));
leta3=1/output^2;