function output=combination_model(p1,mu1,mu2,leta1,leta2)

Ntrials=500;
data=zeros(Ntrials,3);
for i=1:Ntrials
    distribParams1(1) = mu1;
    distribParams1(2) = leta1;
    data(i,1) = randraw('vonmises', distribParams1, 1);
    
    distribParams2(1) = mu2;
    distribParams2(2) = leta2;
    data(i,2) = randraw('vonmises', distribParams2, 1);
    
    data(i,3)=p1*data(i,1)+(1-p1)*data(i,2);
    
end

output=circ_std(data(:,3));
%leta3=1/output^2;