function output=OptimalWeightVonmises(mu1,leta1,mu2,leta2)

% dist=sqrt((center1(1)-center2(1))^2+(center1(2)-center2(2))^2);
% mu1=0;
% mu2=dist;
distribParams1(1) = mu1;
distribParams1(2) =leta1;
distribParams2(1) = mu2;
distribParams2(2) =leta2;

Ntrials=500;
w1=0:0.02:1;
w2=1-w1;
for i=1:length(w1)
    w=w1(i);
    data=zeros(Ntrials,3);
    for n=1:Ntrials
    data(n,1) = randraw('vonmises', distribParams1, 1);
    data(n,2) = randraw('vonmises', distribParams2, 1);
    data(n,3)=data(n,1)*w+data(n,2)*(1-w);
    end
%     vari(i,1)=circ_var(data(:,1));
%     vari(i,2)=circ_var(data(:,2));
     vari(i,3)=circ_var(data(:,3));
end

% find the lowest vari(3)
index=find(vari(:,3)==min(vari(:,3)));
opt_w=w1(index);
output=opt_w;
end