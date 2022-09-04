function output=variance2d(xzs)
xs=xzs(:,1);
zs=xzs(:,2);
center=[mean(xs),mean(zs)];
temp=zeros(1,length(xs));
for i=1:length(xs),
    temp(i)=(xs(i)-center(1))^2+(zs(i)-center(2))^2;
end
%output=sqrt(sum(temp)/(length(temp)-1));
output=sqrt(sum(temp)/(length(temp)-1));
output=output.^2;
end