function output=MeanSquare(xs,zs,target)
center=target;
for i=1:length(xs),
    dist(i)=(xs(i)-center(1))^2+(zs(i)-center(2))^2;
end
output=sqrt(sum(dist)/length(dist));
end