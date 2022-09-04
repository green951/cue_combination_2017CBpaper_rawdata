function output=variance1d(xs)
center=mean(xs);
for i=1:length(xs),
    dist(i)=(xs(i)-center)^2;
end
%output=sqrt(sum(dist)/(length(dist)-1));
output=sqrt(sum(dist)/(length(dist)-1));
output=output.^2;
end