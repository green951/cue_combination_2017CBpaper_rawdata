function [output,outliers]=BoxPlotOutlier(responses, quar)

centers=mean(responses);
[r c]=size(responses);
for i=1:r
    %dist(i)=sqrt(responses(i,1)^2+responses(i,2)^2);
     dist(i)=sqrt((responses(i,1)-centers(1))^2+(responses(i,2)-centers(2))^2);
end

[dist_sort,index]=sort(dist);
quartile3=dist_sort(round(length(dist)*.75));
quartile1=dist_sort(round(length(dist)*.25));
band=(quartile3-quartile1)*quar;
upper=quartile3+band;
lower=quartile1-band;

index=1;
for i=1:r
    %if dist(i)>lower && dist(i)<upper
    % if (dist(i)>lower || dist(i)==lower)  && (dist(i)<upper || dist(i)==upper)
    if (dist(i)<upper || dist(i)==upper)
        output(index,:)=responses(i,:);
        index=index+1;
    end
end

[r1,c1]=size(output);
outliers=r-r1;
end
