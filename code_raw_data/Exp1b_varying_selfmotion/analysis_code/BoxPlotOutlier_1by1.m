function [output,outliers]=BoxPlotOutlier_1by1(responses, quar)

[r c]=size(responses);

if r>3
    index=1;
    for i=1:r
        responses_temp=responses;
        responses_temp(i,:)=[]; % take the current response point out of consideration when calculating outlier bounds
        centers=mean(responses_temp);
        dist=sqrt((responses(:,1)-centers(1)).^2+(responses(:,2)-centers(2)).^2);
        
        dist_sort=sort(dist);
        quartile3=dist_sort(round(length(dist)*.75));
        quartile1=dist_sort(round(length(dist)*.25));
        band=(quartile3-quartile1)*quar;
        
        
        if band ~= 0
            upper=quartile3+band;
            lower=quartile1-band;
        else
            upper=max(dist);
            lower=min(dist);
        end
        
        %if (dist(i)>lower || dist(i)==lower)  && (dist(i)<upper || dist(i)==upper)
        if  (dist(i)<upper || dist(i)==upper)
            output(index,:)=responses(i,:);
            index=index+1;
        end
    end
    
    [r1,c1]=size(output);
    outliers=r-r1;
    
else
    output=responses;
    outliers=0;
end
