function [output1,output2,output3,output4,trial_num]=BoxPlotOutlier4quad(response1, response2,response3, response4,quar)

  %      [response{1,k},response{2,k},response{3,k},response{4,k}]=BoxPlotOutlier4quad(response0{1,k},response0{2,k},response0{3,k},response0{4,k},quar);
 responses{1}=response1;
 responses{2}=response2;
 responses{3}=response3;
  responses{4}=response4;

  for q=1:length(responses)
      centers=mean(responses{q});
      [r(q) c]=size(responses{q});
      
      for i=1:r(q)
          dist(q,i)=sqrt((responses{q}(i,1)-centers(1))^2+(responses{q}(i,2)-centers(2))^2);
      end
  end


%find the boxplot boundary
dist=[dist(1,:),dist(2,:),dist(3,:),dist(4,:)];
dist_sort=sort(dist);

quartile3=dist_sort(round(length(dist)*.75));
quartile1=dist_sort(round(length(dist)*.25));
band=(quartile3-quartile1)*quar;
upper=quartile3+band;
lower=quartile1-band;


for q=1:length(responses)
    index=1;
    for i=1:r(q)
        if dist(q,i)<upper
            output{q}(index,:)=responses{q}(i,:);
            index=index+1;
        end
    end
    
    trial_num(q)=length(output{q}(:,1));  % remaining trial in each quadrant
end
output1=output{1};
output2=output{2};
output3=output{3};
output4=output{4};

end
