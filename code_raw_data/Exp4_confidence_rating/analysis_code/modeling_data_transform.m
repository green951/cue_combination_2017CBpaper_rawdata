
% transform data for model fitting

 %all 24 subjects
 clear all;
close all;
  IDs=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
 Names={'KunQian','XuanLi','YishuangWang','Oliver','TaoLi2', 'LanDu' , 'WenzheLiu' ,'Ulf','JingGuo','XiaoleiYan','GuangtaoXuan','Lisa','Henning','WeiDing','XiaoyangXu','Fabian','Antonia','Katharina','Robert','Katharina','Jasmin','Sandra','Ralf','Tessa'};

 load('modeling_data.mat')
 
 %for i=1:1
 for i=1:length(IDs)
     cen_left=output_cen_left{i};
     cen_right=output_cen_right{i};
     
     % landmark defined target location in the conflict condition
     targ_lm_conf(1,:)=targ_lm_defined{1,4}(i,:); % left side
     targ_lm_conf(2,:)=targ_lm_defined{2,4}(i,:); % right side
     
    % responses=output_responses;
     
     % calculate distribution distances
     % the left side
     
     cen=cen_left;
     % combination condition
     cond=3;
     d_land(1,1)  =sqrt(      (cen(cond,1)-cen(1,1))^2 +  (cen(cond,2)-cen(1,2))^2           );
     d_motion(1,1)=sqrt(      (cen(cond,1)-cen(2,1))^2 +  (cen(cond,2)-cen(2,2))^2           );
     
     % conflict condition
     cond=4;
     d_land(1,2)  =sqrt(      (targ_lm_conf(1,1)+cen(1,1)-cen(cond,1))^2 +  (targ_lm_conf(1,2)+cen(1,2)-cen(cond,2))^2           );
     d_motion(1,2)=sqrt(      (cen(cond,1)-cen(2,1))^2 +  (cen(cond,2)-cen(2,2))^2           );
     
     % right side
     % combination condition
     cen=cen_right;
     cond=3;
     d_land(2,1)  =sqrt(      (cen(cond,1)-cen(1,1))^2 +  (cen(cond,2)-cen(1,2))^2           );
     d_motion(2,1)=sqrt(      (cen(cond,1)-cen(2,1))^2 +  (cen(cond,2)-cen(2,2))^2           );
     
     % conflict condition
     cond=4;
     d_land(2,2)  =sqrt(      (targ_lm_conf(2,1)+cen(1,1)-cen(cond,1))^2 +  (targ_lm_conf(2,2)+cen(1,2)-cen(cond,2))^2           );
     d_motion(2,2)=sqrt(      (cen(cond,1)-cen(2,1))^2 +  (cen(cond,2)-cen(2,2))^2           );
     
     % average the two sides
     
     % combination condition
     d_land_ave(1)  =(d_land(1,1)+ d_land(2,1))/2;
     d_motion_ave(1)=(d_motion(1,1)+ d_motion(2,1))/2;
     % conflict condition
     d_land_ave(2)  =(d_land(1,2)+ d_land(2,2))/2;
     d_motion_ave(2)=(d_motion(1,2)+ d_motion(2,2))/2;
     
     % scale the distances in combination condition
     % since distances are larger in conflict than combination condition
     % but both observed weights are equally important
     sf=(d_land_ave(2)+d_motion_ave(2))/((d_land_ave(1)+d_motion_ave(1)));
     d_land_ave(1) =d_land_ave(1) *sf;     
     d_motion_ave(1)=d_motion_ave(1)*sf;
     
     % now we can average the distances of combination and conflict condition
     d_land_aveave=(d_land_ave(1)+d_land_ave(2))/2;
     d_motion_aveave=(d_motion_ave(1)+d_motion_ave(2))/2;
  
     
     % transform responses
     for k=1:4
         responses=output_responses{i,k};
         responses_new{i,k}= sqrt( responses(:,1).^2+responses(:,2).^2) .* sign(responses(:,1));
         responses_new{i,k}=responses_new{i,k}-mean(responses_new{i,k}); % mean centered
         
         if k==1
             responses_new{i,k}=responses_new{i,k}-d_land_aveave; % shift to the left by the distance 
         elseif k==2
             responses_new{i,k}=responses_new{i,k}+d_motion_aveave; % shift to the right by the distance 
         end
         
     end
    
     
     % calculate bayesian statistics using the transformed responses
     
     resp_bay=var(responses_new{i,1}) * var(responses_new{i,2})/(var(responses_new{i,1}) + var(responses_new{i,2}));
     resp_bay=sqrt(resp_bay);
     
     d1=abs(mean(responses_new{i,1}));
     d2=abs(mean(responses_new{i,2}));
     weight_obs=d2/(d1+d2);
     
     weight_bay=var(responses_new{i,2})/(var(responses_new{i,1}) + var(responses_new{i,2}));
     
     bay_variability(i,:)=[std(responses_new{i,1}),std(responses_new{i,2}),std(responses_new{i,3}),std(responses_new{i,4}),resp_bay];
     bay_weight(i,:)=[weight_bay,weight_obs];
     
          
 end    
     
filename='modeling_data_transformed_exp4.mat';
save(filename,'output_cen_left','output_cen_right','output_responses','responses_new');
     
     