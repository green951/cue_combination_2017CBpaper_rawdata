clear all;


% closer, unstable
IDs1=[4 5 6 7 8 9 10 11 12 13 26 27 30 29 28 31 32 33 34 35];
Names1={'Rachel' 'Emma' 'Katherine' 'Savannah' 'Taylor' 'Ziwei' 'Mohamed' 'Chris' 'Joseph' 'Jacob' 'Marjorie' 'Cal' 'Ellen' 'Vibhuti' 'Amy' 'Ekom' 'Jonathan' 'Alexander' 'Qiaowei' 'Sparsh'};
 % closer, stable
IDs2=[1 2 3 4 5 6 7 8 9 10 11 23 28 25 24 26 27 29 30 31];
Names2={'Kenzie' 'Madison' 'Marg22' 'Alexis' 'Holly' 'Steven' 'Ryan' 'Joshua' 'Paul' 'Bradley' 'Thomas' 'Eunice' 'Dana' 'Angie' 'Michael' 'Ricks' 'Isaiah' 'Zachary' 'Evan' 'David'};



group=input('group number, 1-unstable,2-stable >>>>>');
dimen=input('enter the analysis dimension (1-2d; 2-x ; 3-z)');

quar=3;

if group==1
    IDs=IDs1;
    Names=Names1;
elseif group==2
    IDs=IDs2;
    Names=Names2;
end
    


%
%
%%%%%%%%%%%%%%%%%%%%%%%%%
index=0;
index1=0;
subjects=[];
idx1=0;
idx2=0;
idx3=0;
idx4=0;

if group==1
    data_by_landpos={};    % group responses by landpos across subjects
    for l=1:10
        data_by_landpos{l}=[];
    end
end
        
for i=1:length(IDs)
    id=IDs(i);
    name=Names{i};
    if distance==1
     
        cd('/Volumes/my_book/toshiba_ext_2_files/experiments/CueCombination/unstable_landmarks/closer_landmarks');

    if group==1
        filename=['VR_CueCombinationUnstable_Response' num2str(id) '_' name '.dat'];
    elseif group==2
        filename=['VR_CueCombinationStable_Response' num2str(id) '_' name '.dat'];
    end
    

    data0=textread(filename);
    
   % cd  E:\Xiaoli--2014-11-16\Green\Documents\Documents_thinkpad_05_28_2013\Documents\experiments\CueCombination\unstable_landmarks
    %cd('/Users/chenx/Documents/experiments/CueCombination/unstable_landmarks');
    cd('/Volumes/my_book/toshiba_ext_2_files/experiments/CueCombination/unstable_landmarks');
    gender(i,1)=data0(1,2);
    condition=data0(:,5);
    x_origin=data0(:,10);
    z_origin=data0(:,11);
  
    if group==2  % stable group
        landpos(i)=data0(1,6);
        sym(i)=syms(landpos(i)+1);
    elseif group==1  % unstable group
        landpos=data0(:,6);
       
    end
    
   [r,c]=size(data0);
    if group==1
            if side==1
                index=find(x_origin>0); % side, left of right
            elseif side==2
                index=find(x_origin<0);
            elseif side==3
                index=find(x_origin >0 | x_origin <0);
            end
        data=data0(index,:);
    elseif group==2
        if sym(i)==1
            if side==1
                index=find(x_origin>0);
                data=data0(index,:);
            elseif side==2
                index=find(x_origin<0);
                data=data0(index,:);
            elseif side==3
                data=data0;
            end 
        elseif sym(i)==-1
            if side==1
                index=find(x_origin>0);
                data=data0(index,:);
            elseif side==2
               index=find(x_origin<0);
                data=data0(index,:);
            elseif side==3
               data=data0;
            end
           
        elseif sym(i)==0 || sym(i)==3
            data=data0;
        end
    end        
    %gender(i)=data(1,2);
    condition=data(:,5);
    x_origin=data(:,10);
    z_origin=data(:,11);
    jitterx=data(:,12);
    jitterz=data(:,13);
    x_resp=data(:,16);
    z_resp=data(:,17);    
    rt=data(:,18);  

    
    [respx,respz]=CoordinateTransform(x_resp,z_resp,x_origin+jitterx,z_origin+jitterz);
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculating variability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for k=1:4
     con_index_left{k}=find(condition==k & x_origin<0);
     con_index_right{k}=find(condition==k & x_origin>0);
     
     %consider both sides when identify outliers
     [response0_left{k},response0_right{k},idx_left{k},idx_right{k}]=BoxPlotOutlier2sides([respx(con_index_left{k}),respz(con_index_left{k})],[respx(con_index_right{k}),respz(con_index_right{k})],quar);
     
      if dimen==1
                response_left{k}=response0_left{k};
                response_right{k}=response0_right{k};
            elseif dimen==2
                response_left{k}=response0_left{k}(:,1);
                response_right{k}=response0_right{k}(:,1);
            elseif dimen==3
                response_left{k}=response0_left{k}(:,2);
                response_right{k}=response0_right{k}(:,2);
      end
            
      cen_left{k,1}(i,:)=mean(response_left{k});
      cen_right{k,1}(i,:)=mean(response_right{k});
            
     centers_left{i,k}=ones(length(con_index_left{k}),1)*mean(response_left{k});
     centers_right{i,k}=ones(length(con_index_right{k}),1)*mean(response_right{k});
     
     index_left{k}=con_index_left{k}(idx_left{k});
     index_right{k}=con_index_right{k}(idx_right{k});
     index_2sides{k}=[index_left{k};index_right{k}];
     
     if dimen==1
      % covariance between 2 dimensions
        [r,p]=corr(response_left{k}(:,1),response_left{k}(:,2),'type','Spearman');
        covar_left(i,k)=r;
        [r,p]=corr(response_right{k}(:,1),response_right{k}(:,2),'type','Spearman');
        covar_right(i,k)=r;
        
%          % pooled across subjects with centroids subtracted
%         m1=mean(response_left{k});
%         resp_left_all{i,k}=[response_left{k}(:,1)-m1(1),response_left{k}(:,2)-m1(2)];
%         m2=mean(response_right{k});
%         resp_right_all{i,k}=[response_right{k}(:,1)-m2(1),response_right{k}(:,2)-m2(2)];
        
        % pooled across subjects with centroids NOT subtracted
        resp_left_all{i,k}=[response_left{k}(:,1),response_left{k}(:,2)];
        resp_right_all{i,k}=[response_right{k}(:,1),response_right{k}(:,2)];
        
     end
        
     %%%% pool responses based on landpos across subjects, only for vision
     %%%% condition
     if group==1
     if k==1
         land_item=0:1:9;
         for l=1:length(land_item)
             temp=intersect(find(landpos==land_item(l)),index_2sides{1});  % find common items
             if isempty(temp)==0
                 for s=1:length(temp)
                     data_by_landpos{land_item(l)+1}(end+1,:)=[respx(temp(s)),respz(temp(s))];
                 end
             end
         end
     end
     end
         
     
     %%% get the centers for vision-only condition
     centers{i,1}=centers_left{i,1};
     centers{i,2}=centers_right{i,1};

     [r1,c1]=size(response_left{k});
     [r2,c2]=size(response_right{k});
     trial_num(i,k)=r1+r2;
     %trial_del(i,k)=delete_num_left{i,k}+delete_num_right{i,k};
        
    % [respx_left,respz_left]=CoordinateTransform(response_left{k}(:,1),response_left{k}(:,2),centers_left{i,k}(:,1),centers_left{i,k}(:,2));
    % [respx_right,respz_right]=CoordinateTransform(response_right{k}(:,1),response_right{k}(:,2),centers_right{i,k}(:,1),centers_right{i,k}(:,2));
     
      if dimen==1
                [respx_left,respz_left]=CoordinateTransform(response_left{k}(:,1),response_left{k}(:,2),centers_left{i,k}(:,1),centers_left{i,k}(:,2));
                [respx_right,respz_right]=CoordinateTransform(response_right{k}(:,1),response_right{k}(:,2),centers_right{i,k}(:,1),centers_right{i,k}(:,2));
      
%                 % for pooling data across subjects for covariance 2
%                 % dimensions calculation
%       resp_left_all{i,k}=[respx_left,respz_left];
%         resp_right_all{i,k}=[respx_right,respz_right];
        
      elseif dimen==2 || dimen==3
                resp_left= response_left{k}-mean(response_left{k});
                resp_right= response_right{k}-mean(response_right{k});
      end
      
            if dimen==1
                response{k}=[respx_left,respz_left;respx_right,respz_right];
                centers{k}=mean(response{k});
                resp_var(i,k)=variance2d(response{k});
                resp_var_left(i,k)=variance2d(response_left{k});
                resp_var_right(i,k)=variance2d(response_right{k});
            elseif dimen==2 || dimen==3
                response{k}=[resp_left;resp_right];
                centers{k}=mean(response{k});
                resp_var(i,k)=var(response{k});
                resp_var_left(i,k)=var(response_left{k});
                resp_var_right(i,k)=var(response_right{k});
            end      
     
     [a,b]=size(response{k});
     trial_con(i,k)=a;

     resp_std(i,k)=sqrt(resp_var(i,k));
     resp_std_left(i,k)=sqrt(resp_var_left(i,k));
     resp_std_right(i,k)=sqrt(resp_var_right(i,k));
     
     %%% calculate bias
      if dimen==1
                a=mean(response_left{k});
                bias_left(i,k)=sqrt(a(1)^2+a(2)^2);
               a=mean(response_right{k});
                bias_right(i,k)=sqrt(a(1)^2+a(2)^2);
                bias(i,k)=(bias_left(i,k)+bias_right(i,k))/2;
            elseif dimen==2
                a=mean(response_left{k});
                bias_left(i,k)=abs(a(1));
                a=mean(response_right{k});
                bias_right(i,k)=abs(a(1));
                bias(i,k)=(bias_left(i,k)+bias_right(i,k))/2;
            elseif dimen==3
                a=mean(response_left{k});
                bias_left(i,k)=abs(a(1));
                a=mean(response_right{k});
                bias_right(i,k)=abs(a(1));
                bias(i,k)=(bias_left(i,k)+bias_right(i,k))/2;
            end
     
 end
 Baypred_var(i,1)=resp_var(i,1)*resp_var(i,2)/(resp_var(i,1)+resp_var(i,2));
 Baypred_std(i,1)=sqrt(Baypred_var(i));
 opt_index(i,1)=Baypred_std(i)/resp_std(i,3);  % optimality index for combination condition
 opt_index(i,2)=Baypred_std(i)/resp_std(i,4); % optimality index for conflict condition
 p_bay(i,1)=resp_var(i,2)/(resp_var(i,2)+resp_var(i,1));
 p_bay_left(i,1)=resp_var_left(i,2)/(resp_var_left(i,2)+resp_var_left(i,1));
 p_bay_right(i,1)=resp_var_right(i,2)/(resp_var_right(i,2)+resp_var_right(i,1));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate landmark proximity in combination condition
  prox_comb_left(i,1)=RelativeProximityDist(centers_left{i,3}(1,:),centers_left{i,1}(1,:),centers_left{i,2}(1,:));
  prox_comb_right(i,1)=RelativeProximityDist(centers_right{i,3}(1,:),centers_right{i,1}(1,:),centers_right{i,2}(1,:));
  prox_comb(i,1)=(prox_comb_left(i,1)+prox_comb_right(i,1))/2;
  resp_proxinte_comb(i,1)=sqrt( resp_var(i,1)*prox_comb(i,1)^2+resp_var(i,2)*(1-prox_comb(i,1))^2);
  
  resp_proxalt_comb(i,1)=alternation_model_normal_2sides(prox_comb_left(i),centers_left{i,1}(1,:),centers_left{i,2}(1,:),resp_std_left(i,1),...
      resp_std_left(i,2),prox_comb_right(i),centers_right{i,1}(1,:),centers_right{i,2}(1,:),resp_std_right(i,1),...
      resp_std_right(i,2),length(con_index_left{3}),length(con_index_right{3}));
  resp_bayalt_comb(i,1)=alternation_model_normal_2sides(p_bay_left(i),centers_left{i,1}(1,:),centers_left{i,2}(1,:),resp_std_left(i,1),...
      resp_std_left(i,2),p_bay_right(i),centers_right{i,1}(1,:),centers_right{i,2}(1,:),resp_std_right(i,1),...
      resp_std_right(i,2),length(con_index_left{3}),length(con_index_right{3}));
  
  % calculate landmark proximity by shifting the vision center to be aligned
  % with the landmark-defined target location in the conflict condition
  targ_lm_left=landmarkdefinedtarget([x_origin(con_index_left{4})+jitterx(con_index_left{4}),z_origin(con_index_left{4})+jitterz(con_index_left{4})],15,dimen);   %find the landmark defined target locations, which are different on different trials
  targ_lm_right=landmarkdefinedtarget([x_origin(con_index_right{4})+jitterx(con_index_right{4}),z_origin(con_index_right{4})+jitterz(con_index_right{4})],15,dimen);   %find the landmark defined target locations, which are different on different trials
  
    targ_lm_defined_left(i,:)=mean(targ_lm_left);
    targ_lm_defined_right(i,:)=mean(targ_lm_right);

 targlmleft{i}=targ_lm_left;
    targlmright{i}=targ_lm_right;
        
        
  prox_conf_left(i,1)=RelativeProximityDist(centers_left{i,4}(1,:),centers_left{i,1}(1,:)+mean(targ_lm_left),centers_left{i,2}(1,:));
  prox_conf_right(i,1)=RelativeProximityDist(centers_right{i,4}(1,:),centers_right{i,1}(1,:)+mean(targ_lm_right),centers_right{i,2}(1,:));
  prox_conf(i,1)=(prox_conf_left(i,1)+prox_conf_right(i,1))/2;
  resp_proxinte_conf(i,1)=sqrt( resp_var(i,1)*prox_conf(i,1)^2+resp_var(i,2)*(1-prox_conf(i,1))^2);
  
  resp_proxalt_conf(i,1)=alternation_model_normal_2sides(prox_conf_left(i),centers_left{i,1}(1,:)+mean(targ_lm_left),centers_left{i,2}(1,:),resp_std_left(i,1),...
      resp_std_left(i,2),prox_conf_right(i),centers_right{i,1}(1,:)+mean(targ_lm_right),centers_right{i,2}(1,:),resp_std_right(i,1),...
      resp_std_right(i,2),length(con_index_left{4}),length(con_index_right{4}));
  
  resp_bayalt_conf(i,1)=alternation_model_normal_2sides(p_bay_left(i),centers_left{i,1}(1,:)+mean(targ_lm_left),centers_left{i,2}(1,:),resp_std_left(i,1),...
      resp_std_left(i,2),p_bay_right(i),centers_right{i,1}(1,:)+mean(targ_lm_right),centers_right{i,2}(1,:),resp_std_right(i,1),...
      resp_std_right(i,2),length(con_index_left{4}),length(con_index_right{4}));
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % concatenating
    xs=respx(find(condition==1));
     zs=respz(find(condition==1));
   for h=1:length(xs)
       idx1=idx1+1;
       resp_all1(idx1,:)=[xs(h),zs(h)];
   end
   
   xs=respx(find(condition==2));
     zs=respz(find(condition==2));
   for h=1:length(xs)
       idx2=idx2+1;
       resp_all2(idx2,:)=[xs(h),zs(h)];
   end
   
   xs=respx(find(condition==3));
     zs=respz(find(condition==3));
   for h=1:length(xs)
       idx3=idx3+1;
       resp_all3(idx3,:)=[xs(h),zs(h)];
   end
   
   xs=respx(find(condition==4));
     zs=respz(find(condition==4));
   for h=1:length(xs)
       idx4=idx4+1;
       resp_all4(idx4,:)=[xs(h),zs(h)];
   end
    
end

 dataouput=[resp_std,Baypred_std,resp_bayalt_comb,resp_bayalt_conf,resp_proxinte_comb,resp_proxinte_conf,resp_proxalt_comb,resp_proxalt_conf,...
     prox_comb,prox_conf,p_bay];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% t test %%%%%%%%%%%%%%%%%%
     
%%%% variability %%%%
std_all=[resp_std,Baypred_std,resp_bayalt_comb,resp_bayalt_conf,resp_proxinte_comb,...
    resp_proxinte_conf,resp_proxalt_comb,resp_proxalt_conf];

%std_all([3,7],:)=[]; % outliers excluded
%variability reduction?
pair{1}=[1 3];% 1 vision vs comb
pair{2}=[2 3]; % 2 motion vs comb
pair{3}=[1 4]; % 3 vision vs conf
pair{4}=[2 4];% 4 motion vs conf
% optimal reduction? bayesian integation
pair{5}=[3 5];% 5 comb vs bay integration prediction
pair{6}=[4 5];% 6 motion vs bay integration prediction
% optimal alternation? bayesian alternation?
pair{7}=[3 6]; % comb
pair{8}=[4 7]; % conf
% proximity_based integration?
pair{9}=[3 8]; % comb
pair{10}=[4 9]; % conf
% proximity_based alternation?
pair{11}=[3 10]; % comb
pair{12}=[4 11]; % conf

index=0;
for i=1:12
    index=index+1;
    [h(index,1) p_main(index,1) ci{index} stats]=ttest(std_all(:,pair{i}(1)),std_all(:,pair{i}(2)));  % 1 vision vs comb
    t(index,1)=stats.tstat;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% proximity vs rr %%%
% prox_comb=data(:,12);
% prox_conf=data(:,13);
% p_bay=data(:,14);
prox_all=[prox_comb,prox_conf,(prox_comb+prox_conf)/2,p_bay];

%prox_all([3,7],:)=[]; % outliers excluded
%prox_all(:,:,2)=[prox_comb(:,2),prox_conf(:,2),(prox_comb(:,2)+prox_conf(:,2))/2,p_bay(:,2)];
pair{1}=[1 4]; % comb vs rr
pair{2}=[2 4]; % conf vs rr
pair{3}=[3 4]; % comb&conf vs rr

 index=0;
    for i=1:3
        index=index+1;
        [h(index,1) p_prox(index,1) ci{index} stats]=ttest(prox_all(:,pair{i}(1)),prox_all(:,pair{i}(2)));  % 7 comb vs bay integration prediction
        t_prox(index,1)=stats.tstat;
    end

%correlation between proximity and rr
pair{1}=[1 4];
pair{2}=[2 4];
pair{3}=[3 4];
index=0;

    index=0;
    for i=1:3
        index=index+1;
        [r0,p0]=corrcoef(prox_all(:,pair{i}(1)),prox_all(:,pair{i}(2)));
        r(index,1)=r0(1,2);
        p_r(index,1)=p0(1,2);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %% calculate covariance betwen dimensions
% data pooled across subjects after center-subtraction
 %resp_left_all{i,k}=[respx_left,respz_left];
 %resp_right_all{i,k}=[respx_right,respz_right];
 
 [is,ks]=size(resp_left_all);
 
 respLeftAll=[];
 for k=1:ks
     index=1;
     for i=1:is         
         temp=resp_left_all{i,k}         ;
         [rs,cs]=size(temp);        
         for r=1:rs             
             respLeftAll{k}(index,1)=temp(r,1);
             respLeftAll{k}(index,2)=temp(r,2);
             
             index=index+1;
         end
     end
 end
 
 for k=1:ks
     [r,p]=corr(respLeftAll{k}(:,1),   respLeftAll{k}(:,2)  ,'type','Spearman');
     covar_left_rs(k)=r;
     covar_left_ps(k)=p;
 end


  [is,ks]=size(resp_right_all);
 
 respRightAll=[];
 for k=1:ks
     index=1;
     for i=1:is         
         temp=resp_right_all{i,k}   ;      
         [rs,cs]=size(temp);        
         for r=1:rs             
             respRightAll{k}(index,1)=temp(r,1);
             respRightAll{k}(index,2)=temp(r,2);
             
             index=index+1;
         end
     end
 end

 for k=1:ks
     [r,p]=corr(respRightAll{k}(:,1),   respRightAll{k}(:,2)  ,'type','Spearman');
     covar_right_rs(k)=r;
     covar_right_ps(k)=p;
 end
 
  % pool targ_lm_left
 targlmLeftAll=[];
 index=1;
 for i=1:is
     temp=targlmleft{i}         ;
     [rs,cs]=size(temp);
     for r=1:rs
         targlmLeftAll(index,1)=temp(r,1);
         targlmLeftAll(index,2)=temp(r,2);
         index=index+1;
     end
 end
 
   % pool targ_lm_right
 targlmRightAll=[];
 index=1;
 for i=1:is
     temp=targlmright{i}         ;
     [rs,cs]=size(temp);
     for r=1:rs
         targlmRightAll(index,1)=temp(r,1);
         targlmRightAll(index,2)=temp(r,2);
         index=index+1;
     end
 end
 



