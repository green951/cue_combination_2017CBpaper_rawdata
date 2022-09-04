clear all;
close all;

currentFolder = pwd;


% Names={'Rebecca'};
 dimen=input('enter the analysis dimension (1-2d; 2-x ; 3-z)');
dayChose=input('enter the day to analyze (1, 2, 3)');


 %all 24 subjects
  IDs=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
 Names={'KunQian','XuanLi','YishuangWang','Oliver','TaoLi2', 'LanDu' , 'WenzheLiu' ,'Ulf','JingGuo','XiaoleiYan','GuangtaoXuan','Lisa','Henning','WeiDing','XiaoyangXu','Fabian','Antonia','Katharina','Robert','Katharina','Jasmin','Sandra','Ralf','Tessa'};

%  % 22 subjects with scans
%   IDs=[1,2,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
%   Names={'KunQian','XuanLi','Oliver','TaoLi2', 'LanDu'  ,'Ulf','JingGuo','XiaoleiYan','GuangtaoXuan','Lisa','Henning','WeiDing','XiaoyangXu','Fabian','Antonia','Katharina','Robert','Katharina','Jasmin','Sandra','Ralf','Tessa'};


 
  
%  IDs=[20,21]; 
%  Names={'Katharina','Jasmin'};
% %  
%    IDs=[19,20,21,22];
%   Names={'Robert','Katharina','Jasmin','Sandra'};

  


quar=30;
%%%%%%%%%%%%%%%%%%%%%%%%%
index=0;
index1=0;
subjects=[];
idx1=0;
idx2=0;
idx3=0;
idx4=0;
for i=1:length(IDs),
    id=IDs(i);
    name=Names{i};
  
      cd([currentFolder '/data/sub_' num2str(id)]) 
  
      
        filename=['CueCombination_Confidence_Response_' num2str(id) '_' name '.dat'];
  
    data0=textread(filename);
    day=data0(:,20);
    if dayChose==1
        data=data0(find(day==1),:);
    elseif dayChose==2
        data=data0(find(day==2),:);
    elseif dayChose==3
        data=data0;
    end
    
      cd(currentFolder) 
    gender(i,1)=data(1,2);
  

 
    trialNum=data(:,3);
    elapsedTime=data(:,4);
    condition=data(:,5);
 
    post1=data(:,6);
    post2=data(:,7);
    post3=data(:,8);
    
    targx=data(:,9);
    targz=data(:,10);
    jitterx=data(:,11);
    jitterz=data(:,12);
    
    startx=data(:,13);
    startz=data(:,14);
    respx_org=data(:,15);% before coordinate transformation
    respz_org=data(:,16); % before coordinate transformation
    rt=data(:,17);
     conf=data(:,18);
     confRt=data(:,19);
    % day=data(:,20);
     
%      % convert conf to rank score
      [a,confRank]=sort(conf);
      
      %%% range of correct distances
      cor_dists=sqrt((targx+jitterx-startx).^2+(targz+jitterz-startz).^2);
      
      cor_dist={};
      for n=1:4
      index=find(condition==n);
      cor_dist{n}=cor_dists(index);
      dist_min(i,n)=min(cor_dist{n});
      dist_max(i,n)=max(cor_dist{n});
      dist_std(i,n)=std(cor_dist{n});
      end
      
  
     landx=0;
     landz=4.5;
     
     
     respx=respx_org;
     respz=respz_org;
     
    %[respx,respz]=CoordinateTransform(respx_org,respz_org,targx+jitterx,targz+jitterz);
    
    % define the coordinate relative to the landmark
   % originx=landx;
   % originz=landz;
    %[respx,respz]=CoordinateTransform_new(respx_org,respz_org,targx+jitterx,targz+jitterz,originx,originz);
    
    dist_to_landmark=sqrt((respx-landx).^2+(respz-landz).^2) ;
    dist_to_landmark_error= abs(dist_to_landmark - sqrt( (landx-targx-jitterx).^2+(landz-targz-jitterz).^2 ));
    
    for k=1:4
        index=find(condition==k);
        distlandmarkerror(i,k)=mean(dist_to_landmark_error);
    end
    
    
%     disterr_all=sqrt( respx.^2+respz.^2  );
%     for n=1:2
%         for k=1:4
%             index=find(condition==k & environment==n);
            
    
    
    %for n=1:2
    condN=4;
        for k=1:condN
            % % %             % % environments analyzed separately
            con_index_left{k}=find(condition==k  & targx<0);
            con_index_right{k}=find(condition==k  & targx>0);
            resp_time(i,k)=mean(rt(   [con_index_left{k};con_index_right{k}] ));    
            %
           
            %              %consider both sides when identify outliers
            
%             [response0_left{k},response0_right{k},conf_left{k},conf_right{k},confRank_left{k},confRank_right{k},delete_num_left{i,k},delete_num_right{i,k},left_ind,right_ind]=BoxPlotOutlier2sides_distErr([respx(con_index_left{k}),respz(con_index_left{k})],...
%                 [respx(con_index_right{k}),respz(con_index_right{k})],conf(con_index_left{k}),conf(con_index_right{k}),confRank(con_index_left{k}),confRank(con_index_right{k}),quar,con_index_left{k},con_index_right{k});
%           
            [response0_left{k},response0_right{k},conf_left{k},conf_right{k},confRank_left{k},confRank_right{k},delete_num_left{i,k},delete_num_right{i,k},left_ind,right_ind]=BoxPlotOutlier2sides([respx(con_index_left{k}),respz(con_index_left{k})],...
                [respx(con_index_right{k}),respz(con_index_right{k})],conf(con_index_left{k}),conf(con_index_right{k}),confRank(con_index_left{k}),confRank(con_index_right{k}),quar,con_index_left{k},con_index_right{k});
          
            
            conf2sides{k}=[conf_left{k};conf_right{k}];
            confRank2sides{k}=[confRank_left{k};confRank_right{k}];
            
            RespToTarg{k}=[response0_left{k};response0_right{k}];
            ErrToTarg{k}=[];
            
            for j=1:length(RespToTarg{k})
                x=RespToTarg{k}(j,1);
                z=RespToTarg{k}(j,2);
                if dimen==1
                    ErrToTarg{k}(end+1,1)=sqrt(x^2+z^2);
                elseif dimen==2
                    ErrToTarg{k}(end+1,1)=abs(x);
                elseif dimen==3
                    ErrToTarg{k}(end+1,1)=abs(z);
                end
            end
                
            
            % calculate accuracy, distance to target location
            temp1=response0_left{k};
            temp2=response0_right{k};
            [r1,c1]=size(temp1);
            [r2,c2]=size(temp2);
           disterr{1}(i,k)=( sum(sqrt (temp1(:,1).^2+temp1(:,2).^2))+sum(sqrt (temp2(:,1).^2+temp2(:,2).^2)))/(r1+r2);  % 2d distance error
           disterr{2}(i,k)=( sum(abs(temp1(:,1)))+sum(abs(temp2(:,1))))/(r1+r2);  % x dimension distance error
           disterr{3}(i,k)=( sum(abs(temp1(:,2)))+sum(abs(temp2(:,2))))/(r1+r2);  % z dimension distance error
            
            
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
                            
             if dimen==1
            % covariance between 2 dimensions
            [r,p]=corr(response_left{k}(:,1),response_left{k}(:,2),'type','Spearman');
            covar_left(i,k)=r;
            [r,p]=corr(response_right{k}(:,1),response_right{k}(:,2),'type','Spearman');
            covar_right(i,k)=r;
            
%               % pooled across subjects with centroids subtracted
%         m1=mean(response_left{k});
%         resp_left_all{i,k}=[response_left{k}(:,1)-m1(1),response_left{k}(:,2)-m1(2)];
%         m2=mean(response_right{k});
%         resp_right_all{i,k}=[response_right{k}(:,1)-m2(1),response_right{k}(:,2)-m2(2)];
        
         % pooled across subjects with centroids NOT subtracted
        resp_left_all{i,k}=[response_left{k}(:,1),response_left{k}(:,2)];
        resp_right_all{i,k}=[response_right{k}(:,1),response_right{k}(:,2)];
             end
             
        
            
            [r1,c1]=size(response_left{k});
            [r2,c2]=size(response_right{k});
            
            trial_num(i,k)=r1+r2;
            trial_del(i,k)=delete_num_left{i,k}+delete_num_right{i,k};
            
            output_cen_left{i,1}(k,:)=mean(response_left{k});
            output_cen_right{i,1}(k,:)=mean(response_right{k});
            
            
            %%%%% calculate bias
            if dimen==1
                a=mean(response_left{k});
                bias_left{k}(i,1)=sqrt(a(1)^2+a(2)^2);
               a=mean(response_right{k});
                bias_right{k}(i,1)=sqrt(a(1)^2+a(2)^2);
                bias{k}(i,1)=(bias_left{k}(i,1)+bias_right{k}(i,1))/2;
            elseif dimen==2
                a=mean(response_left{k});
                bias_left{k}(i,1)=abs(a(1));
                a=mean(response_right{k});
                bias_right{k}(i,1)=abs(a(1));
                bias{k}(i,1)=(bias_left{k}(i,1)+bias_right{k}(i,1))/2;
            elseif dimen==3
                a=mean(response_left{k});
                bias_left{k}(i,1)=abs(a(1));
                a=mean(response_right{k});
                bias_right{k}(i,1)=abs(a(1));
                bias{k}(i,1)=(bias_left{k}(i,1)+bias_right{k}(i,1))/2;
            end
            
            centers_left{i,k}=ones(length(con_index_left{k}),1)*mean(response_left{k});
            centers_right{i,k}=ones(length(con_index_right{k}),1)*mean(response_right{k});
            
            
            
            if dimen==1  %correct intrinsic bias in each side before the pooling
                [respx_left,respz_left]=CoordinateTransform(response_left{k}(:,1),response_left{k}(:,2),centers_left{i,k}(:,1),centers_left{i,k}(:,2));
                [respx_right,respz_right]=CoordinateTransform(response_right{k}(:,1),response_right{k}(:,2),centers_right{i,k}(:,1),centers_right{i,k}(:,2));
          
%                 % for the pooling across subjects, covariance 2 dimensions
%                 resp_left_all{i,k}=[respx_left,respz_left];
%                 resp_right_all{i,k}=[respx_right,respz_right];
            
            elseif dimen==2 || dimen==3
                resp_left= response_left{k}-mean(response_left{k});
                resp_right= response_right{k}-mean(response_right{k});
            end
            
            if dimen==1
                response{k}=[respx_left,respz_left;respx_right,respz_right];
                output_responses{i,k}=response{k};
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
                             
            resp_std(i,k)=sqrt(resp_var(i,k));
            resp_std_left(i,k)=sqrt(resp_var_left(i,k));
            resp_std_right(i,k)=sqrt(resp_var_right(i,k));
            
            RespToCenter{k}=response{k};
            ErrToCenter{k}=[];
            for j=1:length(RespToCenter{k})
                x=RespToCenter{k}(j,1);
                
                if dimen==1
                    z=RespToCenter{k}(j,2);
                    ErrToCenter{k}(end+1,1)=sqrt(x^2+z^2);
                elseif dimen==2 || dimen==3
                    ErrToCenter{k}(end+1,1)=abs(x);
                end
            end
            
            
            %%% calculate retrospective ability for each condition
            % pearson R correlation
            temp=corrcoef(ErrToCenter{k},conf2sides{k});  % error=distance to center
            retroAbil{1}(i,k)=-temp(1,2);
            temp=corrcoef(ErrToTarg{k},conf2sides{k});   % error=distance to target
            retroAbil{2}(i,k)=-temp(1,2);
%             
%             % rank score correlation
%             [a,errRank]=sort(ErrToCenter{k});
%             temp=corrcoef(errRank,confRank2sides{k});  % error=distance to center
%             retroAbil{1}(i,k)=-temp(1,2);
%             
%              [a,errRank]=sort(ErrToTarg{k});
%             temp=corrcoef(errRank,confRank2sides{k});  % error=distance to center
%             retroAbil{2}(i,k)=-temp(1,2);
%                         
             confAve(i,k)=mean(conf2sides{k});
            
            
          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        %%% calculate retrospective ability across conditions, pooled
        
        confAll=[];
        confRankAll=[];
        errAll1=[];
        errAll2=[];
        for k=1:condN
            for n=1:length(conf2sides{k})
                confAll(end+1,1)=conf2sides{k}(n);
                confRankAll(end+1,1)=confRank2sides{k}(n);
                errAll1(end+1,1)=ErrToCenter{k}(n);
                errAll2(end+1,1)=ErrToTarg{k}(n);
                
            end
        end
        
     %%%% correlating raw scores   
        %%% error= distance to center
        temp=corrcoef(errAll1,confAll);
        retroAbil{1}(i,condN+1)=-temp(1,2);
        
        %%% error = distance to target
         temp=corrcoef(errAll2,confAll);
        retroAbil{2}(i,condN+1)=-temp(1,2); 

%    %  %% correlating rank scores   
%        % %% error= distance to center
%         [a,errRank]=sort(errAll1);
%         temp=corrcoef(errRank,confRankAll);
%         retroAbil{1}(i,condN+1)=-temp(1,2);
%         
%       %  %% error = distance to target
%          [a,errRank]=sort(errAll2);
%          temp=corrcoef(errRank,confRankAll);
%         retroAbil{2}(i,condN+1)=-temp(1,2);  
        
        
        
        
        
        Baypred_var(i)=resp_var(i,1)*resp_var(i,2)/(resp_var(i,1)+resp_var(i,2));
        Baypred_std(i)=sqrt(Baypred_var(i));
        p_bay(i)=resp_var(i,2)/(resp_var(i,2)+resp_var(i,1));
        p_bay_left(i)=resp_var_left(i,2)/(resp_var_left(i,2)+resp_var_left(i,1));
        p_bay_right(i)=resp_var_right(i,2)/(resp_var_right(i,2)+resp_var_right(i,1));
        
      
        %%%%% calculate landmark proximity in combination condition
        prox_comb_left(i)=RelativeProximityDist(centers_left{i,3}(1,:),centers_left{i,1}(1,:),centers_left{i,2}(1,:));
        prox_comb_right(i)=RelativeProximityDist(centers_right{i,3}(1,:),centers_right{i,1}(1,:),centers_right{i,2}(1,:));
        prox_comb(i)=(prox_comb_left(i)+prox_comb_right(i))/2;
        resp_proxinte_comb(i)=sqrt( resp_var(i,1)*prox_comb(i)^2+resp_var(i,2)*(1-prox_comb(i))^2);
        
        resp_proxalt_comb(i)=alternation_model_normal_2sides(prox_comb_left(i),centers_left{i,1}(1,:),centers_left{i,2}(1,:),resp_std_left(i,1),...
            resp_std_left(i,2),prox_comb_right(i),centers_right{i,1}(1,:),centers_right{i,2}(1,:),resp_std_right(i,1),...
            resp_std_right(i,2),length(con_index_left{3}),length(con_index_right{3}));
        resp_bayalt_comb(i)=alternation_model_normal_2sides(p_bay_left(i),centers_left{i,1}(1,:),centers_left{i,2}(1,:),resp_std_left(i,1),...
            resp_std_left(i,2),p_bay_right(i),centers_right{i,1}(1,:),centers_right{i,2}(1,:),resp_std_right(i,1),...
            resp_std_right(i,2),length(con_index_left{3}),length(con_index_right{3}));
        
        
        if condN==4
        % calculate landmark proximity by shifting the vision center to be aligned
        % with the landmark-defined target location in the conflict condition
        
        targ_lm_left=landmarkdefinedtarget([targx(con_index_left{4})+jitterx(con_index_left{4}),targz(con_index_left{4})+jitterz(con_index_left{4})],15,dimen);   %find the landmark defined target locations, which are different on different trials
        targ_lm_right=landmarkdefinedtarget([targx(con_index_right{4})+jitterx(con_index_right{4}),targz(con_index_right{4})+jitterz(con_index_right{4})],15,dimen);   %find the landmark defined target locations, which are different on different trials
        
        
        targ_lm=[targ_lm_left;targ_lm_right];
        
        targ_lm_defined{1,k}(i,:)=mean(targ_lm_left);
        targ_lm_defined{2,k}(i,:)=mean(targ_lm_right);
        
         targlmleft{i}=targ_lm_left;
    targlmright{i}=targ_lm_right;
        
        prox_conf_left(i)=RelativeProximityDist(centers_left{i,4}(1,:),centers_left{i,1}(1,:)+mean(targ_lm_left),centers_left{i,2}(1,:));
        prox_conf_right(i)=RelativeProximityDist(centers_right{i,4}(1,:),centers_right{i,1}(1,:)+mean(targ_lm_right),centers_right{i,2}(1,:));
        prox_conf(i)=(prox_conf_left(i)+prox_conf_right(i))/2;
        resp_proxinte_conf(i)=sqrt( resp_var(i,1)*prox_conf(i)^2+resp_var(i,2)*(1-prox_conf(i))^2);
        
        
        
        resp_proxalt_conf(i)=alternation_model_normal_2sides(prox_conf_left(i),centers_left{i,1}(1,:)+mean(targ_lm_left),centers_left{i,2}(1,:),resp_std_left(i,1),...
            resp_std_left(i,2),prox_conf_right(i),centers_right{i,1}(1,:)+mean(targ_lm_right),centers_right{i,2}(1,:),resp_std_right(i,1),...
            resp_std_right(i,2),length(con_index_left{4}),length(con_index_right{4}));
        resp_bayalt_conf(i)=alternation_model_normal_2sides(p_bay_left(i),centers_left{i,1}(1,:)+mean(targ_lm_left),centers_left{i,2}(1,:),resp_std_left(i,1),...
            resp_std_left(i,2),p_bay_right(i),centers_right{i,1}(1,:)+mean(targ_lm_right),centers_right{i,2}(1,:),resp_std_right(i,1),...
            resp_std_right(i,2),length(con_index_left{4}),length(con_index_right{4}));
        
        end
        
%         VARindex_comb(i,1)= (1/resp_var(i,3))  / (1/resp_var(i,1)+1/resp_var(i,2)); % bayinteg_
%         VARindex_comb(i,2)=(1/resp_var(i,3)) / (1/resp_bayalt_comb(i)^2); %bayalt_
%         VARindex_comb(i,3)= (1/resp_var(i,3)) / (1/resp_proxinte_comb(i)^2); %proxinteg_
%         VARindex_comb(i,4)=(1/resp_var(i,3)) / (1/resp_proxalt_comb(i)^2); %proxalt_
%         
%         
%         VARindex_conf(i,1)= (1/resp_var(i,4))  / (1/resp_var(i,1)+1/resp_var(i,2)); %bayinteg_
%         VARindex_conf(i,2)=(1/resp_var(i,4)) / (1/resp_bayalt_conf(i)^2); %bayalt_
%         VARindex_conf(i,3)= (1/resp_var(i,4)) / (1/resp_proxinte_conf(i)^2); %proxinteg_
%         VARindex_conf(i,4)=(1/resp_var(i,4)) / (1/resp_proxalt_conf(i)^2); %proxalt_
%         
%         
%         
%         RPindex(i,1)=prox_comb(i)/p_bay(i); % _comb
%         RPindex(i,2)=prox_conf(i)/p_bay(i); % _conf
%         RPindex(i,3)=(prox_comb(i)/2+prox_conf(i)/2)/p_bay(i); %_average
    
    %%%%%%%%%%% pool data across subjects %%%%%%%%%%%%%%%
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
    
    if condN==4
        xs=respx(find(condition==4));
        zs=respz(find(condition==4));
        for h=1:length(xs)
            idx4=idx4+1;
            resp_all4(idx4,:)=[xs(h),zs(h)];
        end
    end
    
end
%

if condN==4
    dataouput(:,:,1)=[resp_std,Baypred_std',resp_bayalt_comb',resp_bayalt_conf',resp_proxinte_comb',resp_proxinte_conf',resp_proxalt_comb',resp_proxalt_conf',...
        prox_comb',prox_conf',p_bay'];
else
    dataouput(:,:,1)=[resp_std,Baypred_std',resp_bayalt_comb',resp_proxinte_comb',resp_proxalt_comb',...
        prox_comb',p_bay'];
end
   
std_all(:,:,1)=[resp_std,Baypred_std',resp_bayalt_comb',resp_bayalt_conf',resp_proxinte_comb',...
    resp_proxinte_conf',resp_proxalt_comb',resp_proxalt_conf'];
prox_all=[prox_comb',prox_conf',(prox_comb'+prox_conf')/2,p_bay'];


% save files for computational modeling

filename='modeling_data.mat';
save(filename,'output_cen_left','output_cen_right','output_responses','targ_lm_defined');

%% t test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variability %%%%%%%%%%%%%%%%%%%%%%
std_all(:,:,1)=[resp_std,Baypred_std',resp_bayalt_comb',resp_bayalt_conf',resp_proxinte_comb',...
    resp_proxinte_conf',resp_proxalt_comb',resp_proxalt_conf'];
% std_all(:,:,2)=[resp_std(:,:,2),Baypred_std(:,2),resp_bayalt_comb(:,2),resp_bayalt_conf(:,2),resp_proxinte_comb(:,2),...
%     resp_proxinte_conf(:,2),resp_proxalt_comb(:,2),resp_proxalt_conf(:,2)];
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
        [h(index,1) p_t(index,1) ci{index} stats]=ttest(std_all(:,pair{i}(1)),std_all(:,pair{i}(2)));  % 1 vision vs comb
        t(index,1)=stats.tstat;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% proximity vs rr %%%%%%%%%%%%%%%%%%%%%%%%
prox_all=[prox_comb',prox_conf',(prox_comb'+prox_conf')/2,p_bay'];
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
        r_cor(index,1)=r0(1,2);
        p_r(index,1)=p0(1,2);
    end


% manipulation check
% index=1;
% [h_com p_com(index,1) ci_com stats]=ttest(std_all(:,1,1),std_all(:,1,2));
% t_com(index,1)=stats.tstat;
% index=index+1;
% [h_com p_com(index,1) ci_com stats]=ttest(std_all(:,2,1),std_all(:,2,2));
% t_com(index,1)=stats.tstat;
% index=index+1;
% [h_com p_com(index,1) ci_com stats]=ttest(prox_all(:,4,1),prox_all(:,4,2));
% t_com(index,1)=stats.tstat;
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
%%%% compare variability bayesian index (bayesian integration index & bayesian alternation index) to 1
% for n=1:1
%     VARindex(:,1,n)=VARindex_comb(:,1,n);  % combination, bayesian integration
%     VARindex(:,2,n)=VARindex_comb(:,2,n);  % combination, bayesian alternation
%     VARindex(:,3,n)=VARindex_comb(:,3,n);  % combination, proximity integration
%     VARindex(:,4,n)=VARindex_comb(:,4,n);  % combination, proximity alternation
%         
%     VARindex(:,5,n)=VARindex_conf(:,1,n);  % conflict, bayesian integration
%     VARindex(:,6,n)=VARindex_conf(:,2,n);   % conflict, bayesian alternation
%     VARindex(:,7,n)=VARindex_conf(:,3,n);  % conflict, proximity integration
%     VARindex(:,8,n)=VARindex_conf(:,4,n);   % conflict, proximity alternation
% end

% %bayesian integration,combination
% for i=1:8
%     for n=1:1
%         [mean_ind(i,n),Nsub(i,n),p_ind(i,n),t_ind(i,n)]=IndexTest(VARindex(:,i,n),1);
%     end
%     
% end


%%  data analysis
% detect outlier participants

% for i=1:length(IDs)
%     var_all(i)=sum(resp_var(i,:))/4;
% end
%
% [var_sort,index]=sort(var_all);
% quartile3=var_sort(round(length(var_sort)*.75));
% quartile1=var_sort(round(length(var_sort)*.25));
% range=quartile3-quartile1;
% upper=quartile3+1.5*range;
% lower=quartile1-1.5*range;
% %
% index=0;
% id_index=0;
% for i=1:length(var_all)
%     id_index(i)=i;
%     if var_all(i)<upper+0.005 && var_all(i)>lower-0.005 %&& var_all(i)>.1
%         index=index+1;
%         var_cleaned(index)=var_all(i);
%         %id_index(index)=i;
%     end
% end


% % variances
% RichVarVision=RichVarVision(id_index,:);
% RichVarMotion=RichVarMotion(id_index,:);
% RichVarComb=RichVarComb(id_index,:);
% RichVarConf=RichVarConf(id_index,:);
% PoorVarVision=PoorVarVision(id_index,:);
% PoorVarMotion=PoorVarMotion(id_index,:);
% PoorVarComb=PoorVarComb(id_index,:);
% PoorVarConf=PoorVarConf(id_index,:);
%
% % landmark proximity (2d space)
% rich_prox_stand=rich_prox_stand(id_index);
% poor_prox_stand=poor_prox_stand(id_index);
% rich_prox_comb=rich_prox_comb(id_index);
% poor_prox_comb=poor_prox_comb(id_index);
% rich_prox_conf=rich_prox_conf(id_index);
% poor_prox_conf=poor_prox_conf(id_index);
% rich_prox_stand=rich_prox_stand';
% poor_prox_stand=poor_prox_stand';
%
% RichVisionCenter=RichVisionCenter(id_index,:);
% RichMotionCenter=RichMotionCenter(id_index,:);
% PoorVisionCenter=PoorVisionCenter(id_index,:);
% PoorMotionCenter=PoorMotionCenter(id_index,:);
% RichVisionCenterDist=(RichVisionCenter(:,1).^2+RichVisionCenter(:,2).^2).^(1/2);
% RichMotionCenterDist=(RichMotionCenter(:,1).^2+RichMotionCenter(:,2).^2).^(1/2);
% PoorVisionCenterDist=(PoorVisionCenter(:,1).^2+PoorVisionCenter(:,2).^2).^(1/2);
% PoorMotionCenterDist=(PoorMotionCenter(:,1).^2+PoorMotionCenter(:,2).^2).^(1/2);
% % calculate relative reliability (2d space)
% rich_rr=RichVarMotion(1,:) ./ (RichVarVision(1,:) +RichVarMotion(1,:) );
% poor_rr=PoorVarMotion(1,:) ./ (PoorVarVision(1,:) +PoorVarMotion(1,:) );
%
% % calculate predicted combination var (bayesian model)
% rich_pred_bay=RichVarVision(:,1) .*RichVarMotion(:,1) ./(RichVarVision(:,1) +RichVarMotion(:,1));
% poor_pred_bay=PoorVarVision(:,1) .*PoorVarMotion(:,1) ./(PoorVarVision(:,1) +PoorVarMotion(:,1));
%
% % calculate predicted combination var (alternation model)
% rich_pred_alt= rich_prox_stand .* (.46^2+RichVarVision(:,1)) + (1-rich_prox_stand).*(0^2+RichVarMotion(:,1)) -(rich_prox_stand *.46+0).^2;
% poor_pred_alt= poor_prox_stand .* (.46^2+PoorVarVision(:,1)) + (1-poor_prox_stand).*(0^2+PoorVarMotion(:,1)) -(poor_prox_stand *.46+0).^2;
%%
figure(1)
subplot(2,2,1)
xmin=-4;
xmax=4;
ymin=-6;
ymax=2;
plot(resp_all1(:,1),resp_all1(:,2),'k*');
axis([xmin xmax ymin ymax]);
title('vision')
subplot(2,2,2)
plot(resp_all2(:,1),resp_all2(:,2),'k*');
axis([xmin xmax ymin ymax]);
title('motion')
subplot(2,2,3)
plot(resp_all3(:,1),resp_all3(:,2),'k*');
axis([xmin xmax ymin ymax]);
title('comb')
subplot(2,2,4)
plot(resp_all4(:,1),resp_all4(:,2),'k*');
axis([xmin xmax ymin ymax]);
title('conf')

%% plot response centers for each subject

centers{1,1}=VisionCenter(:,:,1);
centers{1,2}=MotionCenter(:,:,1);
centers{1,3}=CombCenter(:,:,1);
centers{1,4}=ConfCenter(:,:,1);

xmin=-1.5;
xmax=1.5;
ymin=-1.5;
ymax=1.5;
xline=xmin:.05:xmax;
yline=ymin:.05:ymax;
for n=1:1
    for k=1:4
        figure(n);
        subplot(2,2,k)
        points=centers{n,k};
        plot(points(:,1),points(:,2),'b*')
        hold on;
        plot(xline,zeros(length(xline),1),'-')
        hold on;
        plot(zeros(length(yline),1),yline,'-')
        axis([xmin xmax ymin ymax]);
        title(['n=' num2str(n) ', k=' num2str(k)]);
    end
    
end



