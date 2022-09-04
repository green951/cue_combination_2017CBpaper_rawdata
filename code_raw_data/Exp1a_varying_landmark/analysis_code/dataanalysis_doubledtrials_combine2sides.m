clear all;

%cd E:\Xiaoli--2014-6-21\Green\Documents\Documents_thinkpad_05_28_2013\Documents\experiments\CueCombination\2environments_doubledtrials_data_flag
%cd E:\Xiaoli--2014-11-16\Green\Documents\Documents_thinkpad_05_28_2013\Documents\experiments\CueCombination\2environments_doubledtrials_data_flag

IDs=[1 2 3 6 7 9 10,8, 13,16,17,18,19,20,14,21,22,23];
Names={'Anna' 'Sophie' 'Aaron' 'Fernando' 'Courthey'  'Caleb' 'Noa' 'Raven'  'Sam','Alex' 'Michael','Ahmad','Emma','Brett','Robert', 'Bradley', 'Bria', 'Susanna'};
dimen=input('enter the analysis dimension (1-2d; 2-x ; 3-z)');

quar=30;
%%%%%%%%%%%%%%%%%%%%%%%%%
index=0;
index1=0;
subjects=[];
idx1=0;
idx2=0;
idx3=0;
idx4=0;
for i=1:length(IDs)
    id=IDs(i);
    name=Names{i};
    if i<0
        filename=['VR_CueCombination_2environment_doubledTrials_Response' num2str(id) '_' name '.dat'];
    else
        filename=['VR_CueCombination_2environmentFlag_doubledTrials_Response' num2str(id) '_' name '.dat'];
    end
    % filename=['VR_desktop_ProximalDistal_retrieve_' num2str(id) '_' name '.dat'];
    data=textread(filename);
    gender(i,1)=data(1,2);
    day=data(:,3);
    index_day=find( day==1 | day==2);
    environment=data(index_day,4);  % 1-rich; 2-poor
    trialNum=data(index_day,5);
    
    condition=data(index_day,7);
    post1=data(index_day,8);
    post2=data(index_day,9);
    post3=data(index_day,10);
    
    targx=data(index_day,11);
    targz=data(index_day,12);
    jitterx=data(index_day,13);
    jitterz=data(index_day,14);
    
    startx=data(index_day,15);
    startz=data(index_day,16);
    respx_org=data(index_day,17);% before coordinate transformation
    respz_org=data(index_day,18); % before coordinate transformation
    rt=data(index_day,19);
    
    
    respx=respx_org; % no transformation
    respz=respz_org;
    %[respx,respz]=CoordinateTransform(respx_org,respz_org,targx+jitterx,targz+jitterz);
    
%     disterr_all=sqrt( respx.^2+respz.^2  );
%     for n=1:2
%         for k=1:4
%             index=find(condition==k & environment==n);
            
    
    
    for n=1:2
        for k=1:4
            % % %             % % environments analyzed separately
            con_index_left{n,k}=find(condition==k & environment==n &targx<0);
            con_index_right{n,k}=find(condition==k & environment==n & targx>0);
            resp_time(i,k,n)=mean(rt(   [con_index_left{n,k};con_index_right{n,k}] ));    
            %
            %              % consider each side separately when identify outliers
            %             [response_left{n,k},delete_num_left{i,k,n}]=BoxPlotOutlier([respx(con_index_left{n,k}),respz(con_index_left{n,k})],quar);
            %             [response_right{n,k},delete_num_right{i,k,n}]=BoxPlotOutlier([respx(con_index_right{n,k}),respz(con_index_right{n,k})],quar);
            
            %              %consider both sides when identify outliers
            
            [response0_left{n,k},response0_right{n,k},delete_num_left{i,k,n},delete_num_right{i,k,n}]=BoxPlotOutlier2sides([respx(con_index_left{n,k}),respz(con_index_left{n,k})],[respx(con_index_right{n,k}),respz(con_index_right{n,k})],quar);
          
            % calculate accuracy
            temp1=response0_left{n,k};
            temp2=response0_right{n,k};
            [r1,c1]=size(temp1);
            [r2,c2]=size(temp2);
           disterr(i,k,n)=( sum(sqrt (temp1(:,1).^2+temp1(:,2).^2))+sum(sqrt (temp2(:,1).^2+temp2(:,2).^2)))/(r1+r2);
            
            
            if dimen==1
                response_left{n,k}=response0_left{n,k};
                response_right{n,k}=response0_right{n,k};
            elseif dimen==2
                response_left{n,k}=response0_left{n,k}(:,1);
                response_right{n,k}=response0_right{n,k}(:,1);
            elseif dimen==3
                response_left{n,k}=response0_left{n,k}(:,2);
                response_right{n,k}=response0_right{n,k}(:,2);
            end
%                                     %%%%%% save file for individual subjects for collapsed analysis
%                                     if dimen==1
%                                     savefile=['sub' num2str(id) '_n' num2str(n), '_k' num2str(k) '.mat'];
%                                     elseif dimen==2
%                                        savefile=['x_sub' num2str(id) '_n' num2str(n), '_k' num2str(k) '.mat'];
%                                     elseif dimen==3
%                                         savefile=['z_sub' num2str(id) '_n' num2str(n), '_k' num2str(k) '.mat'];
%                                     end
%                                     responses_save_left=response_left{n,k};
%                                     responses_save_right=response_right{n,k};
%                                     targleft=[targx(con_index_left{n,k})+jitterx(con_index_left{n,k}),targz(con_index_left{n,k})+jitterz(con_index_left{n,k})];
%                                     targright=[targx(con_index_right{n,k})+jitterx(con_index_right{n,k}),targz(con_index_right{n,k})+jitterz(con_index_right{n,k})];
%                                     save(savefile,'responses_save_left','responses_save_right','targleft','targright');
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [r1,c1]=size(response_left{n,k});
            [r2,c2]=size(response_right{n,k});
            
            trial_num(i,k,n)=r1+r2;
            trial_del(i,k,n)=delete_num_left{i,k,n}+delete_num_right{i,k,n};
            
%             cen_left{i,k,n}=mean(response_left{n,k});
%             cen_right{i,k,n}=mean(response_right{n,k});
            
            cen_left{k,n}(i,:)=mean(response_left{n,k});
            cen_right{k,n}(i,:)=mean(response_right{n,k});
            
            centers_left{n,k}=ones(length(con_index_left{n,k}),1)*mean(response_left{n,k});
            centers_right{n,k}=ones(length(con_index_right{n,k}),1)*mean(response_right{n,k});
            
            
            
            if dimen==1  %correct intrinsic bias in each side before the pooling
                [respx_left,respz_left]=CoordinateTransform(response_left{n,k}(:,1),response_left{n,k}(:,2),centers_left{n,k}(:,1),centers_left{n,k}(:,2));
                [respx_right,respz_right]=CoordinateTransform(response_right{n,k}(:,1),response_right{n,k}(:,2),centers_right{n,k}(:,1),centers_right{n,k}(:,2));
            elseif dimen==2 || dimen==3
                resp_left= response_left{n,k}-mean(response_left{n,k});
                resp_right= response_right{n,k}-mean(response_right{n,k});
            end
            
            if dimen==1
                response{n,k}=[respx_left,respz_left;respx_right,respz_right];
                centers{n,k}=mean(response{n,k});
                resp_var(i,k,n)=variance2d(response{n,k});
                resp_var_left(i,k,n)=variance2d(response_left{n,k});
                resp_var_right(i,k,n)=variance2d(response_right{n,k});
            elseif dimen==2 || dimen==3
                response{n,k}=[resp_left;resp_right];
                centers{n,k}=mean(response{n,k});
                resp_var(i,k,n)=var(response{n,k});
                resp_var_left(i,k,n)=var(response_left{n,k});
                resp_var_right(i,k,n)=var(response_right{n,k});
            end
            
            if n==2
            centers_pooled_left(i,:,k)=centers_left{n,k}(1,:);
            centers_pooled_right(i,:,k)=centers_right{n,k}(1,:);
            end
                   
            resp_std(i,k,n)=sqrt(resp_var(i,k,n));
            resp_std_left(i,k,n)=sqrt(resp_var_left(i,k,n));
            resp_std_right(i,k,n)=sqrt(resp_var_right(i,k,n));
            
            %%%%%%%%% calculate bias
             if dimen==1
                a=mean(response_left{n,k});
                bias_left{n}(i,k)=sqrt(a(1)^2+a(2)^2);
               a=mean(response_right{n,k});
                bias_right{n}(i,k)=sqrt(a(1)^2+a(2)^2);
                bias{n}(i,k)=(bias_left{n}(i,k)+bias_right{n}(i,k))/2;
            elseif dimen==2
                a=mean(response_left{n,k});
                bias_left{n}(i,k)=abs(a(1));
                a=mean(response_right{n,k});
                bias_right{n}(i,k)=abs(a(1));
                bias{n}(i,k)=(bias_left{n}(i,k)+bias_right{n}(i,k))/2;
            elseif dimen==3
                a=mean(response_left{n,k});
                bias_left{n}(i,k)=abs(a(1));
                a=mean(response_right{n,k});
                bias_right{n}(i,k)=abs(a(1));
                bias{n}(i,k)=(bias_left{n}(i,k)+bias_right{n}(i,k))/2;
            end
            
            %               %%%%%%%% save file for temporal analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %              savefile=['temporal_sub' num2str(id) '_n' num2str(n), '_k' num2str(k) '.mat'];
            %             responses_save=response{n,k};
            %             save(savefile,'responses_save');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        Baypred_var(i,n)=resp_var(i,1,n)*resp_var(i,2,n)/(resp_var(i,1,n)+resp_var(i,2,n));
        Baypred_std(i,n)=sqrt(Baypred_var(i,n));
        p_bay(i,n)=resp_var(i,2,n)/(resp_var(i,2,n)+resp_var(i,1,n));
        p_bay_left(i,n)=resp_var_left(i,2,n)/(resp_var_left(i,2,n)+resp_var_left(i,1,n));
        p_bay_right(i,n)=resp_var_right(i,2,n)/(resp_var_right(i,2,n)+resp_var_right(i,1,n));
        
        
    end
    
    
    
    for n=1:2
        % calculate landmark proximity in combination condition
        prox_comb_left(i,n)=RelativeProximityDist(centers_left{n,3}(1,:),centers_left{n,1}(1,:),centers_left{n,2}(1,:));
        prox_comb_right(i,n)=RelativeProximityDist(centers_right{n,3}(1,:),centers_right{n,1}(1,:),centers_right{n,2}(1,:));
        prox_comb(i,n)=(prox_comb_left(i,n)+prox_comb_right(i,n))/2;
        resp_proxinte_comb(i,n)=sqrt( resp_var(i,1,n)*prox_comb(i,n)^2+resp_var(i,2,n)*(1-prox_comb(i,n))^2);
        
        resp_proxalt_comb(i,n)=alternation_model_normal_2sides(prox_comb_left(i,n),centers_left{n,1}(1,:),centers_left{n,2}(1,:),resp_std_left(i,1,n),...
            resp_std_left(i,2,n),prox_comb_right(i,n),centers_right{n,1}(1,:),centers_right{n,2}(1,:),resp_std_right(i,1,n),...
            resp_std_right(i,2,n),length(con_index_left{n,3}),length(con_index_right{n,3}));
        resp_bayalt_comb(i,n)=alternation_model_normal_2sides(p_bay_left(i,n),centers_left{n,1}(1,:),centers_left{n,2}(1,:),resp_std_left(i,1,n),...
            resp_std_left(i,2,n),p_bay_right(i,n),centers_right{n,1}(1,:),centers_right{n,2}(1,:),resp_std_right(i,1,n),...
            resp_std_right(i,2,n),length(con_index_left{n,3}),length(con_index_right{n,3}));
        
        
        % calculate landmark proximity by shifting the vision center to be aligned
        % with the landmark-defined target location in the conflict condition
        
        targ_lm_left{n}=landmarkdefinedtarget([targx(con_index_left{n,4})+jitterx(con_index_left{n,4}),targz(con_index_left{n,4})+jitterz(con_index_left{n,4})],15,dimen);   %find the landmark defined target locations, which are different on different trials
        targ_lm_right{n}=landmarkdefinedtarget([targx(con_index_right{n,4})+jitterx(con_index_right{n,4}),targz(con_index_right{n,4})+jitterz(con_index_right{n,4})],15,dimen);   %find the landmark defined target locations, which are different on different trials
        
        targ_lm_left_out{n}(i,:)=mean(targ_lm_left{n});
        targ_lm_right_out{n}(i,:)=mean(targ_lm_right{n});
        
        targ_lm{n}=[targ_lm_left{n};targ_lm_right{n}];
        
        prox_conf_left(i,n)=RelativeProximityDist(centers_left{n,4}(1,:),centers_left{n,1}(1,:)+mean(targ_lm_left{n}),centers_left{n,2}(1,:));
        prox_conf_right(i,n)=RelativeProximityDist(centers_right{n,4}(1,:),centers_right{n,1}(1,:)+mean(targ_lm_right{n}),centers_right{n,2}(1,:));
        prox_conf(i,n)=(prox_conf_left(i,n)+prox_conf_right(i,n))/2;
        resp_proxinte_conf(i,n)=sqrt( resp_var(i,1,n)*prox_conf(i,n)^2+resp_var(i,2,n)*(1-prox_conf(i,n))^2);
        
        
        
        resp_proxalt_conf(i,n)=alternation_model_normal_2sides(prox_conf_left(i,n),centers_left{n,1}(1,:)+mean(targ_lm_left{n}),centers_left{n,2}(1,:),resp_std_left(i,1,n),...
            resp_std_left(i,2,n),prox_conf_right(i,n),centers_right{n,1}(1,:)+mean(targ_lm_right{n}),centers_right{n,2}(1,:),resp_std_right(i,1,n),...
            resp_std_right(i,2,n),length(con_index_left{n,4}),length(con_index_right{n,4}));
        resp_bayalt_conf(i,n)=alternation_model_normal_2sides(p_bay_left(i,n),centers_left{n,1}(1,:)+mean(targ_lm_left{n}),centers_left{n,2}(1,:),resp_std_left(i,1,n),...
            resp_std_left(i,2,n),p_bay_right(i,n),centers_right{n,1}(1,:)+mean(targ_lm_right{n}),centers_right{n,2}(1,:),resp_std_right(i,1,n),...
            resp_std_right(i,2,n),length(con_index_left{n,4}),length(con_index_right{n,4}));
        
        
        
        VARindex_comb(i,1,n)= (1/resp_var(i,3,n))  / (1/resp_var(i,1,n)+1/resp_var(i,2,n)); % bayinteg_
        VARindex_comb(i,2,n)=(1/resp_var(i,3,n)) / (1/resp_bayalt_comb(i,n)^2); %bayalt_
        VARindex_comb(i,3,n)= (1/resp_var(i,3,n)) / (1/resp_proxinte_comb(i,n)^2); %proxinteg_
        VARindex_comb(i,4,n)=(1/resp_var(i,3,n)) / (1/resp_proxalt_comb(i,n)^2); %proxalt_
        
        
        VARindex_conf(i,1,n)= (1/resp_var(i,4,n))  / (1/resp_var(i,1,n)+1/resp_var(i,2,n)); %bayinteg_
        VARindex_conf(i,2,n)=(1/resp_var(i,4,n)) / (1/resp_bayalt_conf(i,n)^2); %bayalt_
        VARindex_conf(i,3,n)= (1/resp_var(i,4,n)) / (1/resp_proxinte_conf(i,n)^2); %proxinteg_
        VARindex_conf(i,4,n)=(1/resp_var(i,4,n)) / (1/resp_proxalt_conf(i,n)^2); %proxalt_
        
        
        
        RPindex(i,1,n)=prox_comb(i,n)/p_bay(i,n); % _comb
        RPindex(i,2,n)=prox_conf(i,n)/p_bay(i,n); % _conf
        RPindex(i,3,n)=(prox_comb(i,n)/2+prox_conf(i,n)/2)/p_bay(i,n); %_average
    end
    
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
%
% savefile=['rich_poor_2sides_' num2str(quar) '.mat'];
%
% %save(savefile,'resp_std','Bay_intepred_std','prox_comb','prox_conf','rr');
% save(savefile,'resp_std','Baypred_std', 'resp_bayalt_comb', 'resp_bayalt_conf','resp_proxinte_comb','resp_proxinte_conf',...
%     'resp_proxalt_comb','resp_proxalt_conf','prox_comb','prox_conf','p_bay');

dataouput(:,:,1)=[resp_std(:,:,1),Baypred_std(:,1),resp_bayalt_comb(:,1),resp_bayalt_conf(:,1),resp_proxinte_comb(:,1),resp_proxinte_conf(:,1),resp_proxalt_comb(:,1),resp_proxalt_conf(:,1),...
    prox_comb(:,1),prox_conf(:,1),p_bay(:,1)];
dataouput(:,:,2)=[resp_std(:,:,2),Baypred_std(:,2),resp_bayalt_comb(:,2),resp_bayalt_conf(:,2),resp_proxinte_comb(:,2),resp_proxinte_conf(:,2),resp_proxalt_comb(:,2),resp_proxalt_conf(:,2),...
    prox_comb(:,2),prox_conf(:,2),p_bay(:,2)];

% t test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% variability %%%%%%%%%%%%%%%%%%%%%%
std_all(:,:,1)=[resp_std(:,:,1),Baypred_std(:,1),resp_bayalt_comb(:,1),resp_bayalt_conf(:,1),resp_proxinte_comb(:,1),...
    resp_proxinte_conf(:,1),resp_proxalt_comb(:,1),resp_proxalt_conf(:,1)];
std_all(:,:,2)=[resp_std(:,:,2),Baypred_std(:,2),resp_bayalt_comb(:,2),resp_bayalt_conf(:,2),resp_proxinte_comb(:,2),...
    resp_proxinte_conf(:,2),resp_proxalt_comb(:,2),resp_proxalt_conf(:,2)];
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

for n=1:2  % n=1, no rotation; n=2, rotation
    index=0;
    for i=1:12
        index=index+1;
        [h(index,n) p(index,n) ci{index} stats]=ttest(std_all(:,pair{i}(1),n),std_all(:,pair{i}(2),n));  % 1 vision vs comb
        t(index,n)=stats.tstat;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% proximity vs rr %%%%%%%%%%%%%%%%%%%%%%%%
prox_all(:,:,1)=[prox_comb(:,1),prox_conf(:,1),(prox_comb(:,1)+prox_conf(:,1))/2,p_bay(:,1)];
prox_all(:,:,2)=[prox_comb(:,2),prox_conf(:,2),(prox_comb(:,2)+prox_conf(:,2))/2,p_bay(:,2)];
pair{1}=[1 4]; % comb vs rr
pair{2}=[2 4]; % conf vs rr
pair{3}=[3 4]; % comb&conf vs rr
for n=1:2
    index=0;
    for i=1:3
        index=index+1;
        [h(index,n) p_prox(index,n) ci{index} stats]=ttest(prox_all(:,pair{i}(1),n),prox_all(:,pair{i}(2),n));  % 7 comb vs bay integration prediction
        t_prox(index,n)=stats.tstat;
    end
end
%correlation between proximity and rr
pair{1}=[1 4];
pair{2}=[2 4];
pair{3}=[3 4];
index=0;
for n=1:2
    index=0;
    for i=1:3
        index=index+1;
        [r0,p0]=corrcoef(prox_all(:,pair{i}(1),n),prox_all(:,pair{i}(2),n));
        r(index,n)=r0(1,2);
        p_r(index,n)=p0(1,2);
    end
end

% manipulation check
index=1;
[h_com p_com(index,1) ci_com stats]=ttest(std_all(:,1,1),std_all(:,1,2));
t_com(index,1)=stats.tstat;
index=index+1;
[h_com p_com(index,1) ci_com stats]=ttest(std_all(:,2,1),std_all(:,2,2));
t_com(index,1)=stats.tstat;
index=index+1;
[h_com p_com(index,1) ci_com stats]=ttest(prox_all(:,4,1),prox_all(:,4,2));
t_com(index,1)=stats.tstat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% compare variability bayesian index (bayesian integration index & bayesian alternation index) to 1
for n=1:2
    VARindex(:,1,n)=VARindex_comb(:,1,n);  % combination, bayesian integration
    VARindex(:,2,n)=VARindex_comb(:,2,n);  % combination, bayesian alternation
    VARindex(:,3,n)=VARindex_comb(:,3,n);  % combination, proximity integration
    VARindex(:,4,n)=VARindex_comb(:,4,n);  % combination, proximity alternation
        
    VARindex(:,5,n)=VARindex_conf(:,1,n);  % conflict, bayesian integration
    VARindex(:,6,n)=VARindex_conf(:,2,n);   % conflict, bayesian alternation
    VARindex(:,7,n)=VARindex_conf(:,3,n);  % conflict, proximity integration
    VARindex(:,8,n)=VARindex_conf(:,4,n);   % conflict, proximity alternation
end

%bayesian integration,combination
for i=1:8
    for n=1:2
        [mean_ind(i,n),Nsub(i,n),p_ind(i,n),t_ind(i,n)]=IndexTest(VARindex(:,i,n),1);
    end
    
end


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



