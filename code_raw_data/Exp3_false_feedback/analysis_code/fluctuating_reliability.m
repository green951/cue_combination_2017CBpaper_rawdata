
% changing reliability
clear all;
close all;
%
v_mu=0;  % 0, no bias
s_mu=0;  % 0, no bias
v_std0=1;  % initially equivalent variability between two cues
s_std0=1;

% evaluate cue reliability using the instantenous estimation
% reliability function
%rates=0:0.1:2;
rate=0.1;
iterations=10;

Ntrials=[100,300,500,700,1000];

wind_len=[1];

for q=1:iterations
    

    for n=1:length(Ntrials)
      Ntrial=Ntrials(n);
        v=zeros(1,Ntrial);
        s=zeros(1,Ntrial);
        comb=zeros(1,Ntrial);
        conf=zeros(1,Ntrial);
        
        for i=1:Ntrial,
            v_std(i,n)=v_std0 *(i^ (-rate/1));  % vision positive
            %v_std(i,n)=v_std0 * (1- i/Ntrial);
            
            if v_std(i,n)<.2
                v_std(i,n)=.2;
            end
            
            s_std(i,n)=s_std0;
            if i>2
             w_opt(i,n)=var(s(1:i-1))/(var(s(1:i-1))+var(v(1:i-1)));
            else
                w_opt(i,n)=s_std(i,n)^2/(s_std(i,n)^2+v_std(i,n)^2);
            end
%             
%             if i>1+wind_len(n)
%                 w_opt(i,n)=var(s(i-1-wind_len(n):i-1))/(var(s(i-1-wind_len(n):i-1))+var(v(i-1-wind_len(n):i-1)));
%             elseif i>2
%                 w_opt(i,n)=var(s(1:i-1))/(var(s(1:i-1))+var(v(1:i-1)));
%             else
%                 w_opt(i,n)=s_std(i,n)^2/(s_std(i,n)^2+v_std(i,n)^2);
%             end

            v(i)=randn(1,1)*v_std(i,n)+v_mu;
            s(i)=randn(1,1)*s_std(i,n)+s_mu;
            
            comb(i)=v(i)*w_opt(i,n)+s(i)*(1-w_opt(i,n));
            conf(i)=v(i)*w_opt(i,n)+s(i)*(1-w_opt(i,n));
        end
        
        std_opt(q,n)=sqrt(var(v)*var(s)/(var(v)+var(s)));
        std_act(q,n)=std(comb);
        pred_err(q,n)=std_opt(n)/std_act(n);
    end
    
    
    
    
end


figure(1)
plot(wind_len,mean(pred_err));
xlabel('integration window')
ylabel('prediction error')
%
figure(2)
plot(1:1:Ntrial, v_std);
xlabel('trial number')
ylabel('instantenous reliability')


