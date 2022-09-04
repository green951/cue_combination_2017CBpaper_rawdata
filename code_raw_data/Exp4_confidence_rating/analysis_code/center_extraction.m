%a=cen_left;
a=cen_right;

% cen_left{i,n,k}
     
%  c_rich=cen_left{:,:,1};
%  c_poor=cen_left{:,:,2};
        
  for n=1:2   
      for k=1:4
          tempp_left=[0,0];
          tempp_right=[0,0];
          for i=1:18 
              tempp_left=tempp_left + a{i,k,n};
              tempp_right=tempp_right + a{i,k,n};
          end
          cen_output_left(k,:,n)=tempp_left/18;
          cen_output_right(k,:,n)=tempp_right/18;
      end
  end 
 
            
            
            
        