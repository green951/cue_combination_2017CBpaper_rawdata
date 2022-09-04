% correlations between rr and rp

weight=[0.58 	0.43 	0.70 	0.57 
0.74 	0.80 	0.76 	0.78 		
0.73 	0.36 	0.37 	0.37 
0.97 	0.88 	0.89 	0.89 
0.79 	0.65 	0.66 	0.66 		
0.51 	0.70 	0.72 	0.71 
0.70 	0.75 	0.61 	0.68 
0.87 	0.71 	0.82 	0.76 
0.57 	0.65 	0.51 	0.58 
0.73 	0.76 	0.84 	0.80 
0.85 	0.61 	0.77 	0.69 
0.63 	0.72 	0.58 	0.65 
0.49 	0.75 	0.61 	0.68 
0.98 	0.97 	0.96 	0.97 
0.72 	0.76 	0.85 	0.80 
0.82 	0.66 	0.64 	0.65 
0.82 	0.74 	0.83 	0.78 
0.74 	0.81 	0.84 	0.83 
];  %  1st column, rr; 2nd column, rp_comb; 3rd column, rp_conf; 4th column, rp_averaged

condition=input('Please enter the condition (1, comb; 2, conf;3, averaged)');

if condition==1  % rp, dependent variable
    y=weight(:,2);
elseif condition==2
    y=weight(:,3);
elseif condition==3
    y=weight(:,4);
end

x=weight(:,1);  % rr, predictor
X=[x,ones(length(x),1)];
a=regress(y,X);
y_pred=X*a;
y_res=y-y_pred;

boxplot(y_res,'whisker',3)




        
        