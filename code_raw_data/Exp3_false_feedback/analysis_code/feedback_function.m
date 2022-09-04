

clear all;
% negative feedback


e=2.71828;
%theta=4;
theta=2.5;
%leta=1;
leta=0.9;
x1=0:0.01:5;
figure(1)

y=zeros(1,length(x1));  % negative feedback
for i=1:length(x1)
    y(i)=theta*(1-e ^ (-leta*x1(i)));
end

figure(1)
plot(x1,y);
hold on
plot(x1,x1,'k-');

% positive feedback
x2=0:0.01:theta;
z=zeros(1,length(x2));
for i=1:length(x2)
    z(i)=(-1/leta) * log (1-x2(i)/theta);
end

hold on

plot(x2,z,'r-')

%find the cut point

% for positive feedback
hold on;
z_lin=x1 .*.8;
plot(x1,z_lin,'r-')

% for negative feedback
hold on;
y_lin=x1 .* 1.2;
plot(x1,y_lin,'b-')
