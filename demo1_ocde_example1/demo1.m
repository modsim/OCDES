function demo1
% Solve the OCDE in form of 
% dx=f(x,v), x=x0
% min_v obj(x,v), s.t. heq(x,v)=0, hieq(x,v)>=0

%demo example:
%dx1=0.5+v1*v2, x1(0)=pi/4
%dx2=x1, x2(0)=0
%min v1^2+v2^2
%s.t. 1.7*sin(x1)+v2-exp(-v1) >=0,
%     0.2*cos(x2)+2-v2>=0.

clear;
clc;
currentFolder=pwd;
cd ..\ 

addpath(strcat(pwd,'\OCDEtool'));
cd(currentFolder);

%%Start: Defined by user
syms x1 x2  v1 v2
x=[x1;x2]; %state of ODE
v=[v1;v2]; % var. for inner opt.

f=[0.5+v1*v2;
    x1]; %excellular equation, dx=f(x,v)
obj=v'*v; %objective function, sym 
heq=[]; %equality constraint, heq=0
hieq=[sin(x1)*1.7+v2-exp(-v1);
      0.2*cos(x2)+2-v2;
      ]; % inequality constraints, hieq>=0
x0=[pi/4;0]; %initial condition of x
%x0=[3.9521;11.1138];
v0=[]; %initial guess of v0
opt_init.tol_act=1e-4; %tolerance for active inequality constraint
opt_init.optimoptions=optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500);
opt_sol.integrator='ode15s';
opt_sol.opt_integrator=odeset('RelTol', 1e-10, 'AbsTol', 1e-10,'InitialStep',1e-14); %option for integrator
opt_sol.tol_feasible=1e-4; %feasibility tolerence
opt_sol.MaxNoUptActiveSet=100; %maximum number of updating active set
tstart=0;
tfinal=20;
%End: Defined by user

%%Main part
[tout,yout,teout,yeout,ieq_actout,miu_ieq, miu_eq]=sOCDE_main(f,obj,heq, hieq,x,x0,...
    v,v0,tstart, tfinal, opt_init,opt_sol);


%%Plot simulation results
figure
h1=plot(tout,yout(:,1),'red');
hold on
h2=plot(tout,yout(:,2),'blue');
legend('x_1(t)','x_2(t)')
xlabel('t');
plot(teout(1:end-1),yeout(1:end-1,1),'o');
plot(teout(1:end-1),yeout(1:end-1,2),'o');

figure
h1=plot(tout,yout(:,3),'red');
hold on
h2=plot(tout,yout(:,4),'blue');
legend('v_1(t)','v_2(t)')
xlabel('t');
plot(teout(1:end-1),yeout(1:end-1,3),'o');
plot(teout(1:end-1),yeout(1:end-1,4),'o');

figure
h1=plot(tout,yout(:,5),'red');
hold on
h2=plot(tout,yout(:,6),'blue');
legend('\mu_1 (t)','\mu_2(t)')
xlabel('t');
plot(teout(1:end-1),yeout(1:end-1,5),'o');
plot(teout(1:end-1),yeout(1:end-1,6),'o');


figure
matlabFunction(hieq,'File','hieq_gen','Vars',{[x;v;miu_ieq;miu_eq]});
hieq_fun=@(allvar)hieq_gen(allvar);
hieq_v=[];
for i=1:length(tout)
hieq_v(:,i)=hieq_fun(yout(i,:)');
end
h1=plot(tout,hieq_v(1,:),'red');
hold on
h2=plot(tout,hieq_v(2,:),'blue');
hieq_e=[];
for i=1:length(teout)
hieq_e(:,i)=hieq_fun(yeout(i,:)');
end
for i=1:length(teout)-1 
    if abs(hieq_e(1,i))<1e-4
plot(teout(i),0,'o');
    end
    if abs(hieq_e(2,i))<1e-4
plot(teout(i),hieq_e(2,i),'o');
    end    
end
legend('hieq(1)(t)','hieq(2)(t)')
xlabel('t');



end

% function dx=f_demo1(x,v)
% x1=v(1);
% x2=v(2);
% dx=[0.5+x1*x2;
%     x(1)];
% end