function demo2
% Solve the OCDE in form of 
% dx=f(x,v), x=x0
% min_v obj(x,v), s.t. heq(x,v)=0, hieq(x,v)>=0

%demo example:


clear;
clc;
currentFolder=pwd;
cd ..\

addpath(strcat(pwd,'\OCDEtool'));
cd(currentFolder);

%%Start: Defined by user
syms x1 x2  z1 z2
x=[x1;x2]; %state of ODE
v=[z1;z2]; % var. for inner opt.

f=[sin(z1*z2)+1;
   cos(x1*z1+z2)]; %excellular equation, dx=f(x,v)
obj=v'*v; %objective function, sym 
heq=[]; %equality constraint, heq=0
hieq=[1.8*sin(x1)+z2-exp(-z1);
      0.2*cos(x2)+2-z2;
      ]; % inequality constraints, hieq>=0
x0=[pi/4;0]; %initial condition of x
%x0=[3.9521;11.1138];
z0=[]; %initial guess of v0
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
    v,z0,tstart, tfinal, opt_init,opt_sol);


%%Plot simulation results
figure
h1=plot(tout,yout(:,1));
hold on
h2=plot(tout,yout(:,2),'--');
l=legend('x_1(t)','x_2(t)');
xlabel('t');
plot(teout(1:end-1),yeout(1:end-1,1),'o');
plot(teout(1:end-1),yeout(1:end-1,2),'o');
set(l,'FontSize',14);
set(gca,'fontsize',14);

figure
h1=plot(tout,yout(:,3));
hold on
h2=plot(tout,yout(:,4),'--');
l=legend('z_1(t)','z_2(t)');
xlabel('t');
plot(teout(1:end-1),yeout(1:end-1,3),'o');
plot(teout(1:end-1),yeout(1:end-1,4),'o');
set(l,'FontSize',14);
set(gca,'fontsize',14);


figure
h1=plot(tout,yout(:,5));
hold on
h2=plot(tout,yout(:,6),'--');
l=legend('\mu_1 (t)','\mu_2(t)')
xlabel('t');
plot(teout(1:end-1),yeout(1:end-1,5),'o');
plot(teout(1:end-1),yeout(1:end-1,6),'o');
set(l,'FontSize',14);
set(gca,'fontsize',14);

figure
matlabFunction(hieq,'File','hieq_gen','Vars',{[x;v;miu_ieq;miu_eq]});
hieq_fun=@(allvar)hieq_gen(allvar);
hieq_v=[];
for i=1:length(tout)
hieq_v(:,i)=hieq_fun(yout(i,:)');
end
h1=plot(tout,hieq_v(1,:));
hold on
h2=plot(tout,hieq_v(2,:),'--');
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
l=legend('l_1(y(t))','l_2(y(t))');
xlabel('t');
set(l,'FontSize',14);
set(gca,'fontsize',14);


end