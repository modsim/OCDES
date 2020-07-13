function demo_poly
% Solve the OCDE in form of 
% dx=f(x,v), x=x0
% min_v obj(x,v), s.t. heq(x,v)=0, hieq(x,v)>=0

%demo example:
%dx=Ax+Bv
%min 1/2 vT*M*v
%s.t. L*z >= Kx

%This example fulfil the sufficient condition of reviewer. It has a
%solution

clear;
clc;
close all;
currentFolder=pwd;
cd ..\
addpath(strcat(pwd,'\OCDEtool'));
cd(currentFolder);





%%Start: Defined by user
syms T p4 z
p1=1;
p2=0.5;
p3=-2.5;

x=[T;p4];
v=z;

f=[1;
   -2/(exp(10*T-5)+exp(-10*T+5))^2];
obj=p1*z^4+p2*z^3+p3*z^2+p4*z; %excellular equation, dx=f(x,v)
obj_fun=@(p4,zv)p1.*zv.^4+p2.*zv.^3+p3.*zv.^2+p4.*zv;
heq=[]; %equality constraint, heq=0
hieq=[z+1.5;
      1.5-z];% inequality constraints, hieq>=0
x0=[0;-4.5398e-05]; %initial condition of x
v0=[-1]; %initial guess of v0
opt_init.tol_act=1e-4; %tolerance for active inequality constraint
opt_init.optimoptions=optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500);
opt_sol.integrator='ode15s';
opt_sol.opt_integrator=odeset('RelTol', 1e-10, 'AbsTol', 1e-10,'InitialStep',1e-14); %option for integrator
opt_sol.tol_feasible=1e-4; %feasibility tolerence
opt_sol.MaxNoUptActiveSet=100; %maximum number of updating active set
tstart=0;
tfinal=1;
%End: Defined by user

%%Main part
[tout1,yout1,teout,yeout,ieq_actout,miu_ieq, miu_eq]=sOCDE_main(f,obj,heq, hieq,x,x0,...
    v,v0,tstart, tfinal, opt_init,opt_sol);
gv_out1=obj_fun(p4_fun(tout1),yout1(:,3));
v0=[1];
[tout2,yout2,teout,yeout,ieq_actout,miu_ieq, miu_eq]=sOCDE_main(f,obj,heq, hieq,x,x0,...
    v,v0,tstart, tfinal, opt_init,opt_sol);
gv_out2=obj_fun(p4_fun(tout2),yout2(:,3));



%Plot
figure
t_step=0:0.05:1;
z_step=-1.5:0.05:1.5;
[t_mesh,z_mesh]=meshgrid(t_step,z_step);
gv=zeros(length(z_step),length(t_step));
for i=1:size(t_mesh,1) 
    for j=1:size(t_mesh,2)
    pv=p4_fun(t_mesh(i,j));
    zv=z_mesh(i,j);
    gv(i,j)=p1*zv^4+p2*zv^3+p3*zv^2+pv*zv;
    end
end
surface(t_mesh,z_mesh,gv);
xlabel('t'); ylabel('z'); zlabel('objective function')
hold on
plot3(tout1,yout1(:,3),gv_out1,'Linewidth',6,'Color','red')
plot3(tout2,yout2(:,3),gv_out2,'Linewidth',6,'Color','blue')
set(gca,'FontSize', 14)
end
    
function y=p4_fun(t)
    y=0.5*(1-tanh(10*(t-0.5)))-1;
end