function [tout,yout,teout,yeout,ieq_actout]=OCDE_sol(ODE_fun, ...
    KKT_fun,x0, v0, miu_ieq0, miu_eq0, ieq_act0, tstart, tfinal, ...
    event_fun,Jac_fun,option)
%This function integrates the OCDE by using KKT optimality conditions. It %solves the following quasi-DAE system
%dx=ODE_fun(x,v);
%0=KKT_fun([x;v;miu_ieq;miu_eq], ieq_act)
%x=x, z=[v;miu_ieq;miu_eq], y=[x;z]

%Inputs:
%ODE_fun: the right hand side of dx=f(x,v) (function handle)
%KKT_fun([x;v;miu_ieq;miu_eq], ieq_act): return the right hand side of 0=KKT. %(function handle).
%x0:initial value of x (real vector) 
%v0: initial value of v, i.e. the optimal solution of inner NLP (real vector)
%miu_ieq0: initial value of Lagrange multipliers miu_ieq (real vector)
%miu_eq0: initial value of Lagrange multipliers miu_ieq (real vector)
%ieq_act0: initial value of the active set index ieq_act (real vector)
%tstart: start time of integration (real scalar)
%tfinal: end time of integeration (real scalar)
%event_fun([x;v;miu_ieq;miu_eq]): function handle to detect the switching of %active set. It is define to be miu_ieq-hieq(x,v). (function handle)
%Jac_fun(t,[x;v;miu_ieq;miu_eq]), ieq_act): Jacobian matrix df(x,z)/dxz. If empty, it is approximated by finite difference methods.
%option: options

%Outputs:
%tout: integration steps of t (real vector)
%yout: integration results of [x;v;miu_ieq;miu_eq] (real matrix)
%teout: the time, when the switching of active set happens (real vector)
%yeout: states of [x;v;miu_ieq;miu_eq], when the switching of active set happens (real matrix)
%ieq_actout: record the different values of active sets along the solution trajectory (real matrix)


y0 = [x0;v0;miu_ieq0;miu_eq0];
n_y0=length(y0);
nx=length(x0);
nv=length(v0);
nmiu_ieq=length(miu_ieq0);
nmiu_eq=length(miu_eq0);
nz=n_y0-nx;

x_id=[1:nx];
v_id=[nx+1:nx+nv];

if nmiu_ieq>0
miu_ieq_id=[nx+nv+1:nx+nv+nmiu_ieq];
else
miu_ieq_id=[];
end
if nmiu_eq>0
miu_eq_id=[nx+nv+nmiu_ieq+1:nx+nv+nmiu_ieq+nmiu_eq];
else
miu_eq_id=[];    
end


f_fun=@(t,y, ieq_act)[ODE_fun(y(x_id),y(v_id));KKT_fun(y, ieq_act)];
Mass=zeros(n_y0,n_y0);
for i=1:nx
    Mass(i,i)=1;
end


if isempty(Jac_fun)    
    Jac_fun=@(t,y,ieq_act_s)gradient_FD(@(y)f_fun(t,y,ieq_act_s),y,1e-6);
end

refine = 4;  

if isfield(option, 'integrator') && ~isempty(option.integrator)
integrator=option.integrator;
else
integrator='ode15s';
end  

if isfield(option, 'opt_integrator') && ~isempty(option.opt_integrator)
option_ode = odeset('OutputFcn', @odeplot,...
                  'OutputSel',1,'Refine',refine, 'Mass', Mass,...
                  'MStateDependence','none', 'RelTol', 1e-6, 'AbsTol', 1e-6,...
                  'InitialStep',1e-14);
else
option_ode = odeset(option.opt_integrator, 'OutputFcn', @odeplot,...
                 'OutputSel',1,'Refine',refine, 'Mass', Mass,...
                 'MStateDependence','none');
end

if isfield(option, 'tol_feasible') && ~isempty(option.tol_feasible)
tol_feasible=option.tol_feasible;
else
tol_feasible=1e-6;
end    

if isfield(option, 'MaxNoUptActiveSet') && ~isempty(option.MaxNoUptActiveSet)
MaxNoUptActiveSet=option.MaxNoUptActiveSet;
else
MaxNoUptActiveSet=100;
end  

tout = tstart;
yout = y0';
teout = [tstart];
yeout = [y0'];
te=tstart;
ye=y0'; 
  
  
ieq_actout=[ieq_act0];

ieq_act_s=ieq_act0;
for i = 1:MaxNoUptActiveSet
    
%Check LICQ of new system
  vJ=Jac_fun(te,ye',ieq_act_s);
  vJ_s=vJ(nx+1:end,nx+1:end);
  grad_heq=vJ_s(nv+1:nv+nmiu_eq,1:nv);
  grad_hieq_all=vJ_s(nv+nmiu_eq+1:end,1:nv);
  grad_hieq=[];
  for j=1:length(ieq_act_s)
      if ieq_act_s(j)==1
          grad_hieq=[grad_hieq;
                     grad_hieq_all(j,:);    
          ];
      end
  end
  grad=full([grad_heq;grad_hieq]);
  if rank(grad)<nmiu_eq+sum(ieq_act_s)
     warning('LICQ condition fails.\n');
     isLICQ=0;
  else
     isLICQ=1;     
  end    
    
    
  
% Check SOSC  
  ieq_act_1=ones(nmiu_ieq,1);
  vJ=Jac_fun(te,ye',ieq_act_1);
  D2Lag_vv=vJ(nx+1:nx+nv,nx+1:nx+nv); %Hessian matrix
  
  dheq=vJ(nx+nv+1:nx+nv+nmiu_eq,nx+1:nx+nv);
  dhieq=vJ(nx+nv+nmiu_eq+1:end,nx+1:nx+nv);
  miu_ieq_s=ye(nx+nv+1:nx+nv+nmiu_ieq);
  isSOSC=checkSOSC(D2Lag_vv,dheq,dhieq,ieq_act_s,miu_ieq_s);
  if ~isSOSC
  fprintf('SOSC fails at the checked point.\n');
  end
  
  %Check the index of DAE
  vJ=Jac_fun(te,ye',ieq_act_s);
  vJ_s=vJ(nx+1:end,nx+1:end);
  if isempty(null(vJ_s))
  isIndex1=1;
  else
  isIndex1=0;    
  fprintf('DAE index is higher than 1.\n');
  end
  
  
  %check consistant initial point
  fv=f_fun(te,ye',ieq_act_s);
  fv_alg=abs(fv(nx+1:end));
  if any(fv_alg>=1e-6)
     warning('Non-consistant initial point, norm(f_algebraic)=%e\n', norm(fv_alg));
  end
    
  if i>=2
  fprintf('No. of updating active set: %d\n',i-1);   
  end
  event_m_fun=@(t,y)events_mother(event_fun(y),ieq_act_s);
  option_ode = odeset(option_ode, 'Events', event_m_fun);
  option_ode = odeset(option_ode, 'Jacobian', @(t,y)Jac_fun(t,y,ieq_act_s)); 
    
  % Solve until the first terminal event.
  [t,y] =feval(integrator,@(t,y)f_fun(t,y,ieq_act_s),[tstart tfinal],y0,option_ode);
  te=t(end);
  ye=y(end,:);   
  hold on
  
  %check feasibility mu_j>=0, hieq>=0
  id_ifea=[];
  for j=1:length(t)
  d_mh=event_fun(y(j,:)');
  eps=tol_feasible;
      for k=1:length(d_mh)
          if ieq_act_s(k)==1 && d_mh(k)<=-eps 
             id_ifea=[id_ifea;[j,k]]; 
             warning('Infeasible mu at t(%d), mu_%d=%e\n',j,k,d_mh(k));             
          elseif ieq_act_s(k)==0 && d_mh(k)>=eps
             id_ifea=[id_ifea;[j,k]]; 
             warning('Infeasible hieq>=0 at t(%d), hieq_%d=%e\n',j, k, -d_mh(k));   
          end                    
      end  
  end
  
  % Accumulate output.  
  nt = length(t);
  tout = [tout; t(2:nt)];
  yout = [yout; y(2:nt,:)];
  
  % Set the new initial conditions, with .9 attenuation.
  y0 = y(end,:)';
  tstart = t(nt); 
  if t(end)>=tfinal
     break;
  end
 
  
  % Update the new active set
  ieq_act_o=ieq_act_s; %ieq_act_o: act ieq at the last phase
  value_s=event_m_fun(t(end-1),y(end-1,:)');
  value_e=event_m_fun(t(end),y(end,:)');
  Value=[value_s';value_e'];
  switch_j=[];
  for j=1:size(Value,2)
      v_st=Value(1,j);
      v_end=Value(end,j);
      isact_j=ieq_act_o(j);
      if v_st*v_end>0 %v_st and v_end do not change the sign
         continue;
      end
      
      
      if isact_j % if active, diff changes from + to -
           if v_st-v_end>0 
              ieq_act_s(j)=0;
              switch_j=[switch_j;j]; %we do not allow multiple switching
           else
              warning('Drection error.\n');
           end           
      else %if non-active, diff changes from - to +
           if v_st-v_end<0 
              ieq_act_s(j)=1;
              switch_j=[switch_j;j]; %we do not allow multiple switching
           else
              warning('Drection error.\n');
           end
      end    
  end  
  
  %store event time
  %store DAE(ieq_actout)
  teout = [teout; te];    % Events at tstart are never reported.
  yeout = [yeout; ye];
  %ieout = [ieout; ie];  
  ieq_actout=[ieq_actout,ieq_act_s];
  
  %Check wether there are multiple switch possiblities
  if length(switch_j)>=2
      warning('Multiple switch happens.\n');
  end
  
   
  %Check tranversaility condition
  [dx,dz]=getdxdz(@(t,y)f_fun(t,y,ieq_act_s),t(end),y0(1:nx),y0(nx+1:end),@(t,y)Jac_fun(t,y,ieq_act_s));
  for j=1:length(switch_j)
      sj=switch_j(j);
      if ieq_act_s(sj)==1 %activated ->check the sign of du_j/dt
         du_j=dz(nv+sj);
         if du_j<=0
            warning('Transversality condition failed. Continuing simulation leads to infeasibility of mu>=0.\n');
         end
      else  %deactivated ->check the sign of du_j/dt
         deventdy=gradient_FD(@(y)event_fun(y),y0,1e-6);
         deventdy_sj=deventdy(sj,:); %d(mu_sj-hieq_sj(x,z))/dy
         dhieq_dx=-deventdy_sj(1:nx);
         e=zeros(1,nz);
         e(nv+sj)=1;
         dhieq_dz=-(deventdy_sj(nx+1:end)-e);
         dldt=dhieq_dx*dx+dhieq_dz*dz;
         if dldt<=0
            warning('Transversality condition failed. Continuing simulation leads to infeasibility of l_j(x,z)>=0.\n');
         end
      end
  end
   
  
  
end

end

function [dx,dz]=getdxdz(F_fun,t0,x0,z0,Jac_fun)
%F_fun(t,[x0;z0])=[dx;0]=[f(.);g(.)]
%Jac_fun(t,[x0;z0])
nx=length(x0);
nz=length(z0);
F=F_fun(t0,[x0;z0]);
J=Jac_fun(t0,[x0;z0]);



dx=F(1:nx);
gx=J(nx+1:end,1:nx);
gz=J(nx+1:end,nx+1:end);

dz=-(gz)^(-1)*gx*dx;
end

function [value,isterminal,direction] = events_mother(diff_miu_hieq, ieq_act_s)
% Locate the time when height passes through zero in a 
% decreasing direction and stop integration.
n_ieq=length(diff_miu_hieq);
value = diff_miu_hieq;     % Detect height = 0
isterminal = 1*ones(n_ieq,1);   % Stop the integration
% direction = 0*ones(n_ieq,1);    % detect all crossings (negtive+positive direction)
direction=zeros(n_ieq,1);
for i=1:n_ieq
    if ieq_act_s(i)==1 %active
       direction(i)=-1; %detect from + to -
    else  %nonactive
       direction(i)=1; %detect from - to +
    end
end


end
%=============================================================================
%  OCDES: Solving Optimization-Constrained Dynamic Systems by optimality
%         conditions
%  
%  Copyright (c) 2020: Xiao Zhao, Forschungszentrum Juelich GmbH, Juelich,
%                      Germany. 
%
%  This code can only be used for academic purposes.                               
%  All rights reserved.
%=============================================================================
