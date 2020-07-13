function [tout,yout]=sOCDE_EF_main(f, obj, heq, hieq,x,x0,v,v0,tstart, tfinal,opt_opt,opt_ode)
%This function solves OCDE by EF-based approach

%Input:
%f: state equations of the dynamic part (sym)
%obj: objective function (sym)
%heq: left hand side of equality constraints (sym)
%hieq: left hand side of inequality (>=) constraints (sym)
%x:  state varaibles(sym)
%x0: initial point of x (real column vector). x0 must be given and fixed.
%v:  inner optimization variables (sym)
%v0: initial guess of inner optimization variables (real column vector). v0 can be empty.
%opt_opt: options for optimization
%opt_sol: options for integrator ode45s

%Outputs:
%tout: integration steps of t 
%yout: integration results of x


matlabFunction(obj,'File','obj_gen','Vars',{x,v});
obj_fun=@(x,v)obj_gen(x,v);

if ~isempty(hieq)
matlabFunction(hieq,'File','hieq_gen','Vars',{x,v});
hieq_fun=@(x,v)hieq_gen(x,v);
nonlcon_fun_cell{1}=@(x,v)-hieq_fun(x,v);
else
hieq_fun=[];
nonlcon_fun_cell{1}=[];
end  

if ~isempty(heq)
matlabFunction(heq,'File','heq_gen','Vars',{x,v});
heq_fun=@(x,v)heq_gen(x,v);
nonlcon_fun_cell{2}=@(x,v)heq_fun(x,v);
else
heq_fun=[];
nonlcon_fun_cell{2}=[];
end
  
%do this here: nonlcon_fun_cell={@(x)-hieq_fun(x),@(x)heq_fun(x)};
nonlcon_fun=@(x,v)MO_implicitFun(nonlcon_fun_cell,x,v);

global v0_ %initialization point of P(x)
v0_=v0;

v_opt_fun=@(x)solveParaOP(x,obj_fun,nonlcon_fun,opt_opt);

matlabFunction(f,'File','f_gen','Vars',{x,v});

f_fun=@(t,x)f_gen(x,v_opt_fun(x));

[tout,yout] =ode45(@(t,x)f_fun(t,x),[tstart tfinal],x0,opt_ode);

end


function [v_sol,fval,exitflag,output,lambda,grad,hessian]=solveParaOP(x,obj_fun,nonlcon_fun, opt_opt)
%sovel parametrized P(x)=min_v obj(x,v)
%Input
%obj_fun(x,v)
%nonlcon_fun=@(x,v){cieq(x,v),ceq(x,v)} with cieq<=0, ceq==0

Prob.objective=@(v)obj_fun(x,v); 
Prob.Aineq=[];
Prob.bineq=[]; 
Prob.Aeq=[]; 
Prob.beq=[];
Prob.lb=[];
Prob.ub=[];
Prob.nonlcon=@(v)nonlcon_fun(x,v); 
Prob.solver='fmincon';

global v0_
if ~isempty(v0_)
Prob.x0=v0_;
else
Prob.x0=zeros(length(v),1)    ;
end

if isempty(opt_opt)
Prob.options = optimoptions(@fmincon,'Algorithm','sqp');
else
Prob.options=opt_opt;
end

[v_sol,fval,exitflag,output,lambda,grad,hessian]=fmincon(Prob);

if exitflag~=1
    y=[];
    msg = 'Error occurred: Unable to solve the inner NLP at the start time. Try to provide good initial guesses.';
    %error(msg);    
    fprintf('%s',msg);
else
    v0_=v_sol;
end

end


