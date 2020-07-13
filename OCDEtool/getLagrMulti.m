function y=getLagrMulti(obj, heq, hieq, x, x0,opt_opt)
%This function returns solution of a NLP and its Lagrange Multipliers
%Input:
%min obj(v), heq(v)==0 hieq(v)>=0
%
%
%Output:
%y.y_opt: solution of NLP
%y.lam.hieq: Lagrange Multipliers of inequalities
%y.lam.heq: Lagrange Multipliers of equalities

matlabFunction(obj,'File','obj_gen','Vars',{x});
obj_fun=@(x)obj_gen(x);
if ~isempty(hieq)
matlabFunction(hieq,'File','hieq_gen','Vars',{x});
hieq_fun=@(x)hieq_gen(x);
nonlcon_fun_cell{1}=@(x)-hieq_fun(x);
else
hieq_fun=[];
nonlcon_fun_cell{1}=[];
end  

if ~isempty(heq)
matlabFunction(heq,'File','heq_gen','Vars',{x});
heq_fun=@(x)heq_gen(x);
nonlcon_fun_cell{2}=@(x)heq_fun(x);
else
heq_fun=[];
nonlcon_fun_cell{2}=[];
end
  
%nonlcon_fun_cell={@(x)-hieq_fun(x),@(x)heq_fun(x)};
nonlcon_fun=@(x)MO_implicitFun(nonlcon_fun_cell,x);

Prob.objective=obj_fun; 
Prob.Aineq=[];
Prob.bineq=[]; 
Prob.Aeq=[]; 
Prob.beq=[];
Prob.lb=[];
Prob.ub=[];
if ~isempty(x0)
Prob.x0=x0;
else
Prob.x0=zeros(length(x),1)    ;
end
Prob.nonlcon=nonlcon_fun; 

if isempty(opt_opt)
Prob.options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500);
else
Prob.options=opt_opt;
end

Prob.solver='fmincon';
[x_sol,fval,exitflag,output,lambda,grad,hessian]=fmincon(Prob);

if exitflag~=1
    y=[];
    msg = 'Error occurred: Unable to solve the inner NLP at the start time. Try to provide good initial guesses.';
    error(msg);    
else
    y.y_opt=x_sol;
end

y.lam.heq=lambda.eqnonlin ;
y.lam.hieq=lambda.ineqnonlin ;

end

%=============================================================================
%  OCDES – Solving Optimization-Constrained Dynamic Systems by optimality
%         conditions
%  
%  Copyright (c) 2020: Xiao Zhao, Forschungszentrum Juelich GmbH, Juelich,
%                      Germany. 
%
%  This code can only be used for academic purposes.                               
%  All rights reserved.
%=============================================================================