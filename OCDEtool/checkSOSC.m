function [isSOSC,x_sol]=checkSOSC(H,dheq,dhieq,ieq_act,miu_ieq)

x_sol=[];

obj_fun=@(s)s'*H*s;
Aeq=dheq;
beq=zeros(size(dheq,1),1);
A=[];
b=[];

for i=1:numel(ieq_act)
    if ieq_act(i)==0
       continue; 
    end
    
    if abs(miu_ieq(i))<1e-3
       fprintf('Checking SOSC: SCC is voilated.\n');
       A=[A;
          -dhieq(i,:)];
       b=[b;
          0];
    else
       Aeq=[Aeq;
            dhieq(i,:)];
       beq=[beq;
           0]; 
    end
        
        
end


lb=[];
ub=[];
options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500, 'Display','off');
exitflag=0;
fval=Inf;
MaxTry=10;
for i=1:MaxTry
x0=rand(size(dheq,2),1);
[x_i,fval_i,exitflag_i] = fmincon(obj_fun,x0,A,b,Aeq,beq,lb,ub,...
    @(s)nonlcon_SOSC(s),options);

if exitflag_i==-2 %no feasible point
   fval=fval_i;
   exitflag=-2;   
end

if exitflag_i==1 && fval_i<fval
    x_sol=x_i;
    fval=fval_i;
    exitflag=exitflag_i;
end
end



if exitflag==-2
    isSOSC=1;
    fprintf('Check SOSC: No feasible point for Qp found.\n ');  
    return;
end

 
if exitflag==1
   if fval<1e-6
      isSOSC=0;
   else
      isSOSC=1;
   end
   return;
end

isSOSC=99;
fprintf('Check SOSC fails. \n');


end



function [c,ceq]=nonlcon_SOSC(s)
ceq=s'*s-1;
c=[];

end