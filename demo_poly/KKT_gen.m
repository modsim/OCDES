function KKT = KKT_gen(in1,in2)
%KKT_GEN
%    KKT = KKT_GEN(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    20-Jul-2019 19:29:19

ieq_act_1 = in2(1,:);
ieq_act_2 = in2(2,:);
miu_ieq_1 = in1(4,:);
miu_ieq_2 = in1(5,:);
p4 = in1(2,:);
z = in1(3,:);
KKT = sparse([1,2,3],[1,1,1],[-miu_ieq_1+miu_ieq_2+conj(p4)-conj(z).*5.0+conj(z).^2.*(3.0./2.0)+conj(z).^3.*4.0,miu_ieq_1.*(ieq_act_1-1.0)+ieq_act_1.*(z+3.0./2.0),miu_ieq_2.*(ieq_act_2-1.0)-ieq_act_2.*(z-3.0./2.0)],3,1);
