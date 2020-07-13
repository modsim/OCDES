function KKT = KKT_gen(in1,in2)
%KKT_GEN
%    KKT = KKT_GEN(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    20-Oct-2019 11:17:27

cs1 = in1(2,:);
cs2 = in1(3,:);
ieq_act_1 = in2(1,:);
ieq_act_2 = in2(2,:);
ieq_act_3 = in2(3,:);
miu_eq_1 = in1(16,:);
miu_eq_2 = in1(17,:);
miu_eq_3 = in1(18,:);
miu_eq_4 = in1(19,:);
miu_eq_5 = in1(20,:);
miu_eq_6 = in1(21,:);
miu_eq_7 = in1(22,:);
miu_ieq_1 = in1(13,:);
miu_ieq_2 = in1(14,:);
miu_ieq_3 = in1(15,:);
v1 = in1(7,:);
v2 = in1(8,:);
v3 = in1(9,:);
v4 = in1(10,:);
v5 = in1(11,:);
v6 = in1(12,:);
vupt1 = in1(5,:);
vupt2 = in1(6,:);
KKT = sparse([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[miu_eq_1+miu_eq_6+(vupt1.*(1.0./2.0)+conj(vupt1).*(1.0./2.0))./conj(sqrt(v1.*conj(v1)+v2.*conj(v2)+v3.*conj(v3)+v4.*conj(v4)+v5.*conj(v5)+v6.*conj(v6)+vupt1.*conj(vupt1)+vupt2.*conj(vupt2))),miu_eq_2+miu_eq_7+(vupt2.*(1.0./2.0)+conj(vupt2).*(1.0./2.0))./conj(sqrt(v1.*conj(v1)+v2.*conj(v2)+v3.*conj(v3)+v4.*conj(v4)+v5.*conj(v5)+v6.*conj(v6)+vupt1.*conj(vupt1)+vupt2.*conj(vupt2))),-miu_eq_1+miu_eq_2-miu_eq_4-miu_ieq_1+(v1.*(1.0./2.0)+conj(v1).*(1.0./2.0))./conj(sqrt(v1.*conj(v1)+v2.*conj(v2)+v3.*conj(v3)+v4.*conj(v4)+v5.*conj(v5)+v6.*conj(v6)+vupt1.*conj(vupt1)+vupt2.*conj(vupt2))),-miu_eq_2+miu_eq_3+miu_eq_5+(v2.*(1.0./2.0)+conj(v2).*(1.0./2.0))./conj(sqrt(v1.*conj(v1)+v2.*conj(v2)+v3.*conj(v3)+v4.*conj(v4)+v5.*conj(v5)+v6.*conj(v6)+vupt1.*conj(vupt1)+vupt2.*conj(vupt2))),-miu_eq_3+miu_eq_4+miu_eq_5+(v3.*(1.0./2.0)+conj(v3).*(1.0./2.0))./conj(sqrt(v1.*conj(v1)+v2.*conj(v2)+v3.*conj(v3)+v4.*conj(v4)+v5.*conj(v5)+v6.*conj(v6)+vupt1.*conj(vupt1)+vupt2.*conj(vupt2))),-miu_eq_1+miu_eq_4-miu_ieq_2+(v4.*(1.0./2.0)+conj(v4).*(1.0./2.0))./conj(sqrt(v1.*conj(v1)+v2.*conj(v2)+v3.*conj(v3)+v4.*conj(v4)+v5.*conj(v5)+v6.*conj(v6)+vupt1.*conj(vupt1)+vupt2.*conj(vupt2))),-miu_eq_4-miu_ieq_3+(v5.*(1.0./2.0)+conj(v5).*(1.0./2.0))./conj(sqrt(v1.*conj(v1)+v2.*conj(v2)+v3.*conj(v3)+v4.*conj(v4)+v5.*conj(v5)+v6.*conj(v6)+vupt1.*conj(vupt1)+vupt2.*conj(vupt2))),-miu_eq_5+(v6.*(1.0./2.0)+conj(v6).*(1.0./2.0))./conj(sqrt(v1.*conj(v1)+v2.*conj(v2)+v3.*conj(v3)+v4.*conj(v4)+v5.*conj(v5)+v6.*conj(v6)+vupt1.*conj(vupt1)+vupt2.*conj(vupt2))),-v1-v4+vupt1,v1-v2+vupt2,v2-v3,-v1+v3+v4-v5,v2+v3-v6,vupt1-(cs1.*(1.9e1./5.0))./(cs1+1.0),vupt2-(cs2.*(1.1e1./1.0e1))./(cs2+1.0),ieq_act_1.*v1+miu_ieq_1.*(ieq_act_1-1.0),ieq_act_2.*v4+miu_ieq_2.*(ieq_act_2-1.0),ieq_act_3.*v5+miu_ieq_3.*(ieq_act_3-1.0)],18,1);
