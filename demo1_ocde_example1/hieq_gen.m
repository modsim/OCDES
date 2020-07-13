function hieq = hieq_gen(in1)
%HIEQ_GEN
%    HIEQ = HIEQ_GEN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    14-Jul-2020 00:17:02

v1 = in1(3,:);
v2 = in1(4,:);
x1 = in1(1,:);
x2 = in1(2,:);
hieq = [v2-exp(-v1)+sin(x1).*(1.7e1./1.0e1);-v2+cos(x2).*(1.0./5.0)+2.0];