function f = f_gen(in1,in2)
%F_GEN
%    F = F_GEN(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    14-Jul-2020 00:16:56

v1 = in2(1,:);
v2 = in2(2,:);
x1 = in1(1,:);
f = [v1.*v2+1.0./2.0;x1];