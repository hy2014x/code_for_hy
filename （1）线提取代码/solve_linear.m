function [t,num]=solve_linear(a,b)
if(abs(a)<10^-4) %�������a==0
%if(a~=0)
    num=0;
    t=0;
    return;
else
    num=1;
    t=-b/a;
    return;
end
end