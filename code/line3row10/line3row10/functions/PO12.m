function [STq,pt,ss]=PO12(q,t)
pt=((abs(q)-(t).^(2/3)*3/2)>0);%pt(pt>0)=1;pt(pt<=0)=0;
temp=1+cos(2/3*acos(-3^(1.5)/4*t.*(abs(q)).^(-1.5)));
STq=(2/3)*q.*temp.*pt;
ss=abs(STq);