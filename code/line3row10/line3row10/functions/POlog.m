% log-sum_th - threshold function for Log-sum regularization 
function [STq,pt,ss]=POlog(q,t,delta)
tao=2*t.^(0.5)-delta;
if max(tao<0)
%     temp=1;
    tao(tao<0)=0;
end
pt=((abs(q)-tao)>0);%pt(pt>0)=1;pt(pt<=0)=0;
u=0.5*(abs(q)-delta+(((abs(q)+delta).^2-4*t).*pt).^(0.5));
STq=u.*sign(q).*pt;
if ~isreal(STq)
    temp=1;
end
ss=abs(STq);