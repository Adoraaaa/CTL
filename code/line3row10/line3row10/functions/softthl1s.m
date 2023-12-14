% l1_softth - soft threshold function for L1 regularization
function [STq,pt]=softthl1s(q,t)
pt=((abs(q)-t)>0);
STq=(q-t.*q./abs(q)).*pt;
%ss=abs(STq);