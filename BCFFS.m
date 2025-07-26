%%
% function to compute f(0), f(1), fdash(0) and fdash(1)
function [f0,f1,df0,df1] = BCFFS(a0,abmatrix)
nffs=size(abmatrix,1);
f0=a0/2;
f1=a0/2;
df0=0;
df1=0;
    for i = 1:nffs
        acoeff=abmatrix(i,1);
        bcoeff=abmatrix(i,2);
        f0=f0 + acoeff ; % tau = 0
        f1 = f1 + acoeff*cos(i*pi) + bcoeff*sin(i*pi); % tau = 1
        df0 = df0 + bcoeff*i*pi;
        df1 = df1 - acoeff*i*pi*sin(i*pi)+ bcoeff*i*pi*cos(i*pi);
    end
end