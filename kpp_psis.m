function y=kpp_psis(xi)
y(xi>=0)=1+5.*xi(xi>=0);
y(xi<0 & xi>=-1)=(1-16.*xi(xi<0 & xi>=-1)).^(-1/2);
y(xi<-1)=(-28.86-98.96.*xi(xi<-1)).^(-1/3);
end