function y=heaviside(x)
y=zeros(size(x));
y(x>=0)=1;
