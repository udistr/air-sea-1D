    function y=kpp_psim(xi)
        y(xi>=0)=1+5.*xi(xi>=0);
        y(xi<0 & xi>=-0.2)=(1-16.*xi(xi<0 & xi>-0.2)).^(-1/4);
        y(xi<-0.2)=(1.26-8.38.*xi(xi<-0.2)).^(-1/3);
    end