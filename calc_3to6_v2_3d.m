nloop=0;
cpl=1;
moist=1;
opt=3;
for sf=0:nloop

    dt=5; %seconds
    ndays=60;
    
    n=ndays*24*60*60/dt;
    
    H=50; %surface layer thickness
    dz=1;
    vlev=100;
    vleva=50;
    rhow=1029; %kg/m^3
    cp=4184; %J/K/kg
    cpa=1005; %J/K/kg
    sb = 5.670e-8;
    d2k=273.15;
    Av=2500000;
    gamma_blk=0.01;
    cvapor_fac=640380;
    cvapor_exp=5107.4;
    cvapor_exp_ice=5897.8;
    saltsat=0.980;
    p0=101300; %reference pressure
    gravity_mks=9.81;
    humid_fac=0.606;
    tref=25; %reference temperature
    sref=25;
    S=sref; % salinity not implemented
    alpha=2e-4; %expansion coefficients for temperature
    beta=7.4e-4; %expansion coefficients for salt
    emissivity=0.97; % ocean emissivity
    kappa=1e-5; % m^2/kg air absorption coefficient
    Rgas=287.05;
    karman=0.4;
    
    ZA=H*(1:vleva);
    Z(2:vlev+1)=-cumsum(ones(1,vlev)*dz);
    kT0=ones(1,vlev)*1e-4; % SST vertical mixing
    kT=kT0; % SST vertical mixing
    kA=1e-4; % wind vertical mixing
    kq=1e-4; % moist vertical mixing
    kTH=1e-4; % air temperature vertical mixing
    
    R=0.62;
    g1=0.6;
    g2=20;
    
    f=R.*exp(Z./g1)+(1-R).*exp(Z./g2);
    f(Z<-200)=0;
    f=-diff(f,1);
    
    albedo=0.06;
    lat=30;
    lon=0;
    time=(1:n).*dt./3600;
    jstart=70;% julian day, borial summer
    jadd=fix((time-time(1))/24);
    jd=jstart+jadd;
    
    AOI=angle_of_incidence(lat,lon,jd,time); AOI(AOI<0)=0;
    SW=(1362*(1-albedo).*repmat(f',1,n).*repmat(AOI,vlev,1))';
    %mean(SW(:,1))
    %SW=(repmat((sin(2*pi/n*(1:n)*ndays+pi/2)+1)/2*1362*(1-albedo),vlev,1).*repmat(f',1,n))'.*cos(lat/180*pi);

    if sf>0
        if opt==1
            SST=circshift(SST2,60*sf*3);
        elseif opt==2
            frac=(1-(1:nloop)/nloop);
            SST=mean(SST2)+(SST2-mean(SST2))*(frac(sf));
        end
    end
    
    T=zeros(n+1,vlev);
    UO=zeros(n+1,vlev);
    SST=zeros(1,n+1);
    U=zeros(n+1,vleva);
    q=zeros(n+1,vleva);
    THETA=zeros(n+1,vleva);
    rhoa=zeros(n+1,vleva);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %U0=ones(n,1)*5;
    %THETA0=(sin(2*pi/n*(1:n)*ndays/5-pi/2)+21);
    TA=17.-0.007.*H.*(1:vleva);
    P=zeros(1,vleva);
    P(1)=101300;
    P(2:vleva)=P(1).*exp(-cumsum(H./(Rgas.*(TA(2:end)+d2k)./gravity_mks)));
    rhoa(1,:)=P./(Rgas.*(TA+d2k));
    THETA(1,:)=(TA+d2k).*(P(1)./P).^(2/7)-d2k;
    THETA(1,1:20)=THETA(1,1);
    THETA0=ones(n,1)*THETA(1,end)+1;%+gamma_blk*H+0.01;%+0.51;%+SW*kappa/cpa*dt;
    T0=ones(n+1,vlev).*THETA(1,end);%repmat((16+(vlev-1:-1:0)/vlev*8),n+1,1);
    T(:,:)=T0;
    T(:,1)=T0(:,1); T(1,:)=T0(1,1)+15*((vlev:-1:1)/vlev-1);
    SST(:)=T(:,1);
    ssq = saltsat*cvapor_fac*exp(-cvapor_exp./(TA(1)+d2k))./rhoa(1,1);
    RH=0.5;
    q0=RH*ssq;
    %Tv=(THETA(i-1,:)+d2k).*(1+humid_fac*q(i-1));
    q(1,:)=saltsat*cvapor_fac*exp(-cvapor_exp./(TA+d2k))./rhoa(1,:).*0.8;
    U(1,:)=(0.2/karman).*(log((ZA-H/2)/1e-4));
    %U(1,:)=U(1,end);
    U0=repmat((0.2/karman).*(log((ZA(end)+H/2)/1e-4)),1,n+1);
    bulk_in.atemp=squeeze(THETA(1,1)+d2k); bulk_in.aqh=squeeze(q(1,1));
    bulk_in.speed=squeeze(abs(U(1,1))); bulk_in.sst=squeeze(SST(1));
    tmp0=exf_bulk_largeyeager04_stress(bulk_in,H/2,H/2,H/2,H/2,rhoa(1,1));
    huol0=tmp0.xi;
    %U0=(sin(2*pi/n*(1:n)*ndays/5+2*pi/n*3600/dt*sf)*8+U(1,end));
    %UO(:)=0.;
    %U0(3600*24*5/dt:end)=U(1,end); one pertubation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    LH=zeros(1,n); SH=zeros(1,n); TAU=zeros(1,n); TAUO=zeros(1,n); 
    LW=zeros(1,n); E=zeros(1,n); ssq=zeros(1,n);
    xi=zeros(1,n); ce=zeros(1,n); ch=zeros(1,n);
    Km=zeros(n,vleva); Kh=zeros(n,vleva); 
    %Km2=zeros(n,vleva); Kh2=zeros(n,vleva);
    rd=zeros(1,n); rh=zeros(1,n); re=zeros(1,n);
    Qnet=zeros(1,n); Kto=zeros(n,vlev);Kmo=zeros(n,vlev);
    ustar=zeros(1,n); qstar=zeros(1,n); tstar=zeros(1,n);
    psimh=zeros(1,n); psixh=zeros(1,n); Ri=zeros(n,vleva);
    
    for i=2:n+1
        %i
        if rem(i,1440*60/dt)==0
            display(['Day ' num2str(i/(1440*60/dt))])
        end
        
        bulk_in.atemp=squeeze(THETA(i-1,1)+d2k);
        bulk_in.aqh=squeeze(q(i-1,1));
        bulk_in.speed=squeeze(abs(U(i-1,1)-UO(i-1,1)));
        bulk_in.sst=squeeze(SST(i-1));
        tmp1=exf_bulk_largeyeager04_stress(bulk_in,H/2,H/2,H/2,H/2,rhoa(i-1,1));
        LH(i-1)=-tmp1.hl;
        SH(i-1)=-tmp1.hs;
        %tau(i-1)=tmp1.tau;
        LW(i-1)=emissivity*sb*(SST(i-1)+d2k).^4;
        E(i-1)=-tmp1.hl/Av;
        ssq(i-1)=tmp1.ssq;
        xi(i-1)=tmp1.xi;
        ce(i-1)=tmp1.ce;
        ch(i-1)=tmp1.ch;
        rd(i-1)=tmp1.rd;
        rh(i-1)=tmp1.rh;
        re(i-1)=tmp1.re;
        ustar(i-1)=tmp1.ustar;
        qstar(i-1)=tmp1.qstar;
        tstar(i-1)=tmp1.tstar;
        psimh(i-1)=tmp1.psimh;
        psixh(i-1)=tmp1.psixh;
        TAU(i-1)=ustar(i-1).^2./rhoa(i-1,1);
        %TAUO(i-1)=ustar(i-1).^2./rhoa(i-1,1);
        
        Qnet(i-1)=SW(i-1,1)-LH(i-1)-SH(i-1)-LW(i-1);
        %p=p0.*exp(-cumsum(H*gravity_mks./Rgas./([THETA(i-1,:) THETA0(i-1)]+d2k)));
        %rhoa(i-1,:)=(p./(THETA0(i-1)+d2k)./Rgas);
        %tfac=((p0./p).^(2/7));
        
        Tv0=(THETA0(i-1)+d2k)*(1+humid_fac*q0);
        Tv=(THETA(i-1,:)+d2k).*(1+humid_fac*q(i-1));
        %Ri0=(gravity_mks./Tv).*(diff([Tv Tv0],1))*H./(diff([U(i-1,:) U0(i-1)],1)).^2;
        Ri0=(gravity_mks./Tv).*(Tv-Tv0).*ZA./(U(i-1,:)).^2;
        Ri(i-1,:)=min(abs(Ri0),10).*sign(Ri0);
        %[Km2, Kh2]=louis(Ri(i-1,:),abs(ZA),diff([U(i-1,:) U0(i-1)],1),H);
        %Km(i-1,:)=Km2;
        %Kh(i-1,:)=Kh2;
        
        [Km2, Kh2]=...
            holtslag(Tv,q(i-1,1),U(i-1,:),U(i-1,:).*0,LH(i-1),SH(i-1),ustar(i-1),[0 ZA]);
        %Km(i-1,Km2~=0)=Km2(Km2~=0);
        %Kh(i-1,Kh2~=0)=Km2(Kh2~=0);
        
        rho=rhow*(1-alpha*(T(i-1,:)-tref)+beta*(S-sref));
        B=gravity_mks*(rhow-rho)./rhow;
        Rig(1:vlev-1)=-diff(B,1)./dz./(U(i-1).^2./(vlev*dz).^2);
        %[Kmo(i-1,1:vlev-1), Kto(i-1,1:vlev-1)]=Pacanowski(Rig);
        Uo=(vlev:-1:1)/vlev/10; Vo=zeros(size(Uo)); FWflux=zeros(size(U));
        [Kmo(i-1,:), Kto(i-1,:), ~]=...
            large(T(i-1,:),repmat(S,1,vlev),UO(i,:),Vo,Qnet(i-1),0,ustar(i-1),Z,lat);
        
        kT=Kto(i-1,:);
        kM=Kmo(i-1,:);
        
        klim=inf;%repmat(0.5.*dt./H.^2,1,vleva);
        kA=min(Km2,klim);
        kq=min(Kh2,klim);
        kTH=min(Kh2,klim);
        %[i max(kA) max(kq)]
        
        TAUS=sign(U(i-1,1)-UO(i-1,1));
        
        tmp1=U(i-1,1)-TAU(i-1)*TAUS/(H*rhoa(i-1,1))*dt...
            -((kA(1)+kA(2))/2).*(U(i-1,1)-U(i-1,2))/H^2*dt;
        tmp2=U(i-1,end)-kA(end)*(U(i-1,end)-U0(i-1))/H^2*dt;
        U(i,2:end-1)=U(i-1,2:end-1)...
            +(((kA(1:end-2)+kA(2:end-1))/2).*(U(i-1,1:end-2)-U(i-1,2:end-1)))/H^2*dt...
            -(((kA(2:end-1)+kA(3:end  ))/2).*(U(i-1,2:end-1)-U(i-1,3:end  )))/H^2*dt;
        U(i,1)=tmp1; U(i,end)=tmp2;
        
        tmp1=q(i-1,1)-E(i-1)/(H*rho(1))*dt...
            -((kq(1)+kq(2))/2)*(q(i-1,1)-q(i-1,2))/H^2*dt;
        tmp2=q(i-1,end)-kq(end)*(q(i-1,end)-q0)/H^2*dt;
        q(i,2:end-1)=q(i-1,2:end-1)...
            +(((kq(1:end-2)+kq(2:end-1))/2).*(q(i-1,1:end-2)-q(i-1,2:end-1)))/H^2*dt...
            -(((kq(2:end-1)+kq(3:end  ))/2).*(q(i-1,2:end-1)-q(i-1,3:end  )))/H^2*dt;
        q(i,1)=tmp1; q(i,end)=tmp2;
        
        tmp1=THETA(i-1,1)+(SH(i-1))/cpa/(H*rhoa(i-1,1))*dt...
            -((kTH(1)+kTH(2))/2)*(THETA(i-1,1)-THETA(i-1,2))/H^2*dt;
        tmp2=THETA(i-1,end)-kTH(end)*(THETA(i-1,end)-THETA0(i-1))/H^2*dt;
        THETA(i,2:end-1)=THETA(i-1,2:end-1)...
            +(((kTH(1:end-2)+kTH(2:end-1))/2).*(THETA(i-1,1:end-2)-THETA(i-1,2:end-1)))/H^2*dt...
            -(((kTH(2:end-1)+kTH(3:end  ))/2).*(THETA(i-1,2:end-1)-THETA(i-1,3:end  )))/H^2*dt;
        THETA(i,1)=tmp1; THETA(i,end)=tmp2;
        
        TA=(THETA(i,:)+d2k).*(P./P(1)).^(2/7)-d2k;
        rhoa(i,:)=P./(Rgas.*(TA+d2k));
        
        if moist
            qsat=saltsat*cvapor_fac*exp(-cvapor_exp./(TA+d2k))./rhoa(i-1,:);
            dq=(q(i,:)-qsat).*((q(i,:)-qsat)>0);
            dq=min(dq,q(i,:));
            dTHETA=Av.*dq./cpa; %[J/kg][j-1 kg K]
            %[i max(dTHETA) max(dq)]
            ind=q(i,:)>qsat;
            q(i,ind)=q(i,ind)-dq(ind);
            THETA(i,:)=THETA(i,: )+dTHETA;
        end        
        
        if cpl
            SST(i)=SST(i-1)+(SW(i-1,1)-LH(i-1)-SH(i-1)-LW(i-1))/cp/(dz*rhow)*dt-...
                ((kT(1)+kT(2))/2)*(SST(i-1)-T(i-1,2))/dz^2*dt;
            tmp=T(i-1,end)+(kT(    end-1).*(T(i-1,end-1)-T(i-1,end)))/dz^2*dt+...
                SW(i-1,end)/cp/(dz*rhow)*dt;
            T(i,2:end-1)=T(i-1,2:end-1)...
                +(((kT(1:end-2)+kT(2:end-1))/2).*(T(i-1,1:end-2)-T(i-1,2:end-1)))/dz^2*dt...
                -(((kT(2:end-1)+kT(3:end  ))/2).*(T(i-1,2:end-1)-T(i-1,3:end  )))/dz^2*dt...
                +SW(i-1,2:end-1)/cp/(dz*rhow)*dt;
            T(i,1)=SST(i);
            T(i,end)=tmp;
            
            utemp1=UO(i-1,1)...
                +(TAU(i-1))*TAUS/(rho(1))/dz^2*dt...
                -((kM(1)+kM(2))/2)*(UO(i-1)-UO(i-1,2))/dz^2*dt;
            utemp2=UO(i-1,end)+(kM(    end-1).*(UO(i-1,end-1)-UO(i-1,end)))/dz^2*dt;
            UO(i,2:end-1)=UO(i-1,2:end-1)...
                +(((kM(1:end-2)+kM(2:end-1))/2).*(UO(i-1,1:end-2)-UO(i-1,2:end-1)))/dz^2*dt...
                -(((kM(2:end-1)+kM(3:end  ))/2).*(UO(i-1,2:end-1)-UO(i-1,3:end  )))/dz^2*dt;
            UO(i,1)=utemp1;
            UO(i,end)=utemp2;
        end
        if rem(i,3600*6/dt)==0; nsub=6;... %3600*6/dt
                          subplot(nsub,1,1); plot(-Z(2:end),T(i,:)); ylim([19 21]);    title('T');...
                          subplot(nsub,1,2); plot(-Z(2:end),UO(i,:));    title('UO');...
                          subplot(nsub,1,3); plot(ZA,U(i,:));     title('U');...
                          subplot(nsub,1,4); plot(ZA,Tv-273.15); title('Tv');...%ylim([19 21]);...
                          subplot(nsub,1,5); plot(ZA,q(i,:));     title('q');...
                          subplot(nsub,1,6); plot(ZA,Kh2);     title('Kh');...
                          suptitle(num2str(i*dt/3600/24)); drawnow;  
        end
                           
    end
    
    if sf==0
        SST2=SST;
    end
    
    MLH(sf+1)=mean(LH);
    MSH(sf+1)=mean(SH);
    MSST(sf+1)=mean(SST);
    
    figure(1)
    
    subplot(711);
    plot(movmean(U0,24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    %ylim([1 2.5])
    title('U0')
    
    subplot(712);
    plot(movmean(U(:,end),24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    %ylim([1 2.5])
    title('Utop')
    
    subplot(713);
    plot(movmean(U(:,end)'-U0,24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    %ylim([1 2.5])
    title('Utop-U0')
    
    subplot(714);
    plot(movmean(U(:,1),24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    %ylim([1 2.5])
    title('U')
    
    subplot(715);
    plot(movmean(THETA(:,2),24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    title('THETA')
    
    subplot(716);
    plot(movmean(SST,24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    title('SST')
    
    subplot(717);
    plot(movmean(T(:,10:10:end),24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    title('T')
    
    if opt==1
        tt=['SST shift hours = ' num2str(sf*3)];
    elseif opt==2
        if sf==0
            tt=['SST amplitude fraction = 1'];
        else
            tt=['SST amplitude fraction = ' num2str(frac(sf))];
        end
    elseif opt==3
        tt=['DZ = ' num2str(dz)];
    end
   suptitle(tt)
    
    
    figure(2)
    
    subplot(711);
    plot(xi)
    title('xi')
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    
    subplot(712);
    plot(ustar)
    title('ustar')
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    
    subplot(713);
    plot(psimh)
    title('psimh')
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on

    subplot(714);
    plot(psixh)
    title('psixh')
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    
    subplot(715);
    plot(rd.^2)
    title('CD')
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    
    subplot(716);
    plot(rh.*rd)
    title('CH')
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    
    subplot(717);
    plot(re.*rd)
    title('CE')
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    
    figure(3)
    subplot(511);
    plot(movmean(SW(:,1)/cp/(dz*rhow)*dt,24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    %ylim([1 2.5])
    title('SW')
    
    subplot(512);
    plot(movmean(LW/cp/(dz*rhow)*dt,24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    title('LW')
    
    subplot(513);
    plot(movmean(SH/cp/(dz*rhow)*dt,24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    title('SH')
    
    subplot(514);
    plot(movmean(LH/cp/(dz*rhow)*dt,24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    title('LH')
    
    subplot(515);
    plot(movmean((SW(:,1)'-LW-SH-LH)/cp/(dz*rhow)*dt,24*60*60/dt,'Endpoints','discard'))
    set(gca,'XTick',0:24*60*60/dt:n)
    set(gca,'XTickLabels',0:ndays)
    grid on
    title('QNET')
  
    tt=['dt = ' num2str(sf) '  hours'];
    suptitle(tt)
    

    dSST=diff(movmean(mean(reshape(SST(2:end),[60*60/dt ndays*24]),1),24,'Endpoints','discard'),1);
    dU=diff(movmean(mean(reshape(U(2:end,1),[60*60/dt ndays*24]),1),24,'Endpoints','discard'),1);

    [r, lag] = xcorr(dU,dSST,120,'coeff');
    figure(4)
    h=plot(lag,r);
    xlim([-120 120])
    ylim([-1 1])
    title('corr')
    set(gca,'XTick',-120:24:120)
    set(gca,'XTickLabel',-120/24:24/24:120/24)
    
end
%{
figure(2)
subplot(311);
plot(MLH)
subplot(312);
plot(MSH)
subplot(313);
plot(MSH+MLH)
%}



