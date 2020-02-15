function [mystruct_out]=exf_bulk_largeyeager04_stress(mystruct_in,varargin)

%varargin	[hu ht hq]
%		[hu ht hq],dragCoeff

tmp1=fieldnames(mystruct_in);
if length(tmp1)>6; doNetFluxes=1; end

% formulae in short:
%   wind stress = (ust,vst) = rhoA * Cd * Ws * (del.u,del.v)
%   Sensib Heat flux = fsha = rhoA * Ch * Ws * del.T * CpAir
%   Latent Heat flux = flha = rhoA * Ce * Ws * del.Q * Lvap
%                    = -Evap * Lvap
% with Ws = wind speed = sqrt(del.u^2 +del.v^2) ;
%      del.T = Tair - Tsurf ; del.Q = Qair - Qsurf ;
%      Cd,Ch,Ce = drag coefficient, Stanton number and Dalton number
%            respectively [no-units], function of height & stability


%parameters from exf_readparms.F
if nargin==1
    hu=10;
    ht=2; %differ from model version because all CORE data are given at 10m
    hq=2;
else
    hu=varargin{1};
    ht=varargin{2};
    hq=varargin{3};
    zref=varargin{4};
    atmrho=varargin{5};
end

umin=0.1;
karman=0.4;
gravity_mks=9.81;
cen2kel=273.150;
humid_fac=0.606;
cvapor_fac=640380;
cvapor_fac_ice=11637800;
cvapor_exp=5107.4;
cvapor_exp_ice=5897.8;
saltsat=0.980;
%atmrho=1.2;
gamma_blk=0.01;
cdrag_1=0.0027000;
cdrag_2=0.0001420;
cdrag_3=0.0000764;
cstanton_1 = 0.0327;
cstanton_2 = 0.0180;
cdalton = 0.0346;
niter_bulk=2;
psim_fac=5;
atmcp=1005;
flamb=2500000;
rhoConstFresh=999.8;
stefanBoltzmann = 5.670e-8;
ocean_emissivity=5.50e-8 / 5.670e-8;
albedo=0.1;

%-- Set surface parameters :
zwln = log(hu./zref);
ztln = log(ht./zref);
czol = hu*karman*gravity_mks;

%inputs:
atemp_in=mystruct_in.atemp;
aqh_in=mystruct_in.aqh;
sst_in=mystruct_in.sst;
speed_in=mystruct_in.speed;
wspeed=max(speed_in,umin);

%-   Surface Temp.
Tsf = sst_in + cen2kel;

%--- Compute turbulent surface fluxes
%-   Pot. Temp and saturated specific humidity

t0 = atemp_in.*(1 + humid_fac*aqh_in);
tmpbulk = cvapor_fac*exp(-cvapor_exp./Tsf);
ssq = saltsat*tmpbulk/atmrho;
deltap = atemp_in + gamma_blk*ht - Tsf;
delq   = aqh_in - ssq;

%--  initial guess for exchange coefficients:
%    take U_N = del.U ; stability from del.Theta ;
stable = 0.5 + 0.5*sign(deltap);
if nargin<7
    tmpbulk = cdrag_1./wspeed + cdrag_2 + cdrag_3.*wspeed; % CD
else
    %use prescribed drag coefficient
    tmpbulk=varargin{6};
end
rdn = sqrt(tmpbulk);
rhn = (1-stable)*cstanton_1 + stable*cstanton_2;
ren = cdalton;
%--  calculate turbulent scales
ustar=rdn.*wspeed;
tstar=rhn.*deltap;
qstar=ren*delq;

%--- iterate with psi-functions to find transfer coefficients
for iter=1:niter_bulk
    
    huol = ( tstar./t0 + qstar./(1./humid_fac+aqh_in) ).*czol./(ustar.*ustar);
    
    %cph The following is different in Large&Pond1981 code:
    %cph huol = max(huol,zolmin) with zolmin = -100
    tmpbulk = min(abs(huol),10);
    huol   = tmpbulk.*sign(huol);
    htol   = huol.*ht./hu;
    hqol   = huol.*hq./hu;
    stable = 0.5 + 0.5*sign(huol);
    
    %                 Evaluate all stability functions assuming hq = ht.
    % The following is different in Large&Pond1981 code:
    % xsq = max(sqrt(abs(1 - 16.*huol)),1)
    xsq    = sqrt( abs(1 - huol*16) );
    x      = sqrt(xsq);
    psimh = -psim_fac*huol.*stable + (1-stable).* ...
        ( log( (1 + 2*x + xsq).*(1+xsq)*.125 ) - 2*atan(x) + 0.5*pi );
    
    % The following is different in Large&Pond1981 code:
    % xsq = max(sqrt(abs(1 - 16.*htol)),1)
    xsq   = sqrt( abs(1 - htol*16) );
    psixh = -psim_fac*htol.*stable + (1-stable).*( 2*log(0.5*(1+xsq)) );
    
    %-   Shift wind speed using old coefficient
    usn = wspeed./(1 + rdn/karman.*(zwln-psimh) );
    usm = max(usn, umin);
    
    %-   Update the 10m, neutral stability transfer coefficients
    if nargin<7
        tmpbulk = cdrag_1./usm + cdrag_2 + cdrag_3.*usm;
    else
        %use prescribed drag coefficient
        tmpbulk=varargin{6};
    end
    rdn = sqrt(tmpbulk);
    rhn = (1-stable)*cstanton_1 + stable*cstanton_2;
    ren = cdalton;
    
    %-   Shift all coefficients to the measurement height and stability.
    rd = rdn./(1 + rdn.*(zwln-psimh)/karman);
    rh = rhn./(1 + rhn.*(ztln-psixh)/karman);
    re = ren./(1 + ren.*(ztln-psixh)/karman);
    
    %--  Update ustar, tstar, qstar using updated, shifted coefficients.
    ustar = rd.*wspeed;
    qstar = re.*delq;
    tstar = rh.*deltap;
    
    % end of iteration loop
end%for iter=1:niter_bulk

%-   Coeff:
tau   = atmrho*rd.*wspeed;

%-   Turbulent Fluxes
hs = atmcp*tau.*tstar;
hl = flamb*tau.*qstar;
%   change sign and convert from kg/m^2/s to m/s via rhoConstFresh
evap = -tau.*qstar/rhoConstFresh;


%mystruct_out=struct('hl',hl,'hs',hs,'evap',evap,'ch',tau.*rh,
%'ce',tau.*re,'tau',tau,'ssq',ssq,'xi',huol);

mystruct_out=struct('hl',hl,'hs',hs,'evap',evap,...
                    'ch',rd.*rh,'ce',rd.*re,'tau',tau,...
                    'ssq',ssq,'xi',huol, 'rd',rd,'re',re,'rh',rh,...
                    'ustar',ustar,'qstar',qstar,'tstar',tstar,...
                    'psimh',psimh,'psixh',psixh);



