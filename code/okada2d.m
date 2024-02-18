function  OU=okada2d(ALP,X,Y,DEP,AL1,AL2,AW1,AW2,...
    SD,CD,DISL1,DISL2,DISL3)
%C*********************************************************
%C*****                                               *****
%C*****    SURFACE DISPLACEMENT,STRAIN,TILT           *****
%C*****    DUE TO RECTANGULAR FAULT IN A HALF-SPACE   *****
%C*****       CODED BY Y.OKADA ... JAN 1985           *****
%C*****                                               *****
%C*********************************************************
%C
%C***** INPUT
%C*****   ALP     : MEDIUM CONSTANT  MYU/(LAMDA+MYU)=1./((VP/VS)**2-1)
%C*****   X,Y     : COORDINATE OF STATION, a raw vector
%C*****   DEP     : SOURCE DEPTH
%C*****   AL1,AL2 : FAULT LENGTH RANGE, a line vector
%C*****   AW1,AW2 : FAULT WIDTH RANGE
%C*****   SD,CD   : SIN,COS OF DIP-ANGLE
%C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)
%C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION
%C
%C***** OUTPUT
%C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL     )
%C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /
%C*****   U31,U32         : TILT                 UNIT OF X,Y,,,AW )
% disp(DISL2);
% disp(DISL1);
F0=0.;
F1=1.0;
%C-----
U = cell(3,1);
%for i=1:9,U{i}=X.*0; end

%
% raw = size(X,1);             % Get N, the number of points
% col = size(AW1,2);           % Get M, the number of patches of fault
%  X = repmat(X,1,col);       % reconstruc X as N*M
%  Y = repmat(Y,1,col);       % reconstruc Y as N*M
%AW1 = repmat(AW1,raw,1);
%AW2 = repmat(AW2,raw,1);
%AL2 = repmat(AL2,raw,1);
%AL1 = repmat(AL1,raw,1);
%DISL1 = repmat(DISL1,raw,1);
%DISL2 = repmat(DISL2,raw,1);
%DISL3 = repmat(DISL3,raw,1);
%  DEP = repmat(DEP,raw,1);
% FACTOR = repmat(1,col,1);
for i=1:3,U{i}=X.*0; end
P = Y.*CD + DEP.*SD;   % Calculation P, size as X & Y
Q = Y.*SD - DEP.*CD;   % Calculation Q, size as X & Y

for k=1:2
    if k == 1,ET = P-AW1; else ET = P-AW2;  end
    for j=1:2
        if j == 1,XI = X-AL1; else XI = X-AL2; end
        jk = j+k;
        if jk ~= 3, SIGN = F1; else SIGN =-F1; end
        XI = (XI~=F0).*XI+(XI==F0).*0.0001;
        ET = (ET~=F0).*ET+(ET==F0).*0.0001;
        DU = MSRECTG_OKADA(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3);
        for i=1:3,U{i}=U{i}+SIGN.*DU{i};end

    end

    OU = U(1:3);
end
%
function DU=MSRECTG_OKADA(ALP,XI,ET,Q,SD,CD,DISL1,DISL2,DISL3)

%C*********************************************************************
%C*****                                                           *****
%C*****  INDEFINITE INTEGRAL OF SURFACE DISPLACEMENT,STRAIN,TILT  *****
%C*****  DUE TO RECTANGULAR FAULT IN A HALF-SPACE                 *****
%C*****                          CODED BY  Y.OKADA ... JAN 1985   *****
%C*****                                                           *****
%C*********************************************************************
%C
%C***** INPUT
%C*****   ALP     : MEDIUM CONSTANT  MYU/(LAMDA+MYU)=1./((VP/VS)**2-1)
%C*****   XI,ET,Q : FAULT COORDINATE
%C*****   SD,CD   : SIN,COS OF DIP-ANGLE
%C*****          (CD=0.D0, SD=+/-1.D0 SHOULD BE GIVEN FOR VERTICAL FAULT)
%C*****   DISL1,DISL2,DISL3 : STRIKE-, DIP- AND TENSILE-DISLOCATION
%C
%C***** OUTPUT
%C*****   U1, U2, U3      : DISPLACEMENT ( UNIT= UNIT OF DISL    )
%C*****   U11,U12,U21,U22 : STRAIN       ( UNIT= UNIT OF DISL /
%C*****   U31,U32         : TILT                 UNIT OF XI,ET,Q )
%C
F0=0.;
F1=1.;
F2=2.;
PI2=6.283185307179586;
%
XI2=XI.^2;
ET2=ET.^2;
Q2=Q.^2;
R2=XI2+ET2+Q2;
R =sqrt(R2);
%R3=R.*R2;
D =ET.*SD-Q.*CD;
Y =ET.*CD+Q.*SD;
RET=R+ET;
RET(RET<F0)=F0;
RD =R+D;
%RRD=F1./(R.*RD);
%
TT = (Q ~=F0).*atan( XI.*ET./(Q.*R) )+(Q==F0).*F0;
%
RE = (RET ~= F0).*F1./RET+(RET==F0).*F0;
DLE = (RET~=F0).*log(RET)-(RET==F0).*log(R-ET);
RRX=F1./(R.*(R+XI));
RRE=RE./R;
%AXI=(F2.*R+XI).*RRX.*RRX./R;
%AET=(F2.*R+ET).*RRE.*RRE./R;
%if (CD == F0)  GO TO 20
if CD==F0
    RD2=RD.*RD;
    A1=-ALP./F2.*XI.*Q./RD2;
    A3= ALP./F2.*( ET./RD + Y.*Q./RD2 - DLE );
    A4=-ALP.*Q./RD;
    A5=-ALP.*XI.*SD./RD;
    %B1= ALP./F2.*  Q  ./RD2.*(F2.*XI2.*RRD - F1);
    %B2= ALP./F2.*XI.*SD./RD2.*(F2.*Q2 .*RRD - F1);
    %C1= ALP.*XI.*Q.*RRD./RD;
    %C3= ALP.*SD./RD.*(XI2.*RRD - F1);
else
    %C==============================
    %C=====   INCLINED FAULT   =====
    %C==============================
    TD=SD./CD;
    X =sqrt(XI2+Q2);
    %A5=(XI==F0).*XI.*F0+(XI~=F0).*XI.*ALP.*F2./CD.*atan((ET.*(X+Q.*CD)+X.*(R+X).*SD)./(XI.*(R+X).*CD));
    A5=(XI==F0).*F0+(XI~=F0).*ALP.*F2./CD.*atan((ET.*(X+Q.*CD)+X.*(R+X).*SD) ./ ((R+X).*XI.*CD) );
    A4= ALP./CD.*( log(RD) - SD.*DLE );
    A3= ALP.*(Y./RD./CD - DLE) + TD.*A4;
    A1=-ALP./CD.*XI./RD        - TD.*A5;
    %C1= ALP./CD.*XI.*(RRD - SD.*RRE);
    %C3= ALP./CD.*(Q.*RRE - Y.*RRD);
    %B1= ALP./CD.*(XI2.*RRD - F1)./RD - TD.*C3;
    %B2= ALP./CD.*XI.*Y.*RRD./RD      - TD.*C1;
end
%
A2=-ALP.*DLE - A3;
%B3=-ALP.*XI.*RRE - B2;
%B4=-ALP.*( CD./R + Q.*SD.*RRE ) - B1;
%C2= ALP.*(-SD./R + Q.*CD.*RRE ) - C3;
U1 =XI.*0;
U2 =XI.*0;
U3 =XI.*0;
%U11=XI.*0;
%U12=XI.*0;
%U21=XI.*0;
%U22=XI.*0;
%U31=XI.*0;
%U32=XI.*0;
%C======================================
%C=====  STRIKE-SLIP CONTRIBUTION  =====
%C======================================
UN=(DISL1~=F0).*DISL1./PI2+(DISL1==F0).*0;
REQ=RRE.*Q;
U1 =U1 - UN.*( REQ.*XI +   TT    + A1.*SD );
U2 =U2 - UN.*( REQ.*Y  + Q.*CD.*RE + A2.*SD );
U3 =U3 - UN.*( REQ.*D  + Q.*SD.*RE + A4.*SD );
%U11=U11+ UN.*( XI2.*Q.*AET - B1.*SD );
%U12=U12+ UN.*( XI2.*XI.*( D./(ET2+Q2)./R3 - AET.*SD ) - B2.*SD );
%U21=U21+ UN.*( XI.*Q./R3.*CD + (XI.*Q2.*AET - B2).*SD );
%U22=U22+ UN.*( Y .*Q./R3.*CD + (Q.*SD.*(Q2.*AET-F2.*RRE)-(XI2+ET2)./R3.*CD - B4).*SD );
%U31=U31+ UN.*(-XI.*Q2.*AET.*CD + (XI.*Q./R3 - C1).*SD );
%U32=U32+ UN.*( D.*Q./R3.*CD + (XI2.*Q.*AET.*CD - SD./R + Y.*Q./R3 - C2).*SD );
%C===================================
%C=====  DIP-SLIP CONTRIBUTION  =====
%C===================================
UN=(DISL2~=F0).*DISL2./PI2+(DISL2==F0).*0;
SDCD=SD.*CD;
U1 =U1 - UN.*( Q./R             - A3.*SDCD );
U2 =U2 - UN.*( Y.*Q.*RRX + CD.*TT - A1.*SDCD );
U3 =U3 - UN.*( D.*Q.*RRX + SD.*TT - A5.*SDCD );
%U11=U11+ UN.*( XI.*Q./R3            + B3.*SDCD );
%U12=U12+ UN.*( Y .*Q./R3 - SD./R     + B1.*SDCD );
%U21=U21+ UN.*( Y .*Q./R3 + Q.*CD.*RRE + B1.*SDCD );
%U22=U22+ UN.*( Y.*Y.*Q.*AXI - (F2.*Y.*RRX + XI.*CD.*RRE).*SD + B2.*SDCD );
%U31=U31+ UN.*( D .*Q./R3 + Q.*SD.*RRE + C3.*SDCD );
%U32=U32+ UN.*( Y.*D.*Q.*AXI - (F2.*D.*RRX + XI.*SD.*RRE).*SD + C1.*SDCD );
%C=========================================================================
%C===================  TENSILE-FAULT CONTRIBUTION  ========================
%C=========================================================================
UN=(DISL3~=F0).*DISL3./PI2+(DISL3==F0).*0;
SDSD=SD.*SD;
U1 =U1 + UN.*( Q2.*RRE                       - A3.*SDSD );
U2 =U2 + UN.*(-D.*Q.*RRX - SD.*(XI.*Q.*RRE - TT) - A1.*SDSD );
U3 =U3 + UN.*( Y.*Q.*RRX + CD.*(XI.*Q.*RRE - TT) - A5.*SDSD );
%U11=U11- UN.*( XI.*Q2.*AET             + B3.*SDSD );
%U12=U12- UN.*(-D.*Q./R3 - XI2.*Q.*AET.*SD + B1.*SDSD );
%U21=U21- UN.*( Q2.*(CD./R3 + Q.*AET.*SD) + B1.*SDSD );
%U22=U22- UN.*((Y.*CD-D.*SD).*Q2.*AXI - F2.*Q.*SD.*CD.*RRX- (XI.*Q2.*AET - B2).*SDSD );
%U31=U31- UN.*( Q2.*(SD./R3 - Q.*AET.*CD) + C3.*SDSD );
%U32=U32- UN.*((Y.*SD+D.*CD).*Q2.*AXI + XI.*Q2.*AET.*SD.*CD- (F2.*Q.*RRX - C1).*SDSD );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DU{1}=U1;%*FACTOR;
DU{2}=U2;%*FACTOR;
DU{3}=U3;%*FACTOR;
%DU{4}=U11;%*FACTOR;
%DU{5}=U12;%*FACTOR;
%DU{6}=U21;%*FACTOR;
%DU{7}=U22;%*FACTOR;
%DU{8}=U31;%*FACTOR;
%DU{9}=U32;%*FACTOR;


