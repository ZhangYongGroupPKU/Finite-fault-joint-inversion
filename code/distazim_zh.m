function da=distazim_zh(loca,epi)
%

% Supply WGS84 earth ellipsoid axis lengths in kilometers:
a = 6378.137; % definitionally
b = 6356.75231424518; % computed from WGS84 earth flattening coefficient definition
c = sqrt(a*a-b*b);
e=c./a;

[d,a]=distance('gc',epi,loca,[a,e]);
%[d,a]=distance(epi,loca,[a,e]);
da=[d,a];
return

