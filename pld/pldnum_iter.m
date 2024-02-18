function n=pldnum_iter(b)
%n=pldnum_iter(b)
pb=abs(fft(b,2^nextpow2(length(b))));
nx=1-(min(pb)./max(pb)).^2;
n=log(0.0001)./log(nx);