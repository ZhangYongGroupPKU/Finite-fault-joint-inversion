function b=interm(a,x)
%b=interm(a,x)
sa=size(a);
if isscalar(x)
    b0=interp1(0:sa(1)-1,a,0:1/x:sa(1)-1);
    b00=interp1(0:sa(2)-1,b0',0:1/x:sa(2)-1);
    b=b00';
else
    b0=interp1(0:sa(1)-1,a,0:1/x(1):sa(1)-1);
    if sa(2)~=1
    b00=interp1(0:sa(2)-1,b0',0:1/x(2):sa(2)-1);
    b0=b00;
    end
    b=b0';
end