function m_out=m2m(m_in,flag)
%
if nargin==1
    if m_in<20
        %m_out=10.^((m_in+6.066666666667)*3/2);
        m_out=10.^((m_in+6.06).*3./2);
    else
        %m_out=log10(m_in)*2/3-6.066666666667;
        m_out=log10(m_in).*2./3-6.06;
    end
else
    if flag==1
        m_out=10.^((m_in+6.06).*3./2);
    else
        m_out=log10(m_in).*2./3-6.06;
    end
end