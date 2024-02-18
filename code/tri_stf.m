function stf = tri_stf(len1,len2)
% =================================================================================
% stf = tri_stf(len1,len2) 
%     TRI_STF is to generate a triang function with rising time of len1 and
%     descending time of len2。 Get the result as stf ,with the length 
%     len1+len2+1
% -------------------------------------------------------------------------
%                                    Zhang Yong Xu Lisheng and Chen Yuntai
%                                             2006/10/14/15:50  地球物理所
% ==================================================================================
% stf1=Bartlett(2*len1+1);
% stf2=Bartlett(2*len2+1);
% stf=[stf1(1:len1+1);stf2(end-len2+1:end)];
stf=[0:1/len1:1,1-1/len2:-1/len2:0];
stf=stf(:);
return