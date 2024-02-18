function str=num2str_zh(num,dex)
str1=num2str(num,dex);
d=strfind(str1,'e');
if isempty(d)
    str=str1;
    return
end
%str1(d)='\\times';
%str1(d)='x';
str=[str1(1:d-1),'\times10^',str1(d+1),'^',str1(end)];