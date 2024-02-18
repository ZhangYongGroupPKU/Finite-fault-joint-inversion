%根据给定经纬度范围，输出crust1.0数据
function [oLon,oLat,odata] = imp(lLon,rLon,dLat,uLat)
load('crustall.mat');
if lLon > rLon
    warning('Longitude input wrong!');
end
if dLat > uLat
    warning('Latitude input wrong!');
end
n = length(Lon);
flag = zeros(n,1);
head = 0;
tail = 0;
for i = 1:n
    if dLat < Lat(i) && Lat(i) < uLat
        flag(i) = 1;
        if ~head
            head = i;
        end
    end
    if i > 2
        if flag(i-2) && ~flag(i-1)
            tail = i-2;
        end
    end
end
for i = head:tail
    if lLon < Lon(i) && Lon(i) < rLon
        flag(i) = 2;
    end
end
index = find(flag(:)==2);
lindex = length(index);
for i = 1:lindex
    oLon(i) = Lon(index(i));
    oLat(i) = Lat(index(i));
    obd1(i) = bd1(index(i));
    obd2(i) = bd2(index(i));
    obd3(i) = bd3(index(i));
    obd4(i) = bd4(index(i));
    obd5(i) = bd5(index(i));
    obd6(i) = bd6(index(i));
    obd7(i) = bd7(index(i));
    obd8(i) = bd8(index(i));
    obd9(i) = bd9(index(i));
    obd = [obd1;obd2;obd3;obd4;obd5;obd6;obd7;obd8;obd9];
    
    oro1(i) = ro1(index(i));
    oro2(i) = ro2(index(i));
    oro3(i) = ro3(index(i));
    oro4(i) = ro4(index(i));
    oro5(i) = ro5(index(i));
    oro6(i) = ro6(index(i));
    oro7(i) = ro7(index(i));
    oro8(i) = ro8(index(i));
    oro9(i) = ro9(index(i));
    oro = [oro1;oro2;oro3;oro4;oro5;oro6;oro7;oro8;oro9];
    
    ovp1(i) = vp1(index(i));
    ovp2(i) = vp2(index(i));
    ovp3(i) = vp3(index(i));
    ovp4(i) = vp4(index(i));
    ovp5(i) = vp5(index(i));
    ovp6(i) = vp6(index(i));
    ovp7(i) = vp7(index(i));
    ovp8(i) = vp8(index(i));
    ovp9(i) = vp9(index(i));
    ovp = [ovp1;ovp2;ovp3;ovp4;ovp5;ovp6;ovp7;ovp8;ovp9];
    
    ovs1(i) = vs1(index(i));
    ovs2(i) = vs2(index(i));
    ovs3(i) = vs3(index(i));
    ovs4(i) = vs4(index(i));
    ovs5(i) = vs5(index(i));
    ovs6(i) = vs6(index(i));
    ovs7(i) = vs7(index(i));
    ovs8(i) = vs8(index(i));
    ovs9(i) = vs9(index(i));
    ovs = [ovs1;ovs2;ovs3;ovs4;ovs5;ovs6;ovs7;ovs8;ovs9];
    
    odata = [obd,oro,ovp,ovs];
end
end
