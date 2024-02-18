function [lonlabel,latlabel] = tick2label(lontick,lattick)
lonlabel = cell(length(lontick),1);
latlabel = cell(length(lattick),1);
for i = 1:length(lontick)
    a = abs(lontick(i));
    lonlabel{i} = [num2str(a+(360-2*a)*(a>180)),'\circ'];
end
if ismember(0,lontick(2:end-1))
    lonlabel{1} = [lonlabel{1},'W'];
    lonlabel{end} = [lonlabel{end},'E'];
elseif ismember(180,lontick(2:end-1))
    lonlabel{1} = [lonlabel{1},'E'];
    lonlabel{end} = [lonlabel{end},'W'];
elseif mean(lontick)>0;
    lonlabel{end} = [lonlabel{end},'E'];
else
    lonlabel{1} = [lonlabel{1},'W'];
end

for i = 1:length(lattick)
    latlabel{i} = [num2str(abs(lattick(i))),'\circ'];
end
if ismember(0,lattick(2:end-1))
    latlabel{1} = [latlabel{1},'S'];
    latlabel{end} = [latlabel{end},'N'];
elseif mean(lattick)>0;
    latlabel{end} = [latlabel{end},'N'];
else
    latlabel{1} = [latlabel{1},'S'];
end
end
