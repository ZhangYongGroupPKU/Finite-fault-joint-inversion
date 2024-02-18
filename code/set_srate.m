function [srate] = set_srate(s1,s2)
if ismember(0,[s1,s2])
    if s1==s2
        srate = 1;
    else
        srate = max(s1,s2);
    end
else
    if s1==s2
        srate = s1;
    else
        error('The sample rate must be coincident!!');
    end
end
end
