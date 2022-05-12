function x = perc_disc(dat, perc)

% apply "artifact correction" by excluding one or both ends 
% of distribution along 1st dimension in dat.
% perc = how much (in %) of distribution should be discard. 
% perc < 50% lower end of distribution will be excluded
% perc > 50% (100%-perc) of upper end of distribution will be excluded
% use two values for perc to discard both ends

if nargin == 1
    perc = 95;
end

if numel(perc) > 2
    error('2nd argument is a vector with max 2 values')
end

switch numel(perc)
    case 2
        sel = dat > prctile(dat,perc(1),1) & dat < prctile(dat,perc(2),1);
    case 1
        if perc < 50
            sel = dat > prctile(dat,perc,1);
        elseif perc > 50
            sel = dat < prctile(dat,perc,1);
        end
end

dat(~sel) = nan;
x = dat;

end

