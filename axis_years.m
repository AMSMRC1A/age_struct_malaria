function axis_years(ax, tfinal)
% convert x-axis units from days to years
% xticks([0:365:tfinal]);
temp = xticks;
set(ax,'xtickLabel',compose('%d',round(temp/365)));
end