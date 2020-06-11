function param = getParameters_ott(os, timePeriod, model, paramName)
% timePeriod: {'RD', 'cue', 'PE'}

tOfInt = ['mod_' timePeriod];
[~, param_ind] = intersect(os(1).(tOfInt).(model).paramNames, paramName);
if any(param_ind) == 0 % parameter not found
    error('Parameter not found')
elseif length(param_ind) > 1 % more than 1 parameter
    error('More than one parameter found with that name')
end

param = [];
for n = 1:length(os)
    param = [param os(n).(tOfInt).(model).bestParams(param_ind)];
end