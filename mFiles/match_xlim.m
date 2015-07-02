function new_xlim = match_xlim(gobj,new_xlim)

% TODO:
% - if ax is figures, go through and get all the children and make them a
%   set of axes
% - if no input arguments, use get(0,'Children') to get handles



if nargin < 2
    
    tmp_xlim = get(gobj,'XLim');
    
    if iscell(tmp_xlim)
        tmp_xlim = cell2mat(tmp_xlim);
    end
    
    new_xlim = [min(tmp_xlim(:,1)) max(tmp_xlim(:,2))];
    
end

set(gobj,'XLim',new_xlim);



return