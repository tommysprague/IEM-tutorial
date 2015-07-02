function new_ylim = match_ylim(gobj,new_ylim)

% TODO:
% - if ax is figures, go through and get all the children and make them a
%   set of axes
% - if no input arguments, use get(0,'Children') to get handles



if nargin < 2
    
    tmp_ylim = get(gobj,'YLim');
    
    if iscell(tmp_ylim)
        tmp_ylim = cell2mat(tmp_ylim);
    end
    
    new_ylim = [min(tmp_ylim(:,1)) max(tmp_ylim(:,2))];
    
end

set(gobj,'YLim',new_ylim);



return