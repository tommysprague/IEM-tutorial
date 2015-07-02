function new_clim = match_clim(gobj,new_clim)

% TODO:
% - if ax is figures, go through and get all the children and make them a
%   set of axes
% - if no input arguments, use get(0,'Children') to get handles



if nargin < 2
    
    tmp_clim = get(gobj,'CLim');
    
    if iscell(tmp_clim)
        tmp_clim = cell2mat(tmp_clim);
    end
    
    new_clim = [min(tmp_clim(:,1)) max(tmp_clim(:,2))];
    
end

set(gobj,'CLim',new_clim);



return