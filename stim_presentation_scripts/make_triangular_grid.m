% make_triangular_grid.m

% takes in n_rows, n_cols
% <n_rows> is always a single number
% <n_cols> can be a single number, in which case a jittered triangular grid
%          subtending a rectangle is computed
%          if <n_cols> is a vector of length equal to n_rows, then each row
%          contains n_cols(row_num) points, centered 
%          (this can be used to create a hexagonal grid, as shown below
%          TODO)
% all coords are returned in normalized units, centered at (0,0) and
% largest absolute value of x and/or y is 1.
% 
% if n_cols is scalar, assumes the first (negative y) row is left (negative
% x) offset, then second is right (positive x) offset
%
% all language here assumes screen coordinates for PTB (negative up and to
% the left)



function [gridpts, rowid] = make_triangular_grid(n_rows,n_cols)


if length(n_cols)==1
    gridpts = nan(n_rows*n_cols,2);
elseif length(n_cols) == n_rows
    gridpts = nan(sum(n_rows*n_cols),2);
else
    % error: TODO
end



% if even number of rows, then offset will be 


% start w/ these "normalized" coords - can renormalize at the end
ypts = linspace(-1,1,n_rows);

% height = side_length*sqrt(3)/2
%
% so side_length = height (row-spacing) * 2 / sqrt(3)

x_sep = (ypts(2)-ypts(1))*2/sqrt(3);

startidx = 1;

for rr = 1:n_rows
    
    if length(n_cols)==1    
        tmp = (0:(n_cols-1))*x_sep;
    else
        tmp = (0:(n_cols(rr)-1))*x_sep;
    end
    
    tmp = tmp-mean(tmp);
    
    % if even number of x points, subtract x_sep/4
    
    if mod(rr,2) == 1
        tmp = tmp-x_sep/4;
    else
        tmp = tmp+x_sep/4;
    end
    
    if mod(length(tmp),2)==0
        tmp = tmp-x_sep/2;
    end
    
    thisidx = startidx:(startidx+length(tmp)-1);
    
    gridpts(thisidx,1) = tmp;
    gridpts(thisidx,2) = ypts(rr);
    
    startidx = thisidx(end)+1;
    clear thisidx tmp;
    
end

figure(1);clf;
plot(gridpts(:,1),gridpts(:,2),'ro');axis equal;



return
