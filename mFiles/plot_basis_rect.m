% plot_basis_rect.m

% adapted from plot_basis by TCS 10/25/13
%
% plot_basis_rect(b),n_rfX,n_rfY;

function plot_basis_rect(b,n_rfX,n_rfY,resX,resY)

figure; clf; 

nr = n_rfY;
nc = n_rfX;

%res = sqrt(size(b,2));

ridx = 1; cidx = 1;

for bb = 1:size(b,1);
    
    % goes down each column first
    
    subplot(nr,nc,cidx + (ridx-1)*nc);
    %subplot(nr,nc,bb);
    imagesc(reshape(b(bb,:),resY,resX));
    set(gca,'XTick',[],'YTick',[]);
    axis equal;
    axis tight;
    if ridx == nr
        ridx = 1;
        cidx = cidx+1;
    else
        ridx = ridx+1;
    end
    
end

figure;
%imagesc(reshape(sum(b,1),resY,resX));
surf(reshape(sum(b,1),resY,resX));
%set(gca,'XTick',[],'YTick',[]); axis square;


return