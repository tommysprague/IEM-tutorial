%

function c = make_checkerboard(rad,sf,contrast)
% creates a circular checkerboard with parameters, in pixels (after
% conversion using deg2pix)
% returns a cell array of 2 matrices which can be turned into textures for
% use by PTB
% ex: stim = Screen('MakeTexture',w,c{1});

% GAMMA correction - this may need to be revised...
correctGamma = linspace(0,1,255);
LUT = correctGamma*255;

%make matrices x and y with meshgrid to hold pixel locations in terms
%of visual angle.
tmpX  = linspace(-rad,rad,rad*2);
[x,y] = meshgrid(tmpX);

%make a checkerboard image containing -1's and 1's.
chex = sign(sin(2*pi/sf*x).*sin(2*pi/sf*y));
circle = x.^2+y.^2<=(rad)^2; %p.stimSizePix

% first make the standard checkerboards
img1 = chex.*circle;
img2 = -1*img1; % contrast reversal

c{1} = LUT(round(((contrast*img1)+1)*127)+1);
%stims(1)=Screen('MakeTexture', w, tmpImg);

c{2} = LUT(round(((contrast*img2)+1)*127)+1);
%stims(2)=Screen('MakeTexture', w, tmpImg); % other checkerboard


return