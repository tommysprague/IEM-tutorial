% IEM_ori_fMRI.m
%
% Tommy Sprague, 3/13/2015, originally for Bernstein Center IEM workshop
% tommy.sprague@gmail.com or jserences@ucsd.edu
%
% Uses a dataset graciously donated by Mary Smith, results form which are 
% to be presented at VSS 2015 (maryesmith@gmail.com).
%
% Please DO NOT use data from these tutorials without permission. Feel free
% to adapt or use code for teaching, research, etc. If possible, please
% acknowledge relevant publications introducing different components of these
% methods (Brouwer & Heeger, 2009; 2011; Sprague & Serences, 2013; 
% Sprague, Saproo & Serences, 2015)
%
% Previously (IEM_spatial.m) we leared about how to implement the IEM
% analysis pipeline for a "spatial" feature space. This method also works
% for other feature spaces, such as color (Brouwer & Heeger, 2009; 2013, J
% Neuro), orientation (Brouwer & Heeger, 2011, J Neurophys; Scolari et
% al, 2012, J Neuro; Ester et al, 2013, J Cog Neuro; Anderson et al, 2013, 
% J Neuro; Byers & Serences, 2014, J Neurophys), and motion direction
% (Saproo & Serences, 2014, J Neuro). Here, we'll focus on a dataset
% acquired while participants viewed visual stimuli of different
% orientations. On each trial, each hemifield was presented with an
% randomly-chosen orientation (0:160 deg, measured clockwise
% from vertical). On each trial, participants were cued to attend to either
% the left grating or the right grating. Let's try and reconstruct channel
% response functions for stimulus orientation across the attention
% conditions.

close all;

root = load_root();%'/usr/local/serenceslab/tommy/berlin_workshop/';
addpath([root 'mFiles/']);

%% load data
%
% Each data file contains:
% - myLBetas and myRBetas: n_trials x n_voxels signals extracted from each
%   voxel of LH & RH V1, respectively
% - myOrL and myOrR: n_trials x 1 orientation labels for the orientation
%   presented on the left and right HEMISPHERES, respectively (opposite
%   side of screen; arranged this way for convenience)
% - myAttnSide: n_trials x 1, label for which side was attended (L or R, 1
%   or 2, respectively)
%

load([root 'fMRI_ori/s03_data.mat']);

% Each data file contains 8 scanning runs, with 36 trials each. Let's create a
% variable which keeps track of run number for our cross-validation later
% on.

myScanNum = [];
for ii = 1:8
    myScanNum = [myScanNum; ii*ones(36,1)];
end


%% generate orientation channels (orientation filters)
%
% This step of the analysis is very similar to the corresponding step in
% the spatial version of the IEM. However, instead of using spatial filters
% to generate our predicted channel responses during weight estimation,
% we're going to use orientation tuning functions (TFs), modeled as
% sinusoids raised to a high power. 
%
% QUESTION: How many unique orientation channels can we model? Fill in the
% answer below, it will be used for the remainder of the script.

n_ori_chans = 9;

% each orientation channel can be modeled as a "steerable filter" (see
% Freeman and Adelson, 1991), which here is a (co)sinusoid raised to a high
% power (specifically, n_ori_chans-1). This particular shape isn't terribly
% important (I don't think, though others may disagree) - you'll get
% equally far modeling responses as Gaussians/von Mises functions, but
% we'll stick with the steerable filters here to keep things consistent w/
% our lab's prior work (and that of Gijs Brouwer & David Heeger). 

make_basis_function = @(xx,mu) (cosd(xx-mu)).^(n_ori_chans-mod(n_ori_chans,2));

% QUESTION: what is the bandwidth of each basis function (FWHM)? Can you
% write an expression which takes any arbitrary power and computes FWHM? 

% Let's look at our channels.

xx = linspace(1,180,180);
basis_set = nan(180,n_ori_chans);
chan_center = linspace(180/n_ori_chans,180,n_ori_chans);

for cc = 1:n_ori_chans
    basis_set(:,cc) = make_basis_function(xx,chan_center(cc));
end

figure;
subplot(1,2,1);plot(xx,basis_set);
xlabel('Orientation (\circ)');
ylabel('Filter amplitude');
title('Basis set');

% now let's see how well the channels tile the space by summing them:
subplot(1,2,2);plot(xx,sum(basis_set,2));


%% Voxel selection (later)
%
% For this tutorial, let's assume we'll just use all voxels in each region.
% In practice, it can often improve data to select subsets of voxels, such
% as those which have the highest Mutual Information (MI) between
% activation level and feature value, or those which have the highest
% F-score for a 1-factor ANOVA with orientation as the single factor.
% Alternatively, independent localizer scans (such as those which show a
% high-contrast checkerboard pattern subtending the relevant portions of
% the visual field) may be used to select relevant voxels for analysis.
%
% If you would like to try different methods for voxel selection, implement
% them below. You should choose the set of voxels based only on the
% training data - then use those voxels to compute channel weights (using training
% data) and orientation reconstructions (using testing data from voxels
% selected using only training data)




%% Build stimulus mask
%
% This is a very similar process to what we did for the spatial IEM
% dataset. Here, instead of a spatial "mask" covering the portion of the
% screen occupied by the visual stimulus/stimuli, we're creating an "image"
% of the portion of the feature space (orientation) occupied by the visual
% stimulus/stimuli for each set of voxels. Because we're analyzing each
% hemisphere separately here, we'll only consider the relevant orientation.
%
% The stimulus "image" should be the same dimensionality as the basis
% function vectors (here, 180 points). Because we're using square-wave,
% single-orientation gratings (that is, they have a single orientation
% value, and are not superimposed w/ other orientations as could be the
% case w/ plaids or other complex stimuli), we can consider the stimulus
% images to be a vector of 0's, w/ a single value of 1 at the stimulus
% orientation. Note that the basis set vectors and stimulus image can each
% be as high- or low-resolution as you like (down to n_channels values).
% Here, for convenience, we've simply used 180 points, 1 point for each
% degree.

% go from orientation label (1:9) to orientation value (20:180)
myOrL = myOrL-1;myOrR = myOrR-1;
myOrL(myOrL==0)=9; myOrR(myOrR==0)=9;
myOrL = myOrL*20; myOrR = myOrR*20;
% if you like, you could adjust these values to cartesian orientations (0
% at +x, + orientation goes CCW) - but for simplicity I'll leave them as
% they were by default.

stim_mask_L = zeros(size(myOrL,1),length(xx));
stim_mask_R = zeros(size(myOrL,1),length(xx));
for tt = 1:size(stim_mask_L,1)  % loop over trials
    stim_mask_L(tt,myOrL(tt))=1;
    stim_mask_R(tt,myOrR(tt))=1;
end


% let's check our stimulus mask
figure;
subplot(1,2,1);imagesc(stim_mask_L);title('Stimulus mask (left)');xlabel('Orientation (\circ)');ylabel('Trial');
subplot(1,2,2);imagesc(stim_mask_R);title('Stimulus mask (right)');xlabel('Orientation (\circ)');ylabel('Trial');



%% Generate design matrix
%
% We're modeling everything the same here as we did for IEM_spatial.m -
% each voxel is modeled as a linear combination of a set of hypothetical
% neural populations, or information channels, which are modeled as
% orientation TFs (see above). The set of TFs forms our "basis set" (see
% above), which is used to generate predictions for how each channel should
% respond to any given stimulus value. Let's generate those predictions as
% a matrix of trials x channels, called a "design matrix"
%
% Because we've computed "stimulus images", we just need to project each
% stimulus image onto the n_ori_chans basis functions. 
%
% QUESTION: What is special about this particular type of projection? Can
% you come up with a different way to generate a design matrix under
% assumptions that a single orientation is always presented at one of a
% discrete number of values? 

trnX_L = stim_mask_L*basis_set;
trnX_R = stim_mask_R*basis_set;

% As before, we need to be sure our design matrix is full rank (that is,
% all channels are non-colinear, and that we can uniquely estimate the
% contribution of each channel to the signal observed in a given voxel). 

fprintf('Design matrix L: rank = %i\nDesign matrix R: rank = %i\n',rank(trnX_L),rank(trnX_R));

% QUESTION: Which factors determine the rank of the design matrix? What
% happens to the rank as you change the number of orientation channels? The
% number of orientations included in the analysis?


% let's look at the predicted response for a single trial

tr_num = 5;
figure;
subplot(1,4,1);hold on;
plot(xx,basis_set);
plot(xx,stim_mask_L(tr_num,:),'k-');
xlabel('Orientation (\circ)');title(sprintf('Left Stimulus, trial %i',tr_num));
xlim([0 180]);
subplot(1,4,2);hold on;
plot(chan_center,trnX_L(tr_num,:),'k-');
for cc = 1:n_ori_chans
    plot(chan_center(cc),trnX_L(tr_num,cc),'o','MarkerSize',8,'LineWidth',3);
end
xlabel('Channel center (\circ)');title('Predicted channel response');
xlim([-5 180]);

% and each design matrix

subplot(1,4,3);hold on;
imagesc(chan_center,1:size(trnX_L,1),trnX_L);
title('Design matrix (L)');
xlabel('Channel center (\circ)');ylabel('Trial'); axis tight ij;

subplot(1,4,4);hold on;
imagesc(chan_center,1:size(trnX_R,1),trnX_R);
title('Design matrix (R)');
xlabel('Channel center (\circ)'); axis tight ij;



%% Cross-validate and train/test encoding model within each hemisphere
%
% As with any machine learning or model-fitting it's important to be sure
% the data used to estimate a model (the weights for each voxel) is
% entirely independent of the data used to evaluate the model's efficacy,
% or here, to reconstruct stimulus representations. Because we z-scored the
% timeseries of each voxel within each run during preprocessing, data
% within a run is not independent. So here, we'll implement a
% leave-one-run-out cross-validation scheme for encoding model estimation
% and reconstruction. We'll use runs 1-(n-1) for model estimation, then map
% activation patterns from run n into the feature space using the map
% computed from runs 1-(n-1).

vox_prctile = 50; % use top this % of voxels from each hemisphere

ru = unique(myScanNum);
n_runs = length(ru);

chan_resp_L = nan(size(trnX_L));
for rr = 1:n_runs
    
    % identify the training & testing halves of the data
    trnIdx = myScanNum~=ru(rr);
    tstIdx = myScanNum==ru(rr);
    
    % voxel selection - compute ANOVA on each voxel, rank them, select top
    % n% of voxels
    ps = nan(size(myLbetas,2),1);
    for ii = 1:size(myLbetas,2)
        ps(ii) = anova1(myLbetas(trnIdx,ii),myOrL(trnIdx),'off');
    end
    
    % if using all voxels
    %which_vox = ones(size(myLbetas,2))==1;
    
    which_vox = ps <= prctile(ps,50);
    
    trnL = myLbetas(trnIdx,which_vox);
    tstL = myLbetas(tstIdx,which_vox);

    % train the encoding model (remember, each hemisphere is separate here,
    % so must do this for each hemisphere seperately - here, I'm just doing
    % LH, you can add RH)
    w_L = trnX_L(trnIdx,:)\trnL;
    
    chan_resp_L(tstIdx,:) = (inv(w_L*w_L')*w_L*tstL').';
    
end



%% Combine channel response functions across orientations
%
% The easiest way to do this is to circularly shift the channel response
% functions so that their centers are all aligned (say, at channel 4). 

ou = unique(myOrL);

targ_chan = ceil(n_ori_chans/2);

chan_resp_shift_L = nan(size(chan_resp_L));
for ii = 1:size(chan_resp_L,1)
    chan_resp_shift_L(ii,:) =  circshift(chan_resp_L(ii,:),targ_chan-find(ou==myOrL(ii)),2);
end

figure;
hold on;
plot(chan_center,    mean(chan_resp_shift_L(myAttnSide==1,:)),'--');
plot(chan_center,    mean(chan_resp_shift_L(myAttnSide==2,:)),'-');
xlabel('Orientation channel (\circ)');
ylabel('Channel response (BOLD Z-score)');
title('LH coregistered orientation reconstructions');

% now do all the model training/testing, plotting, etc for the other
% hemisphere. if you use cell arrays, there's a more efficient way to set
% this up than splitting into LH/RH variables!