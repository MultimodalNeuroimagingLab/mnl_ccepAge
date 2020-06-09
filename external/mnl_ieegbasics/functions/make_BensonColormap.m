function [cm,cm_angle,cm_eccen] = make_BensonColormap()

% This colormap was made using linspecer:
% Beautiful and distinguishable line colors + colormap
% version 1.4.0.0 (8.25 KB) by Jonathan C. Lansey
% and then made darker/lighter for a cortical surface
% and then shuffled to make sure adjacent maps differ in color
%
% dhermes, 2019 UMC Utrecht

% Benson_Area_Names = {'V1','V2','V3','hV4','V01','V02','V3a','V3b','LO1','LO2','TO1','T02'};

cm = [...
    0.9438    0.3910    0.2668
    0.1974    0.5129    0.7403
    0.5978    0.8408    0.6445
    0.9955    0.8227    0.4828
    0.3686    0.3098    0.6353
    0.8417    0.9409    0.6096
    0.6196    0.0039    0.2588
    0.8    0.8    0.4
    0.3539    0.7295    0.6562
    0.9877    0.6154    0.3391
    0.8206    0.2239    0.3094
    0.9    0.6    1];

% % Benson ELife 2018 colors areas:
% cm = [...
%     1,0,0;...%V1
%     0,1,0;...%V2
%     0,0,1;...%V3
%     0,1,1;...%hV4
%     0,0.5,1;...%V01
%     0,1,0.5;...%VO2
%     1,1,0;...'V3a'
%     0.5,1,0;...'V3b'
%     1,0,1;...'LO1'
%     0.5,0,1;...'LO2'
%     1,0,0.5;...'TO1'
%     0.5,0,0.5];%'TO2'

% Benson 2018 ELife colors angle, but one hemisphere from 1:180
cm_angle = zeros(180,3);
cm_angle(1:45,:) = [zeros(45,1),(0:1/44:1)',(0.5:0.5/44:1)'];
cm_angle(46:90,:) = [zeros(45,1),(1:-.5/44:.5)',(1:-1/44:0)'];
cm_angle(91:135,:) = [(0:1/44:1)',(0.5:0.5/44:1)',zeros(45,1)];
cm_angle(136:180,:) = [(1:-.5/44:.5)',(1:-1/44:0)',zeros(45,1)];

% Benson 2018 ELife colors eccentricity 
% In my rendering functions, the colors are rounded to integers, so in
% order to plot 1.25 as different from 1.5, we have to use a colormap that
% is 4x longer than the data values, here: 90*4 = 360. An eccentricity of 1
% should be plotted with the index of 4 etc..
cm_eccen = zeros(360,3);
cm_eccen(1:5,:) = [zeros(5,1),zeros(5,1),(0:0.5/4:.5)'];%1.25*4 = 5
cm_eccen(6:10,:) = [(0:1/4:1)',zeros(5,1),(.5:0.5/4:1)'];%2.5*4 = 10
cm_eccen(11:20,:) = [(1:-.5/9:.5)',zeros(10,1),(1:-1/9:0)'];%5*4 = 20
cm_eccen(21:40,:) = [(.5:.5/19:1)',(0:1/19:1)',zeros(20,1)];%10*4 = 40
cm_eccen(41:80,:) = [(1:-1/39:0)',(1:-0.5/39:0.5)',zeros(40,1)];%20*4 = 80
cm_eccen(81:160,:) = [zeros(80,1),(.5:0.5/79:1)',(0:1/79:1)'];%40*4 = 160
cm_eccen(161:360,:) = [(0:1/199:1)',ones(200,1),ones(200,1)];%161:360

             