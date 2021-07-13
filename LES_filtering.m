%%
% Filtering data

% Grace Heaton
% University of Illinois at Urbana-Champaign
% Saxton-FOx Group @ UIUC

% Using subet of JHTDB channel flow data
% Subset created by Alan Hong
% Subset stored in data.mat as (U,V,W,X,Y,Z)

close all
clear all variables

load('data.mat') % subset of JHTDB channel data

sigma = 1; % standard dev based on channel half height (h=1 when nd)

%%
% Filter + define large scales 

UL = imgaussfilt3(U,sigma); % filtered U (large scales)
VL = imgaussfilt3(V,sigma); % filtered V (large scales)
WL = imgaussfilt3(W,sigma); % filtered W (large scales)


%%
% Extract small scales

US = U-UL; % small scale U
VS = V-VL; % small scale V
WS = W-WL; % small scale W 


%%
% Visualizations of small scale data

s1 = 1; % W location of 2D slice

figure
contourf(US(:,:,s1))
title('Small-scale U')
colorbar

figure
contourf(VS(:,:,s1))
title('Small-scale V')
colorbar

figure
contourf(WS(:,:,s1))
title('Small-scale W')
colorbar


%%
% Visualizations of filtered data (large, 2D slice)

s2 = 1; % W location of 2D slice
figure
contourf(UL(:,:,s2))
title('Filtered U')
colorbar

figure
contourf(VL(:,:,s2))
title('Filtered V')
colorbar

figure
contourf(WL(:,:,s2))
title('Filtered W')
colorbar

%%
% save the data
save('large_scales.mat','UL','VL','WL','X','Y','Z')
save('small_scales.mat','UL','VL','WL','X','Y','Z')
