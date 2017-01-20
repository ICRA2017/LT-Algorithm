%========================================================================
%   plotWindMap
%   version 1.0 - January 18th, 2017
%   
%   Plots an estimated expected wind map.   
%   Inputs:
%   xMat, yMat: Grid mesh of the x,y coordinates of the prediction
%   locations
%   uMat, vMat: U,V components of the predicted airflow vectors.
%   cmap: color map used to display the expected airflow speed
%   scale_val: scale used to display the airflow arrows.
%========================================================================



function plotWindMap(xMat,yMat,uMat,vMat,cmap,arrow_color,scale_val)
%N_COLORS=64;
[rows cols]=size(uMat);
wnd_spd=sqrt(uMat.*uMat+vMat.*vMat);
angMat=atan2(vMat,uMat);
[c,h]=contourf(xMat,yMat,wnd_spd,255);
set(h,'EdgeColor','none');
colormap(cmap)
colorbar;
if scale_val(1)~=scale_val(2)
    caxis(scale_val);
end
hold on;
qvr=quiver(xMat,yMat,wnd_spd.*cos(angMat),wnd_spd.*sin(angMat),arrow_color);
set(qvr,'linewidth',2.5);

