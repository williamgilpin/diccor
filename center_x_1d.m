function cx = center_x_1d(x);

% CENTER_X_1D calculates the x-position of the peak of the 
%	cross-correlation function x by fitting the nearest 
% 	neighbors (left and right) of the maximum value of x to a parabola.

try;
   
[mi,mj]	= find(x == max(max(x)));
mi			= mi(1);
mj			= mj(1); 

%whos x(mi,mj)
%maxcorx=x(mi,mj);
%user_entry = input('prompt') 

[xobj, yobj]	= 	meshgrid(1 : size(x,2), 1 : size(x,1));

x_fit = xobj(mi, mj-1 : mj+1);
y_fit = x(mi, mj-1 : mj+1);

[a,b] = polyfit(x_fit,y_fit,2);

x_vals = [min(x_fit) : 0.01 : max(x_fit)];
y_vals = polyval(a, x_vals);

cx = x_vals(find(y_vals == max(y_vals))) - size(x,2)/2 - 1;

if length(cx) > 1, cx = mean(cx); end;

% %%%% DELETTTTEEEEEEEEEE this
% [cx,cy] = x_cntrcorr(x,mj,mi,2);
% cx= cx- size(x,2)/2 - 1;
maxcorx=x(mi,mj);

if maxcorx < 0.5;, cx =0; end;

catch;
   cx = NaN;
end;
