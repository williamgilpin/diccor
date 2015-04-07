function swapfig(h)
% allows the current figure to be updated without wating time constructing
% a new figure() object
set(0,'CurrentFigure',h)