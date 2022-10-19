function [X,Y] = get_outer_roi(rad)


N = 100;  %  points of the inner circle
[x_c,y_c] = pol2cart(linspace(0,2*pi,N),3*ones(1,N));


X = [rad rad -rad -rad rad x_c];
Y = [0 rad rad -rad -rad y_c];


% patch(X,Y,1,'FaceColor','#bbbbbb','EdgeColor','none');