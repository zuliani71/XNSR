
function XYZ = sig_rot(XYZ,ALPHA,THETA)
% XYZ = sig_rot(XYZ,THETA,ALPHA) rotates
% the XYZ dataset by the ALPHA AZIMUTH ANGLE,
% and THETA DIP ANGLE. Both ALPHA and THETA
% are expressed in degrees.
% XYZ is a matrix of column vectors.
%
% e.g. sig_rot(XYZ,30,45); will roate the XYZ
% vecotors by a 30 degrees of azimuth and 45
% degrees of DIP
ALPHA=ALPHA/180*pi;
THETA=THETA/180*pi;
X=XYZ(:,1);
Y=XYZ(:,2);
Z=XYZ(:,3);
X_ROT = cos(THETA).*cos(ALPHA).*X - cos(THETA).*sin(ALPHA).*Y + sin(ALPHA).*Z;
Y_ROT = sin(ALPHA).*X + cos(ALPHA).*Y;
Z_ROT = -sin(THETA).*cos(ALPHA).*X + sin(ALPHA).*sin(THETA).*Y + cos(THETA).*Z;
XYZ=[X_ROT,Y_ROT,Z_ROT];