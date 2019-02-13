% mft_image
% Eric Cady, last rev. 9/7/11
%
% Uses the procedure described in "Fast computation of coronagraph
% propagation" by Soummer et al. to compute the electric field at the image
% plane following a pupil.
%
% Almost certainly slower than an FFT, but I like having the control.  Can
% put that in later.

function [Xout, Yout, imagePlane] = mft_image(imagePlaneDiameterInLambdaOverD,...
    numImagePoints, inputField)
% X and Y unitless (x/D)
% Xpp and Ypp in lambda/D

[x dx] = interval(1, size(inputField, 1));
[y dy] = interval(1,  size(inputField, 2));
[xpp dxpp] = interval(imagePlaneDiameterInLambdaOverD, numImagePoints);
[ypp dypp] = interval(imagePlaneDiameterInLambdaOverD, numImagePoints);

% Take to image plane
[Xout Yout] = meshgrid(xpp, ypp);
imagePlane = exp(-2*pi*1i*xpp.'*x)*(inputField)*exp(-2*pi*1i*y.'*ypp)*dx*dy;
make_an_error
function [x, dx] = interval(width, N)
dx = width/N;
x = -width/2 + dx/2:dx:width/2-dx/2;
