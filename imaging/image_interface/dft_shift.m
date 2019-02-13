% Function that shifts using DFT
function img_shft = dft_shift( img_in, dlt_x, dlt_y )

dft_in = fftshift( fft2( img_in ) ) ;
n_in = size( img_in, 1 ) ;
% In case it's an odd number
n_in_hlf = floor( n_in / 2 ) ;
% Basic check
  if size( img_in, 2 ) ~= n_in
  disp( '(DFT_SHIFT) Dimensions of the image are different. It ought to be a square. Stopping.' )
  dbstop if error
  make_an_error
  end

% Define shift in frequency domain
[ xF, yF ] = meshgrid( -n_in_hlf : n_in - n_in_hlf - 1, -n_in_hlf : n_in - n_in_hlf - 1 ) ;

% Perform the shift
dft_in = dft_in .* exp( -1i * 2 *pi .* ( xF * dlt_x + yF * dlt_y) / n_in ) ;

% Find the inverse Fourier Transform (imaginary part is zero or within numerical precison from zero)
img_shft = real( ifft2( ifftshift( dft_in ) ) ) ;

