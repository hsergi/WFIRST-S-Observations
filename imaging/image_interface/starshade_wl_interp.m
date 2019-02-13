function startshade_wl_interp( opt_in )

% Function to analyze the interpolation of WFIRST-S PSF across wavelength 

% Basic check
  if nargin == 0
  disp( '(STARSHADE_WL_INTERP) Assuming default options. Maybe you would like to set some. Continuing ...' )
  opt_in = [] ;
  end

% Get fields from te options
opt_in = get_default_options( opt_in ) ;

% Create the simulation if it does not exist
makeStarshadeImage( opt_in ) ;

% Beware that this file has some 'opt' cleared: Nx_image_pix, diam_image_mas 
load( [  opt_in.save_path opt_in.save_filename '.mat' ] )
% Generating the image field
efDefectImg = UTot2ImagePlane( lambdaIn, opt_in, pupil, telescopeDiameter, UTotL ) ;

% L/D (maximum)
l_ovr_d = lambdaIn( end ) * 1e9 / telescopeDiameter * 0.2036 ;
% Choosing some region around the source
n_l = 10 ;
b1 = opt.y_source_mas - l_ovr_d / 2 * n_l ;
b2 = opt.y_source_mas + l_ovr_d / 2 * n_l ;
a1 = opt.x_source_mas - l_ovr_d / 2 * n_l ;
a2 = opt.x_source_mas + l_ovr_d / 2 * n_l ;
% In pixels
px_per_mas = opt_in.Nx_img / opt_in.diam_img_mas ;
y1 = round( opt_in.Nx_img / 2 + b1 * px_per_mas ) ;
y2 = round( opt_in.Nx_img / 2 + b2 * px_per_mas ) ;
x1 = round( opt_in.Nx_img / 2 + a1 * px_per_mas ) ;
x2 = round( opt_in.Nx_img / 2 + a2 * px_per_mas ) ;

% Some checks, although this should not happen in most cases
  if y1 < 1, y1 = 1 ; end
  if y2 > opt_in.Nx_img, y2 = opt_in.Nx_img ; end
  if x1 < 1, x1 = 1 ; end
  if x2 > opt_in.Nx_img, x2 = opt_in.Nx_img ; end

n_lmbd = numel( lambdaIn ) ;
  for ii = 1 : n_lmbd
  int_src( ii, :, : ) = abs( efDefectImg( x1 : x2, y1 : y2, ii ) ).^2 ;
  end

% Interpolation
% Linear in different step sizes
  for n_fct = 1 : floor( n_lmbd / 2 ) ;
  clear dff_int_src
    for ii = n_fct + 1 : n_fct : n_lmbd - n_fct
    dff_int_src( ( ii - 1 ) / n_fct, :, : ) = 0.5 * ( int_src( ii - n_fct, :, : ) + int_src( ii + n_fct, :, : ) ) - int_src( ii, :, : ) ;
    end

  % Some plot:
  close all
  figure
  setwinsize( gcf, 1300, 750 )
  clf
  % Distribution of subplots
  n_dff = size( dff_int_src, 1 ) ;
  n_1 = sqrt( n_dff ) ;
    if round( n_1 ) == floor( n_1 ), n_1 = floor( n_1 ) ; n_2 = ceil( n_dff / n_1 ) ; else, n_1 = ceil( n_1 ) ; n_2 = floor( n_dff / n_1 ) ; end
  % Case of three subplots fails
    if n_dff == 3, n_1 = 2 ; n_2 = 2 ; end
    for i_dff = 1 : n_dff
    subplot( n_1, n_2, i_dff )
    imagesc( squeeze( dff_int_src( i_dff, :, : ) ) ) ; h = colorbar ; set( h, 'FontSize', 12 ) ;
    title( sprintf( '%3.0fnm', lambdaIn( i_dff ) * 1e9 ), 'FontSize', 14 )
    end
  suptitle( sprintf( 'WAVELENGTH INTERPOLATION RESIDUALS: LINEAR, BASE WITH DLAMBDA=%i nm. X_{SOURCE}=%3.1f, Y_{SOURCE}=%3.1f mas. PIXEL=%2.1f mas', opt_in.delta_lambda_nm * 2 * n_fct, opt_in.x_source_mas, opt_in.y_source_mas, 1 / px_per_mas ) )

  % Storing the image
    if ( opt_in.save_fig )
    img = getframe( gcf ) ;
    imwrite( img.cdata, sprintf( '%spanel_dlambda_%2.1f_x_%3.1f_y_%3.1f_pix_%2.2fmas_n_fct_%i.png', opt_in.save_path_fig, opt.delta_lambda_nm, ...
    opt_in.x_source_mas, opt_in.y_source_mas, 1 / px_per_mas, n_fct ) ) ;
    end
  end

dbstop if error
%make_an_error
