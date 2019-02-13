function plotStarshadeImage( opt_in, x_plt, y_plt )
% Function to plot the results of makeStarshadeImage
% History:
% 04/28/17: created, Sergi Hildebrandt (JPL/Caltech)
% 05/12/17: added interface with opt, Sergi Hildebrandt (JPL/Caltech)

% Get default options (look inside the function for specific definitions)
opt = get_default_options( opt_in ) ;
% Developer version? Here is just this
  if ( opt.developer ), units_image ; else, units ; end

% Main parameters
Nx = opt.Nx_pupil_pix ;
dlt_lmbd = opt.delta_lambda_nm ;
r_src = opt.r_source_mas ;
psi_src = opt.psi_source_deg ;
savePath = opt.save_path ;
  if ~isdir( savePath ), system( sprintf( 'mkdir -p %s', savePath ) ) ; end
% Saved files naming
saveFilename = sprintf( 'starshade_out_Nx_%i_pix_dl_%inm_dr_%3.1f_mas_psi_%3.1f_deg', Nx, dlt_lmbd, r_src, psi_src ) ;
% Case of an ideal pupil
  if strcmp( opt.pupil_file, '0' )
  saveFilename = sprintf( '%s_ideal', saveFilename ) ;
  end
savePathImg = opt.save_path_fig ;
  if ~isdir( savePathImg ), system( sprintf( 'mkdir -p %s', savePathImg ) ) ; end

% loading the stored data
load( sprintf( '%s/%s.mat', savePath, saveFilename ) ) ;

  set(0,'defaultlinelinewidth',1.0);
  set(0,'DefaultAxesFontSize',14);
  figure( 1 )
  clf; setwinsize(gcf,1000,600)
  % Size of the data to be plotted
  sz_img = size( efDefectImg ) ;
  n_lmbd = sz_img( 3 ) ;
  n_1 = floor( sqrt( n_lmbd ) ) ;
  n_2 = ceil( n_lmbd / n_1 ) ;
  % Fraction of the image to be plotted
  frc_img_plt = 0.25 ;
  a1 = sz_img( 1 ) * 0.5 * ( 1 - frc_img_plt ) + 1 ; % +1 to avoid index 0 in matlab if frc_img_plt would be 1
  a2 = sz_img( 1 ) * 0.5 * ( 1 + frc_img_plt ) ;
  b1 = sz_img( 2 ) * 0.5 * ( 1 - frc_img_plt ) + 1 ;
  b2 = sz_img( 2 ) * 0.5 * ( 1 + frc_img_plt ) ;
    for i_lmbd = 1 : n_lmbd
    subplot( n_2, n_1, i_lmbd )
    imagescnan( real( efDefectImg( a1 : a2, b1 : b2, i_lmbd ) ) ) ;
    title( sprintf( '%3.0f nm', lambdaIn( i_lmbd ) / nm ) )
    grid
    colorbar
    end
  nm_img = [ saveFilename '_real' ] ;
  suptitle( strrep( nm_img, '_', '\_' ) ) ;
  saveas( gca, sprintf( '%s/%s', savePathImg, nm_img ), 'png' ) ;

  figure( 2 )
  clf; setwinsize(gcf,1000,600)
  % Size of the images
  sz_img = size( efDefectImg ) ;
  n_lmbd = sz_img( 3 ) ;
  n_1 = floor( sqrt( n_lmbd ) ) ;
  n_2 = ceil( n_lmbd / n_1 ) ;
  % Fraction of the image to be plotted
  frc_img_plt = 0.25 ;
  a1 = sz_img( 1 ) * 0.5 * ( 1 - frc_img_plt ) + 1 ; % +1 to avoid index 0 in matlab if frc_img_plt would be 1
  a2 = sz_img( 1 ) * 0.5 * ( 1 + frc_img_plt ) ;
  b1 = sz_img( 2 ) * 0.5 * ( 1 - frc_img_plt ) + 1 ;
  b2 = sz_img( 2 ) * 0.5 * ( 1 + frc_img_plt ) ;
    for i_lmbd = 1 : n_lmbd
    subplot( n_2, n_1, i_lmbd )
    imagescnan( imag( efDefectImg( a1 : a2, b1 : b2, i_lmbd ) ) ) ;
    title( sprintf( '%3.0f nm', lambdaIn( i_lmbd ) / nm ) )
    grid
    colorbar
    end
  nm_img = [ saveFilename '_imag' ] ;
  suptitle( strrep( nm_img, '_', '\_' ) ) ;
  saveas( gca, sprintf( '%s/%s', savePathImg, nm_img ), 'png' ) ;

  % Sectional plots
  figure( 3 )
  clf; setwinsize(gcf,1200,750)
  subplot(221)
  % x axis
  hold all
    for i_plt = 1 : n_lmbd
    h_plt( i_plt ) = plot( real( efDefectImg( sz_img( 2 ) * 0.5, b1 : b2, i_plt ) ), '+-' ) ;
    lmbd_lgnd{ i_plt } = sprintf( '%3.0f nm', lambdaIn( i_plt ) / nm ) ;
    end
  xlabel( 'Pixel' ) ;
  ylabel( 'Real part' )
  title( 'Horizontal cut' )

  legend( h_plt, lmbd_lgnd, 'SouthWest' )
  legend boxoff

  subplot(222)
  % y axis
  hold all
    for i_plt = 1 : n_lmbd
    h_plt( i_plt ) = plot( real( efDefectImg( a1 : a2, sz_img( 1 ) * 0.5, i_plt ) ), '+-' ) ;
    end
  xlabel( 'Pixel' ) ;
  ylabel( 'Real part' )
  title( 'Vertical cut' )

  subplot(223)
  % x axis
  hold all
    for i_plt = 1 : n_lmbd
    h_plt( i_plt ) = plot( imag( efDefectImg( sz_img( 2 ) * 0.5, b1 : b2, i_plt ) ), '+-' ) ;
    end
  xlabel( 'Pixel' ) ;
  ylabel( 'Imaginary part' )
  title( 'Horizontal cut' )

  subplot(224)
  % y axis
  hold all
    for i_plt = 1 : n_lmbd
    h_plt( i_plt ) = plot( imag( efDefectImg( a1 : a2, sz_img( 1 ) * 0.5, i_plt ) ), '+-' ) ;
    end
  xlabel( 'Pixel' ) ;
  ylabel( 'Imaginary part' )
  title( 'Vertical cut' )
  grid

  nm_img = [ saveFilename '_cuts' ] ;
  suptitle( strrep( nm_img, '_', '\_' ) ) ;

  % Looking at the center of the image
  figure( 4 )
  clear h_plt
  clf; setwinsize(gcf,800,550)  
  % Creating the array of data to plot
    if ~exist( 'x_plt', 'var' )
    % Brightest spot
    tmp = squeeze( real( efDefectImg( :, :, 2 ) ) ) ;
    mx_plt = max( tmp( : ) ) ;
    q_mx = find( tmp == mx_plt ) ; % Avoiding the first image that has some aliasing
    q_mx = q_mx( 1 ) ;
    x_plt = floor( q_mx / sz_img( 1 ) ) + 1 ;
    end
    if ~exist( 'y_plt', 'var' )
    % Brightest spot
    y_plt = mod( q_mx, sz_img( 1 ) )  ;
    end
  
    for i_lmbd = 1 : n_lmbd
    plt_dt( i_lmbd ) = efDefectImg( y_plt, x_plt, i_lmbd ) ;
    % Some displacement
    y_plt_shft = y_plt + 5 ;
    x_plt_shft = x_plt + 5 ;
    plt_dt_shft( i_lmbd ) = efDefectImg( y_plt_shft, x_plt_shft, i_lmbd ) ;
    end
  
  % Real part
  subplot(211)
  plt_dt_tmp = real( plt_dt ) ;
  h_plt( 1 ) = plot( lambdaIn/nm, plt_dt_tmp, 'k+-' ) ;
  grid
  xlabel( 'Wavelength (nm)' ) ;
  ylabel( 'Real part' ) ;
  xlim( [ lambdaIn( 1 ) / nm - 50, lambdaIn( end ) / nm + 50 ] )
  % Interpolating with half the points
  lmbd_1 = lambdaIn( 1 : 2 : n_lmbd ) ;
  plt_dt_1 = plt_dt_tmp( 1 : 2 : n_lmbd ) ;
  plt_intrp = spline( lmbd_1, plt_dt_1, lambdaIn ) ;
  hold all
  h_plt( 2 ) = plot( lambdaIn/nm, plt_intrp, 'g-', 'LineWidth', 1.5 ) ;
  % Interpolating with one third of the points
  lmbd_1 = lambdaIn( 1 : 3 : n_lmbd ) ;
  plt_dt_1 = plt_dt_tmp( 1 : 3 : n_lmbd ) ;
  plt_intrp = spline( lmbd_1, plt_dt_1, lambdaIn ) ;
  h_plt( 3 ) = plot( lambdaIn/nm, plt_intrp, 'r-', 'LineWidth', 1.5 ) ;
  ttl_lgnd{ 1 } = 'Simulation' ;
  ttl_lgnd{ 2 } = 'Interpolation: 1/2' ;
  ttl_lgnd{ 3 } = 'Interpolation: 1/3' ;
  legend( h_plt, ttl_lgnd )
  % Imaginary part
  subplot(212)
  clear h_plt
  plt_dt_tmp = imag( plt_dt ) ;
  h_plt( 1 ) = plot( lambdaIn / nm, plt_dt_tmp, 'k+-' ) ;
  grid
  xlabel( 'Wavelength (nm)' ) ;
  ylabel( 'Imaginary part' ) ;
  xlim( [ lambdaIn( 1 ) / nm - 50, lambdaIn( end ) / nm + 50 ] )
  % Interpolating with half the points
  lmbd_1 = lambdaIn( 1 : 2 : n_lmbd ) ;
  plt_dt_1 = plt_dt_tmp( 1 : 2 : n_lmbd ) ;
  plt_intrp = spline( lmbd_1, plt_dt_1, lambdaIn ) ;
  hold all
  h_plt( 2 ) = plot( lambdaIn/nm, plt_intrp, 'g-', 'LineWidth', 1.5 ) ;
  % Interpolating with one third of the points
  lmbd_1 = lambdaIn( 1 : 3 : n_lmbd ) ;
  plt_dt_1 = plt_dt_tmp( 1 : 3 : n_lmbd ) ;
  plt_intrp = spline( lmbd_1, plt_dt_1, lambdaIn ) ;
  h_plt( 3 ) = plot( lambdaIn/nm, plt_intrp, 'r-', 'LineWidth', 1.5 ) ;
  legend( h_plt, ttl_lgnd )

  nm_img = [ saveFilename '_interp' ] ;
  suptitle( { strrep( nm_img, '_', '\_' ), sprintf( 'Pixel=(%i,%i)', x_plt, y_plt ) } ) ;
  saveas( gca, sprintf( '%s/%s', savePathImg, nm_img ), 'png' ) ;
  disp( sprintf( '(plotStarshadeImage) Images stored in %s', savePathImg ) ) ;  

  % Magnitude 
  close all ; setwinsize(gcf,900,250)
  plt_dt_tmp = abs( plt_dt ).^2 ;
  h_plt( 1 ) = plot( lambdaIn/nm, plt_dt_tmp, 'k+-' ) ;
  grid
  xlabel( 'Wavelength (nm)' ) ;
  ylabel( 'Intensity' ) ;
  xlim( [ lambdaIn( 1 ) / nm - 50, lambdaIn( end ) / nm + 50 ] )
  % Interpolating with half the points
  lmbd_1 = lambdaIn( 1 : 2 : n_lmbd ) ;
  plt_dt_1 = plt_dt_tmp( 1 : 2 : n_lmbd ) ;
  plt_intrp = spline( lmbd_1, plt_dt_1, lambdaIn ) ;
  hold all
  h_plt( 2 ) = plot( lambdaIn/nm, plt_intrp, 'g-', 'LineWidth', 1.5 ) ;
  % Interpolating with one third of the points
  lmbd_1 = lambdaIn( 1 : 3 : n_lmbd ) ;
  plt_dt_1 = plt_dt_tmp( 1 : 3 : n_lmbd ) ;
  plt_intrp = spline( lmbd_1, plt_dt_1, lambdaIn ) ;
  h_plt( 3 ) = plot( lambdaIn/nm, plt_intrp, 'r-', 'LineWidth', 1.5 ) ;
  ttl_lgnd{ 1 } = 'Simulation' ;
  ttl_lgnd{ 2 } = 'Interpolation: 1/2' ;
  ttl_lgnd{ 3 } = 'Interpolation: 1/3' ;
  legend( h_plt, ttl_lgnd )

  nm_img = [ saveFilename '_intensity_interp' ] ;
  suptitle( { strrep( nm_img, '_', '\_' ), sprintf( 'Pixel=(%i,%i)', x_plt, y_plt ) } ) ;
  saveas( gca, sprintf( '%s/%s', savePathImg, nm_img ), 'png' ) ;
  disp( sprintf( '(plotStarshadeImage) Images stored in %s', savePathImg ) ) ;

  % Magnitude with shifted pixels
  close all ; setwinsize(gcf,900,350)
  plt_dt_tmp = abs( plt_dt_shft ).^2 ;
  h_plt( 1 ) = plot( lambdaIn/nm, plt_dt_tmp, 'k+-' ) ;
  grid
  xlabel( 'Wavelength (nm)' ) ;
  ylabel( 'Intensity' ) ;
  xlim( [ lambdaIn( 1 ) / nm - 50, lambdaIn( end ) / nm + 50 ] )
  % Interpolating with half the points
  lmbd_1 = lambdaIn( 1 : 2 : n_lmbd ) ;
  plt_dt_1 = plt_dt_tmp( 1 : 2 : n_lmbd ) ;
  plt_intrp = spline( lmbd_1, plt_dt_1, lambdaIn ) ;
  hold all
  h_plt( 2 ) = plot( lambdaIn/nm, plt_intrp, 'g-', 'LineWidth', 1.5 ) ;
  % Interpolating with one third of the points
  lmbd_1 = lambdaIn( 1 : 3 : n_lmbd ) ;
  plt_dt_1 = plt_dt_tmp( 1 : 3 : n_lmbd ) ;
  plt_intrp = spline( lmbd_1, plt_dt_1, lambdaIn ) ;
  h_plt( 3 ) = plot( lambdaIn/nm, plt_intrp, 'r-', 'LineWidth', 1.5 ) ;
  ttl_lgnd{ 1 } = 'Simulation' ;
  ttl_lgnd{ 2 } = 'Interpolation: 1/2' ;
  ttl_lgnd{ 3 } = 'Interpolation: 1/3' ;
  legend( h_plt, ttl_lgnd )

  nm_img = [ saveFilename '_intensity_interp_shifted' ] ;
  suptitle( { strrep( nm_img, '_', '\_' ), sprintf( 'Pixel=(%i,%i)', x_plt_shft, y_plt_shft ) } ) ;
  saveas( gca, sprintf( '%s/%s', savePathImg, nm_img ), 'png' ) ;
  disp( sprintf( '(plotStarshadeImage) Images stored in %s', savePathImg ) ) ;



dbstop if error
%make_an_error

