function starshade_translation( opt_in )

% Reference source to comapre to
  if ~exist( 'i_src_dft', 'var' )
  i_src_dft = 1 ;
  end

% Some sources (mas from the center)
src_arry = [ 400, 350, 305, 304, 303, 302.5, 302, 301, 300.5, 300, 202.5, 202, 200.5, 200.1, 200, 152.5, 150.5, 150, 102.5, 102, 101.5, 101, 100, 70, 60, 50, 40, 30 ] ;
n_src_arry = numel( src_arry ) ;
x_src( 1 : n_src_arry ) = src_arry ; x_src( n_src_arry + 1 : 2 * n_src_arry ) = 0 ; x_src( 2 * n_src_arry + 1 : 3 * n_src_arry ) = src_arry / sqrt( 2 ) ; x_src( 3 * n_src_arry + 1 : 4 * n_src_arry ) = src_arry * 151 / sqrt( 151^2 + 17^2 ) ;
y_src( 1 : n_src_arry ) = 0 ; y_src( n_src_arry + 1 : 2 * n_src_arry ) = src_arry ; y_src( 2 * n_src_arry + 1 : 3 * n_src_arry ) = src_arry / sqrt( 2 ) ; y_src( 3 * n_src_arry + 1 : 4 * n_src_arry ) = src_arry * 17 / sqrt( 151^2 + 17^2 ) ;
n_src=  numel( x_src ) ;
  if numel( y_src ) ~= n_src
  disp( '(starshade_translation) Number of elements for x and y arrays is different. Stopped.' )
  return
  end
% subgroups
alpha_lst = [ 0, 90, 45 ] ;
n_grp = numel( alpha_lst ) ;
  for i_grp = 1 : n_grp
  grp_lst( i_grp ) = n_src_arry ;
  end
% Appending the last case, which is atan2(17,151)=6.42 deg
alpha_lst( end + 1 ) = 6.42 ;
grp_lst( end + 1 ) = n_src_arry ;
n_grp = n_grp + 1 ;

opt_in.developer = 1 ;
% Only one wavelength
%opt_in.delta_lambda_nm = 200 ;
opt_in.redo = 0 ;
opt_in.save = 1 ;
% Only one wavelength
opt_in = get_default_options( opt_in ) ;
i_1 = 1 ;
i_2 = n_src_arry ; % n_src (except for the example with 0.5 mas/pix which has too much memory overall)
  for i_src = i_1 : i_2
  opt_in.x_source_mas = x_src( i_src ) ;
  opt_in.y_source_mas = y_src( i_src ) ;
  makeStarshadeImage( opt_in ) ;
  end

% Loading each simulation (notation from makeStarshadeImage)
dr_out = 'out_dev/' ;
  for i_src = i_1 : i_2
  % Creating the filename
  opt_in.x_source_mas = x_src( i_src ) ;
  opt_in.y_source_mas = y_src( i_src ) ;
  opt_in = get_default_options( opt_in ) ;
  load( [  opt_in.save_filename '.mat' ] ) 
  % Generating the image field
  efDefectImg = UTot2ImagePlane( lambdaIn, opt_in, pupil, telescopeDiameter, UTotL ) ;
  % For plotting purposes
    if i_src == i_src_dft
    x_source_mas_rf = opt_in.x_source_mas ;
    y_source_mas_rf = opt_in.y_source_mas ;
    end
  n_lmbd = numel( lambdaIn ) ;
    for i_lmbd = 1 : n_lmbd
   int_src( i_src, i_lmbd, :, : ) = abs( squeeze( efDefectImg( :, :, i_lmbd ) ) ).^2 ;
    % Computing the DFT of the reference source
      if i_src == i_src_dft
      dft_rf( i_lmbd, :, : ) = fft2( squeeze( int_src( i_src, i_lmbd, :, : ) ) ) ;
      end
    end
  end

% Testing the translation

% Simple circular shift (assuming the source positions correspond to integer values of pixel numbers)
%% Reference source is the first one
fct_img_1 = opt_in.Nx_img / 400 ;
mas_pr_px = opt_in.diam_img_mas / opt_in.Nx_img ;
clear src_dff
  for i_src_rf = 1 : 1 % temporary n_src ;
  a1 = opt_in.Nx_img / 2 + y_src( i_src_rf ) / mas_pr_px - 50 * fct_img_1 ;
  a2 = opt_in.Nx_img / 2 + y_src( i_src_rf ) / mas_pr_px + 50 * fct_img_1 ;
  b1 = opt_in.Nx_img / 2 + x_src( i_src_rf ) / mas_pr_px - 50 * fct_img_1 ;
  b2 = opt_in.Nx_img / 2 + x_src( i_src_rf ) / mas_pr_px + 50 * fct_img_1 ;
    for i_src = 1 : i_2 - i_1 + 1
    idx_shft_x = ( x_src( i_src_rf ) - x_src( i_src ) ) / mas_pr_px ;
    idx_shft_y = ( y_src( i_src_rf ) - y_src( i_src ) ) / mas_pr_px ;
      for i_lmbd =1 : n_lmbd
      % This is using the integer shifts in pixels
      %tmp = circshift( squeeze( int_src( i_src, i_lmbd, :, : ) ), [ ceil( idx_shft_x ) ceil( idx_shft_y ) ] ) - squeeze( int_src( i_src_rf, i_lmbd, :, : ) ) ;
      tmp = dft_shift( squeeze( int_src( i_src, i_lmbd, :, : ) ), idx_shft_y, idx_shft_x ) - squeeze( int_src( i_src_rf, i_lmbd, :, : ) ) ;
      src_dff( i_src_rf, i_src, i_lmbd, :, : ) = tmp( round( b1 ) : round( b2 ),  round( a1 ) : round( a2 ) ) ;
      end
    end
  end

% Reduced to the standard size
%% Two reductions before & after
%% Commented out for now
if ( 0 )
clear src_dff_0
  for i_src_rf = 1 : n_src ;
    for i_src = 1 : n_src
    idx_circ_x = ceil( ( x_src( i_src_rf ) - x_src( i_src ) ) / mas_pr_px / fct_img_1 ) ;
    idx_circ_y = ceil( ( y_src( i_src_rf ) - r_src( i_src ) ) / mas_pr_px / fct_img_1 ) ;
      for i_lmbd =1 : n_lmbd
      a1 = round( opt_in.Nx_img / 2 / fct_img_1 + y_src( i_src_rf ) / mas_pr_px / fct_img_1 - 50 ) ;
      a2 = round( opt_in.Nx_img / 2 / fct_img_1 + y_src( i_src_rf ) / mas_pr_px / fct_img_1 + 50 ) ;
      b1 = round( opt_in.Nx_img / 2 / fct_img_1 + x_src( i_src_rf ) / mas_pr_px / fct_img_1 - 50 ) ;
      b2 = round( opt_in.Nx_img / 2 / fct_img_1 + x_src( i_src_rf ) / mas_pr_px / fct_img_1 + 50 ) ;
      tmp = circshift( imresize( squeeze( int_src( i_src, i_lmbd, :, : ) ), 1 / fct_img_1 ), [ idx_circ_x idx_circ_y ] ) - imresize( squeeze( int_src( i_src_rf, i_lmbd, :, : ) ), 1 / fct_img_1 ) ;
      % Before
      src_dff_0( 1, i_src_rf, i_src, i_lmbd, 1 : b2 - b1 + 1, 1 : a2 - a1 + 1 ) = tmp( b1 : b2, a1 : a2 ) ;
      % after
      tmp = imresize( squeeze( src_dff( i_src_rf, i_src, i_lmbd, :, : ) ), 1 / fct_img_1 ) ;
      src_dff_0( 2, i_src_rf, i_src, i_lmbd, 1 : size( tmp, 1 ), 1 : size( tmp, 2 ) ) = tmp;
      end
    end
  end
end % if ( 0 )

% Some plotting
opt_in.log10 = 0 ;
%% Reduced resolution
lbl_rd = { 'before shifting', 'after shifting' } ;
if ( 0 )
  for i_rd = 1 : 2
    for i_src_rf = 1 : n_src
    figure( ( i_rd - 1 ) * n_src + i_src_rf )
    setwinsize( gcf, 1300, 750 )
    clf
    i_lmbd = 1 ;
    r_src_rf = sqrt( x_src( i_src_rf ) * x_src( i_src_rf ) + y_src( i_src_rf ) * y_src( i_src_rf ) ) ;
      for i_src_plt = 1 : n_src
      subplot( 4, 4, i_src_plt )
        if ( opt_in.log10 )
        imagesc( log10( abs( squeeze( src_dff_0( i_rd, i_src_rf, i_src_plt, i_lmbd, :, : ) ) ) ) ) ; h = colorbar ; ylabel( h, 'Log10(Intensity)', 'FontSize', 16 ) ;
        else
        imagesc( ( abs( squeeze( src_dff_0( i_rd, i_src_rf, i_src_plt, i_lmbd, :, : ) ) ) ) ) ; h = colorbar ; ylabel( h, 'Intensity (Linear)', 'FontSize', 16 ) ;
        end
      set( gca, 'xtick', [], 'ytick', [] )
      r_src_plt = sqrt( x_src( i_src_plt ) * x_src( i_src_plt ) + y_src( i_src_plt ) * y_src( i_src_plt ) ) ;
      title( sprintf( 'r_{src}=%3.2f mas, r_{ref}=%3.2f', r_src_plt, r_src_rf ) )
      end
    suptitle( sprintf( 'Reduction %s', lbl_rd{ i_rd } ) ) ;
    end
  end % i_rd
end % if ( 0 )

% Some plot of the DFT of the reference source
figure( round( opt_in.Nx_image_pix / 100 ) )
% For now looking at one lambda
i_lmbd = 1 ;
setwinsize( gcf, 1100, 800 )
subplot( 2, 2, 1 )
dft_tmp = abs( squeeze( circshift( dft_rf( 1, :, : ), [ 0, round( opt_in.Nx_image_pix / 2 ), round( opt_in.Nx_image_pix / 2 ) ]  ) )  ) ;
imagesc( dft_tmp ) ; h = colorbar ; ylabel( h, 'LINEAR' ) ;
grid
subplot( 2, 2, 2 )
imagesc( log10( dft_tmp ) ) ; h = colorbar ; ylabel( h, 'LOG10' ) ;
suptitle( sprintf( 'DFT OF REFERENCE SOURCE. PIXEL SIZE=%2.1f mas, LAMBDA=%3.1f nm,   X_{source}=%3.1f pix,    Y_{source}=%3.1f pix', opt_in.diam_img_mas / opt_in.Nx_img, 1e9 * lambdaIn( i_lmbd ), x_source_mas_rf / opt_in.diam_img_mas * opt_in.Nx_img, y_source_mas_rf / opt_in.diam_img_mas * opt_in.Nx_img ) ) ;
set( gca, 'FontSize', 16 ) ;
subplot( 2, 2, 3 : 4 )
h1 = plot( squeeze( dft_tmp( round( opt_in.Nx_img / 2 ), : ) ), 'x' ) ;
hold all
h2 = plot( squeeze( dft_tmp( :, round( opt_in.Nx_img / 2 ) ) ), 'x' ) ;
a1 = round( opt_in.Nx_img / 2 ) - 60 ;
a2 = round( opt_in.Nx_img / 2 ) + 60 ;
%h3 = plot( squeeze( dft_tmp( a1 : a2, a1 : a2 ) ), 'x' ) ;
%h4 = plot( squeeze( dft_tmp( a1 : a2, a2 : a1 ) ), 'x' ) ;

xlim( [ round( opt_in.Nx_img / 2 ) - 60, round( opt_in.Nx_img / 2 ) + 60 ] ) ;
ylabel( 'AMPLITUDE DFT', 'FontSize', 16 )
title( 'LONGITUDINAL CUTS ALONG X/Y AXIS AND ALONG DIAGONALS', 'FontSize', 16 )
grid
% Storing the image
  if ( opt_in.save_fig )
  img = getframe( gcf ) ;
  imwrite( img.cdata, sprintf( '%sDFT_rf_%i_alpha_%2.2f_pix_%2.2fmas%s', opt_in.save_path_fig, ...
           i_src_dft, alpha_lst( 1 ), opt_in.diam_img_mas / opt_in.Nx_img, '.png' ) ) ;
  end

dbstop if error
%make_an_error

%% Input resolution
a1 = sqrt( grp_lst( 1 ) ) ;
  if round( a1 ) == floor( a1 ), a1 = floor( a1 ) ; a2 = ceil( grp_lst( 1 ) / a1 ) ; else, a1 = ceil( a1 ) ; a2 = floor( grp_lst( 1 ) / a1 ) ; end
  for i_src_rf = 1 : 1 %n_src
    % For each subgroup
    idx_plt = 1 ;
    for i_grp = 2 : n_grp
    figure( i_grp + 1  )
    setwinsize( gcf, 1300, 750 )
    clf
    i_lmbd = 1 ;
    r_src_rf = sqrt( x_src( i_src_rf ) * x_src( i_src_rf ) + y_src( i_src_rf ) * y_src( i_src_rf ) ) ;
      for i_src_plt = idx_plt : idx_plt + grp_lst( i_grp ) - 1
      subplot( a1, a2, i_src_plt - idx_plt + 1 )
        if ( opt_in.log10 )
        imagesc( log10( abs( squeeze( src_dff( i_src_rf, i_src_plt, i_lmbd, :, : ) ) ) ) ) ; h = colorbar ; ylabel( h, 'Log10(Intensity)', 'FontSize', 16 ) ;
        else
        imagesc( ( abs( squeeze( src_dff( i_src_rf, i_src_plt, i_lmbd, :, : ) ) ) ) ) ; h = colorbar ; set( h, 'FontSize', 12 ) ; 
          if i_src_plt == idx_plt, ylabel( h, 'Intensity (Linear)', 'FontSize', 14 ) ; end
        end
      set( gca, 'xtick', [], 'ytick', [] ) 
      r_src_plt = sqrt( x_src( i_src_plt ) * x_src( i_src_plt ) + y_src( i_src_plt ) * y_src( i_src_plt ) ) ;
      %title( sprintf( 'r_{src}=%3.2f mas, r_{ref}=%3.2f', r_src_plt, r_src_rf ) )
      title( sprintf( 'r_{src}=%3.2f mas', r_src_plt ) )
      end
    suptitle( sprintf( 'Angle: %2.1f deg. Pixel size %2.2fx%2.2f mas', alpha_lst( i_grp ), 5 / fct_img_1, 5 / fct_img_1 ) )
    idx_plt = idx_plt + grp_lst( i_grp ) ;
    % Storing the image
      if ( opt_in.save_fig )
      img = getframe( gcf ) ;
      imwrite( img.cdata, sprintf( '%spanel_rf_%i_alpha_%2.2f_pix_%2.2fmas%s', opt_in.save_path_fig, i_src_rf, alpha_lst( i_grp ), 5 / fct_img_1, '.png' ) ) ;
      end
    end
  end

dbstop if error
%make_an_error

% Function that shifts according using DFT
function img_shft = dft_shift( img_in, dlt_x, dlt_y )

dft_in = fftshift( fft2( img_in ) ) ; 
n_in = size( img_in, 1 ) ;
% In case it's an odd number
n_in_hlf = floor( n_in / 2 ) ;
% Basic check
  if size( img_in, 2 ) ~= n_in
  disp( '(STARSHADE_TRANSLATION) Dimensions of the image are different. It ought to be a square. Stopping.' )
  dbstop if error
  make_an_error
  end

% Define shift in frequency domain
[ xF, yF ] = meshgrid( -n_in_hlf : n_in - n_in_hlf - 1, -n_in_hlf : n_in - n_in_hlf - 1 ) ;

% Perform the shift
dft_in = dft_in .* exp( -1i * 2 *pi .* ( xF * dlt_x + yF * dlt_y) / n_in ) ;

% Find the inverse Fourier Transform (imaginary part is zero or within numerical precison from zero)
img_shft = real( ifft2( ifftshift( dft_in ) ) ) ;

