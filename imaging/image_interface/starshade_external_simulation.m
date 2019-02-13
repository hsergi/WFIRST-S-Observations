function starshade_translation( opt )


% Reference source: it should be one. More if tests are conducted. The reference source comes from the analysis in starshade_translation
r_src = [ 400 ] ;
psi_src = [ 0 ] ;
n_src = numel( r_src ) ;

opt.developer = 1 ;
% Only two wavelengts
opt.delta_lambda_nm = 100 ;
opt.redo = 0 ;
opt.save = 1 ;
% Only one wavelength
opt = get_default_options( opt ) ;
  for i_src = 1 : n_src
  opt.r_source_mas = r_src( i_src ) ;
  opt.psi_source_deg = psi_src( i_src ) ;
  makeStarshadeImage( opt ) ;
  end

% Loading each simulation (notation from makeStarshadeImage)
dr_out = 'out_dev/' ;
% For files created before the new parameter was introduced
Nx_img_tmp = opt.Nx_img ;
diam_img_mas_tmp = opt.diam_img_mas ;
exozodi = opt.image_input.exozodi ;

  for i_src = 1 : n_src
  fl_out = sprintf( 'starshade_out_Nx_%i_pix_dl_%inm_dr_%3.1f_mas_psi_%3.1f_deg', opt.Nx_pupil_pix, opt.delta_lambda_nm, r_src( i_src ), psi_src( i_src ) ) ;
    if opt.Nx_img ~= 400, fl_out = sprintf( '%s_Nx_img_%04i', fl_out, opt.Nx_img ) ; end
    if opt.diam_img_mas ~= 2000, fl_out = sprintf( '%s_diam_%04i', fl_out, opt.diam_img_mas ) ; end
  load( [ fl_out '.mat' ] )

  opt.Nx_img = Nx_img_tmp ;
  opt.diam_img_mas = diam_img_mas_tmp ;
  opt.image_input.exozodi = exozodi ;

  n_lmbd = size( efDefectImg, 3 ) ;
    for i_lmbd = 1 : n_lmbd
    int_src( i_src, i_lmbd, :, : ) = abs( squeeze( efDefectImg( :, :, i_lmbd ) ) ).^2 ;
    end

  % Reading an input image to simulate
    if ( opt.image_input.exozodi == 1 )
    fl_img = 'heic0821c.jpg' ;
    end
    if ( opt.image_input.exozodi == 2 )
    fl_img = 'kalas_fomalhaut_wide.jpg' ;
    end

  img_in = imread( sprintf( 'in/%s', fl_img ) ) ;

  % Reference simulation for a single ray: the reference should have been computed for a ray exactly centered at a pixel
  rf_sm = circshift( squeeze( int_src( 1, 1, :, : ) ), ...
                     [ -r_src( i_src ) * cos( pi / 180 * psi_src( i_src ) ) * opt.Nx_img / opt.diam_img_mas, ...
                       -r_src( i_src ) * sin( pi / 180 * psi_src( i_src ) ) * opt.Nx_img / opt.diam_img_mas ] ) ;
  psf = fft2( rf_sm ) ;

  % Reducing the resolution of the reference simulation to the standard case
  rf_sm_rd = circshift( imresize( squeeze( int_src( 1, 1, :, : ) ), 400 / opt.Nx_img ), ...
                     [ -r_src( i_src ) * cos( pi / 180 * psi_src( i_src ) ), ...
                       -r_src( i_src ) * sin( pi / 180 * psi_src( i_src ) ) ] ) ;
  psf = fft2( rf_sm ) ;
  psf_rd = fft2( rf_sm_rd ) ;

  % Number of colors present
  n_rgb = size( img_in, 3 ) ;

    for i_rgb = 1 : n_rgb
    % choosing the color and re-sizing the image
    img_tmp = imresize( squeeze( img_in( :, :, i_rgb ) ), [ opt.Nx_img, opt.Nx_img ] ) ;
    img_tmp_rd = imresize( squeeze( img_in( :, :, i_rgb ) ), [ 400, 400 ] ) ;
    % testing the identity case: calibration of the circular shift of the Fourier Transform
      if i_rgb == 1
      img_tmp_id = double( img_tmp ) * 0. ;
      img_tmp_id( opt.Nx_img / 2, opt.Nx_img / 2 ) = 1. ;
      img_fft_id = fft2( img_tmp_id ) ;
      % The '+1' provides the best match between rf_sm and img_conv_id (1e-15 differences)
      img_conv_id = circshift( ifft2( img_fft_id .* psf ), [ opt.Nx_img / 2 + 1, opt.Nx_img / 2 + 1 ] ) ;
      img_tmp_id_rd = double( img_tmp_rd ) * 0. ;
      img_tmp_id_rd( 200, 200 ) = 1. ;
      img_fft_id_rd = fft2( img_tmp_id_rd ) ;
      % The '+1' provides the best match between rf_sm and img_conv_id (1e-15 differences)
      img_conv_id_rd = circshift( ifft2( img_fft_id_rd .* psf_rd ), [ 201, 201 ] ) ;
      end
    % The input image case:
    img_fft = fft2( img_tmp ) ;
    img_fft_rd = fft2( img_tmp_rd ) ;
    img_conv( :, :, i_rgb ) = circshift( ifft2( img_fft .* psf ), [ opt.Nx_img / 2 + 1, opt.Nx_img / 2 + 1 ] ) ;
    img_conv_rd( :, :, i_rgb ) = circshift( ifft2( img_fft_rd .* psf_rd ), [ 121, 201 ] ) ;
    end % i_rgb
  rgbImage = cat( 3, squeeze( img_conv( :, :, 1 ) ), squeeze( img_conv( :, :, 2 ) ), squeeze( img_conv( :, :, 3 ) ) ) ;
  rgbImage_rd = cat( 3, squeeze( img_conv_rd( :, :, 1 ) ), squeeze( img_conv_rd( :, :, 2 ) ), squeeze( img_conv_rd( :, :, 3 ) ) ) ;  
  disp( 'RGB image done' ) 
  disp( 'DFT spatially calibrated' )

  % Some display
  ttl_in_rgb = { 'INPUT R', 'INPUT G', 'INPUT B' } ;
  figure( 1 )
  setwinsize( gcf, 1300, 750 )
  clf  
    for i_rgb = 1 : 3
    subplot( 3, 3, i_rgb )
      if i_rgb == 1
      mn_img_in = min( min( squeeze( img_in( :, :, i_rgb ) ) ) ) ;
      mx_img_in = max( max( squeeze( img_in( :, :, i_rgb ) ) ) ) ;
      end
    imagesc( squeeze( img_in( :, :, i_rgb ) ), [ mn_img_in, mx_img_in / double( i_rgb ) ] ) 
    colorbar
    title( ttl_in_rgb{ i_rgb } )
    subplot( 3, 3, 3 + i_rgb )
      if i_rgb == 1
      mn_img = min( min( squeeze( img_conv( :, :, i_rgb ) ) ) ) ;
      mx_img = max( max( squeeze( img_conv( :, :, i_rgb ) ) ) ) ;
      end
    imagesc( squeeze( img_conv( :, :, i_rgb ) ) / mx_img * double( mx_img_in ), [ mn_img_in, mx_img_in / double( i_rgb ) ] )
    colorbar
    title( 'WFIRST 0.5 mas/pix' )
    subplot( 3, 3, 6 + i_rgb )
    imagesc( imresize( squeeze( img_conv( :, :, i_rgb ) ), 400 / opt.Nx_img ) / mx_img * double( mx_img_in ), [ mn_img_in, mx_img_in / double( i_rgb ) ] )
    colorbar
    title( 'WFIRST 5 mas/pix' )
    end % i_rgb
  figure( 2 )
  setwinsize( gcf, 1300, 750 )
  clf
    for i_rgb = 1 : 3
    subplot( 3, 3, i_rgb )
      if i_rgb == 1
      mn_img_in = min( min( squeeze( img_in( :, :, i_rgb ) ) ) ) ;
      mx_img_in = max( max( squeeze( img_in( :, :, i_rgb ) ) ) ) ;
      end
    imagesc( squeeze( img_in( :, :, i_rgb ) ), [ mn_img_in, mx_img_in / double( i_rgb ) ] )
    colorbar
    title( ttl_in_rgb{ i_rgb } )
    subplot( 3, 3, 3 + i_rgb )
      if i_rgb == 1
      mn_img_rd = min( min( squeeze( img_conv_rd( :, :, i_rgb ) ) ) ) ;
      mx_img_rd = max( max( squeeze( img_conv_rd( :, :, i_rgb ) ) ) ) ;
      end
    imagesc( squeeze( img_conv_rd( :, :, i_rgb ) ) / mx_img_rd * double( mx_img_in ), [ mn_img_in, mx_img_in / double( i_rgb ) ] )
    colorbar
    title( 'WFIRST 5 mas/pix (LowRes)' )
    subplot( 3, 3, 3 + i_rgb )
    imagesc( squeeze( img_conv_rd( :, :, i_rgb ) ) / mx_img_rd * double( mx_img_in ), [ mn_img_in, mx_img_in / double( i_rgb ) ] )
    colorbar
    title( 'WFIRST 5 mas/pix (LowRes)' )
    subplot( 3, 3, 6 + i_rgb )
    img_1 = squeeze( img_conv_rd( :, :, i_rgb ) ) ;
    img_2 = imresize( squeeze( img_conv( :, :, i_rgb ) ), 400 / opt.Nx_img ) ;
    imagesc( ( img_1 / max( max( img_1 ) ) - ...
             img_2 / max( max( img_2 ) ) ) * 100, [ -3, 3 ] )
    colorbar
    title( 'WFIRST 5 mas/pix (HighRes) - (LowRes) in %' )
    end % i_rgb

    
  suptitle( 'EXTRNAL IMAGE CONVOLVED WITH WFIRST PSF (NO BLOCKING EFFECT FROM STARSHADE)' ) 


  dbstop if error
  make_an_error
  end % i_src


