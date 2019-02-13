function starshade_radial_interpolation( opt_in, alph_rf_in )

% Function that tests the results of simulating exact PSF for WFIRST-S and compare them with translated, rotated or interpolated PSF
% History
% Suggestion from Stuart Shaklan (JPL).
% 07/09/17: first version. Sergi Hildebrandt (JPL/Caltech). 

% Pixel resolution (must be set in the call)
mas_pr_px = opt_in.diam_img_mas / opt_in.Nx_image_pix ;

% Angle step (degrees, the semi-angle of a WFIRST-S petal):
dlt_alph = 90 / 14 ; 

% Series of angles considered:
alph_lst = [ 0, 1 , 2, 3, 4, 5, 6, 7, 8, 14 ] * dlt_alph ;
n_alph = numel( alph_lst ) ;
% Angle of reference: element index of alph_lst
alph_rf = int8( alph_rf_in ) ;
% Array of angles not considered as reference
alph_lst_2 = alph_lst( find( alph_lst ~= alph_lst( alph_rf ) ) ) ;
n_alph_2 = numel( alph_lst_2 ) ;

% Building the simulations (or skipping them if they already exist)
r_src = ( opt_in.polar_radius_1 : opt_in.step_mas : opt_in.polar_radius_2 ) ;
  for i_alph = 1 : n_alph
  clear opt
  % Get default options
  opt_in.polar_alpha_1 = alph_lst( i_alph ) ;
  % Creating the array of simulations
  % radius 1->radius 2 for angle 1:
  x_src = r_src * cos( opt_in.polar_alpha_1 * pi / 180 ) ;
  y_src = r_src * sin( opt_in.polar_alpha_1 * pi / 180 ) ;
  n_src_arry = numel( x_src ) ;
  i_1 = 1 ;
  i_2 = n_src_arry ; % n_src (except for the example with 0.5 mas/pix which has too much memory overall)
    for i_src = i_1 : i_2
    opt_in.x_source_mas = x_src( i_src ) ;
    opt_in.y_source_mas = y_src( i_src ) ;
    makeStarshadeImage( opt_in ) ;
    end
  end % i_alph

% Loading each simulation (notation from makeStarshadeImage)
dr_out = 'out_dev/' ;
% Loading the reference case
opt_in.polar_alpha_1 = alph_lst( alph_rf ) ;
x_src_rf = ( opt_in.polar_radius_1 : opt_in.step_mas : opt_in.polar_radius_2 ) * cos( opt_in.polar_alpha_1 * pi / 180 ) ;
y_src_rf = ( opt_in.polar_radius_1 : opt_in.step_mas : opt_in.polar_radius_2 ) * sin( opt_in.polar_alpha_1 * pi / 180 ) ;
  for i_src = i_1 : i_2
  % Creating the filename
  opt_in.x_source_mas = x_src_rf( i_src ) ;
  opt_in.y_source_mas = y_src_rf( i_src ) ;
  opt_in = get_default_options( opt_in ) ;
  load( [  opt_in.save_filename '.mat' ] )
  % Generating the image field
  efDefectImg_rf = UTot2ImagePlane( lambdaIn, opt_in, pupil, telescopeDiameter, UTotL ) ;
  efDefectImg_rf_arry( i_src, :, :, : ) = efDefectImg_rf ;
  n_lmbd = numel( lambdaIn ) ;
    for i_lmbd = 1 : n_lmbd
    int_src_rf( i_src, i_lmbd, :, : ) = abs( squeeze( efDefectImg_rf( :, :, i_lmbd ) ) ).^2 ;
    end % i_lmbd
  end % i_src

  for i_alph = 1 : n_alph_2
  opt_in.polar_alpha_1 = alph_lst_2( i_alph ) ;
  x_src = ( opt_in.polar_radius_1 : opt_in.step_mas : opt_in.polar_radius_2 ) * cos( opt_in.polar_alpha_1 * pi / 180 ) ;
  y_src = ( opt_in.polar_radius_1 : opt_in.step_mas : opt_in.polar_radius_2 ) * sin( opt_in.polar_alpha_1 * pi / 180 ) ;
    for i_src = i_1 : i_2
    % Creating the filename
    opt_in.x_source_mas = x_src( i_src ) ;
    opt_in.y_source_mas = y_src( i_src ) ;
    opt_in = get_default_options( opt_in ) ;
    load( [  opt_in.save_filename '.mat' ] )
    % Generating the image field
    efDefectImg = UTot2ImagePlane( lambdaIn, opt_in, pupil, telescopeDiameter, UTotL ) ;
    % For plotting purposes
      if i_src == alph_rf
      x_source_mas_rf = opt_in.x_source_mas ;
      y_source_mas_rf = opt_in.y_source_mas ;
      end
    n_lmbd = numel( lambdaIn ) ;
      for i_lmbd = 1 : n_lmbd
      int_src( i_src, i_lmbd, :, : ) = abs( squeeze( efDefectImg( :, :, i_lmbd ) ) ).^2 ;
      % Comparisons to be investigated for each source
      %% 1) Cartesian translation (worst method)
      dlt_x = ( x_src( i_src ) - x_src_rf( i_src ) ) / mas_pr_px ;
      dlt_y = ( y_src( i_src ) - y_src_rf( i_src ) ) / mas_pr_px ;
      img_rf = squeeze( int_src_rf( i_src, i_lmbd, :, : ) ) ;
      img_shft = dft_shift( img_rf, dlt_y, dlt_x ) ;
      src_tmp = squeeze( int_src( i_src, i_lmbd, :, : ) ) ;
      dff_img = src_tmp - img_shft ;
      % Region around the peak source
      a_1 = round( ( y_src( i_src ) - 150 ) / mas_pr_px + opt_in.Nx_image_pix / 2 ) ;
      a_2 = round( ( y_src( i_src ) + 150 ) / mas_pr_px + opt_in.Nx_image_pix / 2 ) ;
      b_1 = round( ( x_src( i_src ) - 150 ) / mas_pr_px + opt_in.Nx_image_pix / 2 ) ;
      b_2 = round( ( x_src( i_src ) + 150 ) / mas_pr_px + opt_in.Nx_image_pix / 2 ) ;
      % Relative to the peak
      dff_img_tmp =  dff_img( b_1 : b_2, a_1 : a_2 ) / max( max( abs( src_tmp( b_1 : b_2, a_1 : a_2 ) ) ) ) ;
      % Storing the max value of the relative difference for each method
      dff_img_mx( 1, i_lmbd, i_alph, i_src - i_1 + 1 ) = max( abs( dff_img_tmp( : ) ) ) ;
      % Image of the image field for each method
        if ( opt_in.plot ), whitebg( 'k' ) ; end % if removed, change the color in the summary plot at the end for the line with '--' to black, instead of white
        if ( opt_in.plot == 2 )
        close
        figure( 1 )
        setwinsize( gcf, 900, 230 )
        subplot(131)
        imagesc( dff_img_tmp ) ; colorbar ; title( 'EXACT SIM - SHIFTED REFERENCE', 'FontSize', 12 )
        end
      %% Rotation in the image field
      img_rf_rtd = imrotate( img_rf, alph_lst_2( i_alph ) - alph_lst( alph_rf ), opt_in.imrotate, 'Crop' ) ;
      dff_img = src_tmp - img_rf_rtd ;      
      % Relative to the peak
      dff_img_tmp =  dff_img( b_1 : b_2, a_1 : a_2 ) / max( max( abs( src_tmp( b_1 : b_2, a_1 : a_2 ) ) ) ) ;
      % Storing the max value of the relative difference for each method
      dff_img_mx( 2, i_lmbd, i_alph, i_src - i_1 + 1 ) = max( abs( dff_img_tmp( : ) ) ) ;
        if ( opt_in.plot == 2 )
        subplot(132)
        imagesc( dff_img_tmp ) ; colorbar ; title( 'EXACT SIM - ROTATED REFERENCE', 'FontSize', 12 )
        end

      %% Rotation of UTotL

        if ( opt_in.plot == 2 ) 
        suptitle( sprintf( 'RELATIVE DIFFERENCES WITH RESPECT TO PEAK VALUE OF EXACT SIMULATION (R=%3.0f mas, ALPHA=%2.1f deg)', r_src( i_src ), alph_lst_2( i_alph ) ) ) ;
        end
      % Storing the image
        if ( opt_in.save_fig ) && ( opt_in.plot == 2 )
          if ~isdir( opt_in.save_path_fig ), system( [ 'mkdir -p ' opt_in.save_path_fig ] ) ; end
        img = getframe( gcf ) ;
        imwrite( img.cdata, sprintf( '%sradial_%s_%2.2fmas_rf_%i.%s', opt_in.save_path_fig, ...
                          opt_in.save_filename, opt_in.diam_img_mas / opt_in.Nx_img, alph_rf, 'png' ) ) ;
        end
      end % i_lmbd
    end % i_src
  end % i_alph

  % Plot with the overall results 
  for i_alph = 1 : n_alph_2
  alph_lst_str{ i_alph } = sprintf( '%2.1f deg', alph_lst_2( i_alph ) ) ;
  end

  if ( opt_in.plot )
  ttl_mthd = { 'SHIFTING REFERENCE IMAGE', 'ROTATING REFERENCE IMAGE', 'ROTATING REFERENCE UTOT' } ;
    for i_lmbd = 1 : n_lmbd
    close ; clf
    figure( 1 )
    setwinsize( gcf, 1200, 550 )
    % For each method
    n_mthd = 2 ;
      for i_mthd = 1 : 2 % n_mthd
      clear h
      subplot( sprintf( '1%i%i', n_mthd, i_mthd ) )
      hold all
        for i_alph = 1 : n_alph_2
        rgb_tmp = get_rgb_colors( i_alph ) ;
        h( i_alph ) = plot( r_src( i_1 : i_2 ), 100 * squeeze( dff_img_mx( i_mthd, i_lmbd, i_alph, : ) ), 'LineWidth', 1.5, 'Color', rgb_tmp ) ;
        end  
      % Target level of precision
      plot( r_src( i_1 : i_2 ), 1 * ones( 1, i_2 - i_1 + 1 ), 'w--', 'LineWidth', 1.5 ) ; 
      ylabel( 'Relative Difference wrt peak (%)', 'Fontsize', 12 )
      xlabel( 'Radial distance to center of Starshade (mas)', 'Fontsize', 12 )
      xlim( [ r_src( 1 ), r_src( end ) ] ) ;
      legend( h, alph_lst_str ) ;
      title( ttl_mthd{ i_mthd } ) ;
      box on
      end % i_mthd
    legend( h, alph_lst_str ) ;
    suptitle( sprintf( 'Maximum Relative Difference wrt peak in %s. Lambda=%3.2f nm. 1 PIX=%2.1f mas. ROTATION METHOD: %s', '%', ...
                       lambdaIn * 1e9, mas_pr_px, opt_in.imrotate ) ) ; 
    
     if ( opt_in.save_fig )
       if ~isdir( opt_in.save_path_fig ), system( [ 'mkdir -p ' opt_in.save_path_fig ] ) ; end
       img = getframe( gcf ) ;
       imwrite( img.cdata, sprintf( '%sradial_summary_%s_rf_%i_%3.1fnm_%2.2fmas%s', opt_in.save_path_fig, ...
                opt_in.imrotate, alph_rf, 1e9 * lambdaIn( i_lmbd ), opt_in.diam_img_mas / opt_in.Nx_img, '.png' ) ) ;
     end    
    end % i_lmbd
  end % plot

dbstop if error
make_an_error 




