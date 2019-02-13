function starshade_interpolation( opt_in )
% Function that tests the results of imaging starshade exactly and with spatical/wavelength interpolation methods
% History
% 05/14/17: first version. Sergi Hildebrandt (JPL/Caltech)
% 05/22/17: added the analysis of the simulations to evaluate the goodness of interpolation

% Get default options
opt = get_default_options( opt_in ) ;

% 1) Create the reference simulations or study interpolation
  if ( opt.create_reference )
  opt.save = 1 ;

  % Star simulation only (on- or close to on-axis)
    if ( opt.star )
    % Studying some relative positions
    %% Small displacements for the star (pointing errors) centered at the center
    opt.x_source_mas_array = [ -20 : 1 : 20 ] ; % mas
    opt.y_source_mas_array = [ -20 : 1 : 20 ] ; % mas
    end

    %% Medium displacements for the planet centered at some position
    %%% Number of exact simulations to be performed at each side of the central one
    if ( opt.planet )
    n_basis_l = ceil( ( opt.n_basis_interpolation - 1 ) / 2 ) ;
    opt.x_source_mas_array = [ opt.x_planet_mas - opt.step_mas * n_basis_l : opt.step_mas : opt.x_planet_mas + opt.step_mas * n_basis_l ] ;
    opt.y_source_mas_array = [ opt.y_planet_mas - opt.step_mas * n_basis_l : opt.step_mas : opt.y_planet_mas + opt.step_mas * n_basis_l ] ;
    end

    % Stuart suggestion about considering radial basis for the interpolation
    if ( opt.polar )
    % Keeping this for now
    n_basis_l = ceil( ( opt.n_basis_interpolation - 1 ) / 2 ) ;
    % Creating the array of simulations
    % radius 1->radius 2 for angle 1:
    opt.x_source_mas_array = ( opt.polar_radius_1 : opt.step_mas : opt.polar_radius_2 ) * cos( opt.polar_alpha_1 * pi / 180 ) ;
    opt.y_source_mas_array = ( opt.polar_radius_1 : opt.step_mas : opt.polar_radius_2 ) * sin( opt.polar_alpha_1 * pi / 180 ) ;
    % radius 1 to radius 2 for angle 2
    opt.x_source_mas_array = [ opt.x_source_mas_array, ( opt.polar_radius_1 : opt.step_mas : opt.polar_radius_2 ) * cos( opt.polar_alpha_2 * pi / 180 ) ] ;
    opt.y_source_mas_array = [ opt.y_source_mas_array, ( opt.polar_radius_1 : opt.step_mas : opt.polar_radius_2 ) * sin( opt.polar_alpha_2 * pi / 180 ) ] ;
    end

  n_sm_x = numel( opt.x_source_mas_array ) ;
  n_sm_y = numel( opt.y_source_mas_array ) ;
  % For the non-polar simulations
    if ~( opt.polar )
    disp( sprintf( '(starshade_interpolation) Considering %i reference simulations', n_sm_x * n_sm_y ) )
    opt_sm = opt ;
      for i_sm_x = 1 : n_sm_x   
      tic
      disp( sprintf( '(starshade_interpolation) ****** Considering the simulation %i/%i ******', i_sm_x * n_sm_y, n_sm_x * n_sm_y ) )
      opt_sm.x_source_mas = opt.x_source_mas_array( i_sm_x ) ;
        for i_sm_y = 1 : n_sm_y
        opt_sm.y_source_mas = opt.y_source_mas_array( i_sm_y ) ;
        makeStarshadeImage( opt_sm ) ;
        end % i_sm_y
      disp( 'Last block of simulations took: ' ) 
      toc
      end % i_sm_x
    end % For the non-polar case

    % For the polar case
    if ( opt.polar )
    disp( sprintf( '(starshade_interpolation) Considering %i reference simulations', n_sm_x ) )
    opt_sm = opt ;
      for i_sm_x = 1 : n_sm_x
      tic
      disp( sprintf( '(starshade_interpolation) ****** Considering the simulation %i/%i ******', i_sm_x, n_sm_x ) )
      opt_sm.x_source_mas = opt.x_source_mas_array( i_sm_x ) ;
      opt_sm.y_source_mas = opt.y_source_mas_array( i_sm_x ) ;
      makeStarshadeImage( opt_sm ) ;
      disp( 'Last block of simulations took: ' )
      toc
      end % i_sm_x
    end

  end % Create the reference simulations
dbstop if error
disp( '(starshade_interpolation) Simulations created' ) ;
%make_an_error
return
% 2) Analyze the results
%% 2.1) Exploration (non-optimal in terms of speed and memory, but most accurate)
  if ( 1 ) && ~( opt.polar )
    % Collecting all the data
  tic
    for i_sm_x = 1 : n_sm_x 
    opt_sm = opt ;
    opt_sm.x_source_mas = opt.x_source_mas_array( i_sm_x ) ;
      for i_sm_y = 1 : n_sm_y
      opt_sm.y_source_mas = opt.y_source_mas_array( i_sm_y ) ;
      opt_sm = get_default_options( opt_sm ) ;
      % Same naming as in makeStarshadeImage
      %% Main parameters
      Nx = opt_sm.Nx_pupil_pix ;
      dlt_lmbd = opt_sm.delta_lambda_nm ;
      r_src = opt_sm.r_source_mas ;
      psi_src = opt_sm.psi_source_deg ;
      ppl_fl = opt_sm.pupil_file ;
      savePath = opt_sm.save_path ;
      saveFilename = sprintf( 'starshade_out_Nx_%i_pix_dl_%inm_dr_%3.1f_mas_psi_%3.1f_deg', Nx, dlt_lmbd, r_src, psi_src ) 
        if strcmp( ppl_fl, '0' ) == 1, saveFilename = sprintf( '%s_ideal', saveFilename ) ; end
      load( [ savePath saveFilename '.mat' ] )
      % Making sure the options are the ones in the call and not the ones stored:
      opt.super_resolution.res = opt_in.super_resolution.res  ;
      opt.super_resolution.interp_method = opt_in.super_resolution.interp_method ;
      opt.low_resolution.res = opt_in.low_resolution.res ;
      opt.low_resolution.interp_method = opt_in.low_resolution.interp_method ;
      opt.step_mas = opt_in.step_mas ;
        if ( opt.star )
        % Studying some relative positions
        %% Small displacements for the star (pointing errors) centered at the center
        opt.x_source_mas_array = [ -5 : 1 : 5 ] ; % mas
        opt.y_source_mas_array = [ -5 : 1 : 5 ] ; % mas
        end
        %% Medium displacements for the planet centered at some position
        %%% Number of exact simulations to be performed at each side of the central one
        if ( opt.planet )
        n_basis_l = ceil( ( opt_in.n_basis_interpolation - 1 ) / 2 ) ;
        opt.x_source_mas_array = [ opt.x_planet_mas - opt.step_mas * n_basis_l : opt.step_mas : opt.x_planet_mas + opt.step_mas * n_basis_l ] ;
        opt.y_source_mas_array = [ opt.y_planet_mas - opt.step_mas * n_basis_l : opt.step_mas : opt.y_planet_mas + opt.step_mas * n_basis_l ] ;
        end
      % Working with the intensity only
      tmp = abs( efDefectImg( :, :, : ) ).^2 ;
      %% Number of wavelengths considered
      n_lmbd = numel( lambdaIn ) ;
        % First the grid needs some increased resolution
        %% For each wavelength
        for i_lmbd = 1 : n_lmbd
        tmp_spr_rsl( :, : ) = interp2( squeeze( tmp( :, :, i_lmbd ) ), opt.super_resolution.res, opt.super_resolution.interp_method ) ;
        IntDefectImg( :, :, i_lmbd, i_sm_y, i_sm_x ) = tmp_spr_rsl( :, : ) ;
        end
      end % i_sm_y
    end % i_sm_x
 % Size of the super-resolution simulations
 n_2 = size( IntDefectImg, 1 ) ;
 n_1 = size( IntDefectImg, 2 ) ;
 % Mesh where final interpolation will be performed
 [ xq, yq ] = meshgrid( ( 1 : n_sm_x ), ( 1 : n_sm_y ) ) ;
 disp( sprintf( '(starshade_interpolation) Time spent preparing the interpolation data %2.1f seconds', toc ) )
   % Estimating the error performed doing an interpolation
   %% For each wavelength
     for i_lmbd = 1 : n_lmbd
     % for the non-polar case
       if ~( opt.polar )
       % Selecting a subset of simulations as the interpolation basis
       x_basis_arry = ( 1 : opt.low_resolution.res : n_sm_x ) ;
       n_x_basis = numel( x_basis_arry ) ;
       y_basis_arry = ( 1 : opt.low_resolution.res : n_sm_y ) ;
       n_y_basis = numel( y_basis_arry ) ;
       % Associated mesh
       [ x_basis, y_basis ] = meshgrid( x_basis_arry, y_basis_arry ) ;
       % Simulation basis transformed into 1-dimensional arrays in the data
       sim_basis = reshape( squeeze( IntDefectImg( :, :, i_lmbd, y_basis_arry, x_basis_arry ) ), n_2 * n_1, n_y_basis, n_x_basis ) ;
       % Array for the final results
       intrp_sim = NaN( n_1 * n_2, n_sm_y, n_sm_x ) ;
       % This is very slow
       % Going case by case, to center the interpolation operations around the region of importance
       tic
       % Limits chosen for the simulation, plus some times the FWHM 
       l1_d = opt.x_planet_mas - opt.step_mas * n_basis_l - 20 ;
       l1_d = l1_d * opt.super_resolution.res + ceil( n_1 / 2 ) ;
       l1_u = opt.x_planet_mas + opt.step_mas * n_basis_l + 20;
       l1_u = l1_u * opt.super_resolution.res + ceil( n_1 / 2 ) ;
       l2_d = opt.y_planet_mas - opt.step_mas * n_basis_l - 20 ;
       l2_d = l2_d * opt.super_resolution.res + ceil( n_2 / 2 ) ;
       l2_u = opt.y_planet_mas + opt.step_mas * n_basis_l + 20 ;
       l2_u = l2_u * opt.super_resolution.res + ceil( n_2 / 2 ) ;
         for i_2 = l2_d : l2_u
           for i_1 = l1_d : l1_u
           intrp_sim( n_1 * i_1 + i_2 - n_1, :, : ) = interp2( x_basis, y_basis, ...
           squeeze( sim_basis( n_1 * i_1 + i_2 - n_1, :, : ) ), xq, yq, opt.low_resolution.interp_method ) ;
           end
         end
       disp( sprintf( '(starshade_interpolation) Time spent interpolating data %2.1f seconds (%i/%i)', toc, i_lmbd, n_lmbd ) )
       % Re-shaping the result
       intrp_sim_2d( :, :, i_lmbd, 1 : n_sm_y, 1 : n_sm_x ) = reshape( intrp_sim, n_2, n_1, n_sm_y, n_sm_x ) ;
       % reshaping the basis
       sim_basis_2d( :, :, i_lmbd, y_basis_arry, x_basis_arry ) = reshape( sim_basis, n_2, n_1, n_y_basis, n_x_basis ) ;
       end % If it is not the polar case
     end % i_lmbd
  end % Exploration method

% Some figures
% Reduce to original grid size
% Add Fourier
% loop over cases high_res/low_res and lambda (each one a figure)
  if ( 1 ) & ~( opt.polar )
  % Wavelength
  i_lmbd = 1 ;
  % Interpolation basis (odd double indices)
  i_1 = 1 ; i_2 = 1 ;
  % Predictive result (any other case)
  i_1 = 9 ; i_2 = 9 ;
 
  close all
  setwinsize(gcf,800,650)
  tmp = abs( squeeze( intrp_sim_2d( :, :, i_lmbd, i_1, i_2 ) ) ) ;
  q_mx = find( tmp == max( tmp( : ) ) ) ;
  a1 =  round( q_mx / size( tmp, 1 ) ) ;
  a2 = mod( q_mx, size( tmp, 1 ) ) ;
  x_arry = a1 - 40 : a1 + 40 ;
  y_arry = a2 - 40 : a2 + 40 ;
  tmp_1 = squeeze( IntDefectImg( x_arry, y_arry, i_lmbd, i_1, i_2 ) ) ;
  tmp_2 = squeeze( intrp_sim_2d( x_arry, y_arry, i_lmbd, i_1, i_2 ) ) ;
  subplot( 2, 2, 1 )
  imagesc( x_arry, y_arry, tmp_1 ) ; title( 'Exact simulation', 'FontSize', 14 ) ; h = colorbar ; set( h, 'FontSize', 14 ) ; grid
  xlabel( 'Pixel' ) ; ylabel( 'Pixel' ) ;
  subplot( 2, 2, 2 )
  imagesc( x_arry, y_arry, tmp_2 ) ; title( 'Interpolated Simulation', 'FontSize', 14 ) ; h = colorbar ; set( h, 'FontSize', 14 ) ; grid
  xlabel( 'Pixel' ) ; ylabel( 'Pixel' ) ;
  subplot( 2, 2, 3 )
  % Plotting the product of the input signal and the difference
  plt_dt = ( tmp_1 - tmp_2 ) ;
  mn_1 = min( plt_dt( : ) )  ;
  mx_1 = max( plt_dt( : ) )  ;
    if mn_1 == mx_1, mn_1 = mn_1 - 1e-10 ; mx_1 = mx_1 + 1e-10 ; end
  imagesc( x_arry, y_arry, plt_dt, [ mn_1, mx_1 ] ) ; 
  title( 'Difference exact minus interpolated', 'FontSize', 14 ) ; h = colorbar ; set( h, 'FontSize', 14 ) ; grid
  xlabel( 'Pixel' ) ; ylabel( 'Pixel' ) ;
  subplot( 2, 2, 4 )
  % Plotting the product of the input signal and the difference
  plt_dt = tmp_1 .* ( tmp_1 - tmp_2 ) ; %tmp_1 .* ( tmp_1 - tmp_2 / max( tmp_2( : ) ) * max( tmp_1( : ) ) ) ;
  imagesc( x_arry, y_arry, plt_dt, [ mn_1, mx_1 ] ) ; 
  title( '(Difference exact minus interpolated)xexact', 'FontSize', 14 ) ; h = colorbar ; set( h, 'FontSize', 14 ) ; grid ;
  xlabel( 'Pixel' ) ; ylabel( 'Pixel' ) ;
  suptitle( sprintf( 'Occulter: %s, lambda=%3.0f nm (I1=%i/%i=%2d mas, I2=%i/%i=%2d mas)', opt.occulter_name, lambdaIn( i_lmbd ) * 1e9, i_1, numel( opt.x_source_mas_array ), ...
            opt.x_source_mas_array( i_1 ), i_2, numel( opt.y_source_mas_array ), opt.y_source_mas_array( i_2 ) ) ) 

dbstop if error
make_an_error


    if ( 0 )
    setwinsize(gcf,900,350)
    in = squeeze( IntDefectImg( a1 : a2, a1 : a2, i_lmbd, 1, 1 ) ) ; ;
    H=fftshift(fft2(in)); %// Compute 2D Fourier Transform
    x0= 5; %// Define shifts
    y0= 5;

    %// Define shift in frequency domain
    [xF,yF] = meshgrid(-400:398,-400:398);

    %// Perform the shift
    H=H.*exp(-1i*2*pi.*(xF*x0+yF*y0)/799);

    %// Find the inverse Fourier Transform
    IF_image=ifft2(ifftshift(H));

    %// Show the images
    figure;
    subplot(1,3,1);
    imagesc(in); grid
    subplot(1,3,2);
    imagesc(real(IF_image)); grid
    subplot(1,3,3);
    in_2 = squeeze( IntDefectImg( a1 : a2, a1 : a2, i_lmbd, 2, 2) ) ;
    imagesc( in_2 ) ; grid
    end

    % Plotting the interpolation basis
    for i_lmbd = 1 : 1 %n_lmbd 
      setwinsize(gcf,1300,750)
      clf
      hold all
      mx_img = 1 ;
      mn_img = 1e-13 ;
      dlt = 120 ;
      a1 = 400 + 40 - dlt ;
      a2 = 400 + 40 + dlt ;
      for iy = 1 : numel( y_basis_arry )
        for ix = 1 : numel( x_basis_arry )
        subplot( numel( y_basis_arry ), numel( x_basis_arry ), ( iy - 1 ) * numel( x_basis_arry ) + ix )
        imagesc( a1 : a2, a1 : a2, log10( squeeze( sim_basis_2d( a1 : a2, a1 : a2, i_lmbd, y_basis_arry( iy ), x_basis_arry( ix ) ) ) ), [ log10( mn_img ), log10( mx_img ) ] ) ; colorbar
%        imagesc( a1 : a2, a1 : a2, squeeze( sim_basis_2d( a1 : a2, a1 : a2, i_lmbd, y_basis_arry( iy ), x_basis_arry( ix ) ) ), [ mn_img, mx_img ] ) ; colorbar
        title( sprintf( '(%i,%i)', y_basis_arry( iy ), x_basis_arry( ix ) ) )
        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])
        end
      end
      suptitle( { sprintf( 'INTERPOLATION BASIS LOG10 SCALE. Occulter: %s, lambda=%3.0f nm', opt.occulter_name, lambdaIn( i_lmbd ) * 1e9 ) } )
%      suptitle( { sprintf( 'INTERPOLATION BASIS LINEAR SCALE. Occulter: %s, lambda=%3.0f nm', opt.occulter_name, lambdaIn( i_lmbd ) * 1e9 ) } )
    end % i_lmbd

  end % ~opt.polar
% In order to study the interpolation limitations

simula_px = squeeze( IntDefectImg( 399, 399, 1, :, : ) ) ;
sim_basis_px = squeeze( sim_basis_2d( 399, 399, 1, :, : ) ) ;
interp_px = squeeze( intrp_sim_2d( 399, 399, 1, :, : ) ) ;

clf
subplot(311)
imagesc( log10( simula_px ) ) ; colorbar ; title( 'SIMULATIONS' )
subplot(312)
imagesc( log10( sim_basis_px ) ) ; colorbar ; title( 'SIMULATION BASIS' )
subplot(313)
imagesc( log10( abs( interp_px ) ) ) ; colorbar ; title( 'INTERPOLATIONS' )


% 1-D
clf
idx = 9 ;
plot( simula_px( idx, : ), '+-' )
hold all
plot( sim_basis_px( idx, : ), 'x-' )
plot( interp_px( idx, : ), '+-' )

