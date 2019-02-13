function opt = get_default_options( opt )
% Function to define the default options for Starshade imaging
% History:
% 05/12/17: created, Sergi Hildebrandt (JPL/Caltech)

% Developer version or Eric Cady's original version
  if ~isfield( opt, 'developer' )
  opt.developer = 0 ;
  end

% Switch to control some work that needs to be re-done or may be skipped if it already exists
  if ~isfield( opt, 'redo' )
  opt.redo = 0 ;
  end

% Occulter path
  if ~isfield( opt, 'path_occulter' )
  opt.path_occulter = './' ;
  end

% Occulter Name
  if ~isfield( opt, 'occulter_name' )
  opt.occulter_name = 'NI2' ;
  end

% Number of petals
  if ~isfield( opt, 'n_ptl' )
  opt.n_ptl = 24 ; % Default number (S. Shaklan email "NI2 24 petals files" 05/03/18)
  end

% Size of the pupil data in pixels (square)
  if ~isfield( opt, 'Nx_pupil_pix' )
  opt.Nx_pupil_pix = 64 ;
  end

% Size of the secondary if a circular obscuration is fine (no struts). See makeStarshadeImage.m
  if ~isfield( opt, 'secondary_size' )
  opt.secondary_size = 0 ;
  end

% Number of pixels across image plane
  if ~isfield( opt, 'Nx_image_pix' )
  opt.Nx_img = 400 ;
  else
  opt.Nx_img = opt.Nx_image_pix ;
  end

% Diameter of image plane in milliarcseconds
  if ~isfield( opt, 'diam_image_mas' )
  opt.diam_img_mas = 2000 ;
  else
  opt.diam_img_mas = opt.diam_image_mas ;
  end

% Initial wavelength to consider
  if ~isfield( opt, 'lambda_1_nm' )
  opt.lambda_1_nm = 425 ; 
  end

% Final wavelength to consider
  if ~isfield( opt, 'lambda_2_nm' )
  opt.lambda_2_nm = 525 ;
  end

% Step of wavelength to consider
  if ~isfield( opt, 'delta_lambda_nm' )
  opt.delta_lambda_nm = 10 ; % nm
  end

% X/Y positions on the focal plane (instead of r_source_mas and psi_source_deg)
% Matlab: column major. Therefore in an image, X/Y are to be understood as -y/x in usual Cartesain convention (e.g., X=1, Y=0 means (0,-1) in sual x/y axis)
% Convention ox X/Y with respect to r/Psi as in bdwf core function
% s1 = sin(psi1);
% c1 = cos(psi1);
% s2 = sin(psi2);
% c2 = cos(psi2);

% Checks of consistency
  if isfield( opt, 'r_source_mas' ) && ~isfield( opt, 'psi_source_deg' )
  disp( '(get_default_options) r_source_mas set, but psi_source_deg not. Inconsistent. Returning.' )
  return
  end
  if isfield( opt, 'psi_source_deg' ) && ~isfield( opt, 'r_source_mas' )
  disp( '(get_default_options) psi_source_deg set, but r_source_mas not. Inconsistent. Returning.' )
  return
  end

% Separation of the source from the center of the pointing
  if ~isfield( opt, 'r_source_mas' )
  opt.r_source_mas = 0 ; % mas
  end

% Angle of the source with respect the horizontal axis
  if ~isfield( opt, 'psi_source_deg' )
  opt.psi_source_deg = 0 ; % degrees
  end

% Consistency check
  if isfield( opt, 'y_source_mas' ) && ~isfield( opt, 'x_source_mas' )
  disp( '(get_default_options) If opt.y_source_mas is set, opt.x_source_mas must also be set. Returning.' )
  return
  end
  if isfield( opt, 'x_source_mas' ) && ~isfield( opt, 'y_source_mas' )
  disp( '(get_default_options) If opt.x_source_mas is set, opt.y_source_mas must also be set. Returning.' )
  return
  end

% Sentinel for changes in the filename where results are stored
rplc_fl_nm = 0 ;

  if isfield( opt, 'x_source_mas' )
  % Consistency check
    if ~isfield( opt, 'y_source_mas' )
    disp( '(get_default_options) If opt.x_source_mas is set, opt.y_source_mas must also be set. Returning.' ) 
    return
    end
  % Transform into r and psi
  % Different situations where it is requeired to check the consistency between r_source_mas, psi_source_deg, x_source_mas and y_source_mas:
    if ( opt.x_source_mas ) || ( opt.y_source_mas ) 
      if ( ( sqrt( opt.x_source_mas^2 + opt.y_source_mas^2 ) ~= opt.r_source_mas ) || ( atan2( opt.y_source_mas, opt.x_source_mas ) * 180 / pi ~= opt.psi_source_deg ) ) 
      r_source_mas_tmp = sqrt( opt.x_source_mas^2 + opt.y_source_mas^2 ) ;
      psi_source_deg_tmp = atan2( opt.y_source_mas, opt.x_source_mas ) * 180 / pi ; % deg
      disp( sprintf( '(get_default_options) Changing the values of r_source_mas and psi_source_deg from %3.3f, %3.3f to %3.3f, %3.3f', opt.r_source_mas, opt.psi_source_deg, r_source_mas_tmp, psi_source_deg_tmp ) )
      opt.r_source_mas = r_source_mas_tmp ;
      opt.psi_source_deg = psi_source_deg_tmp ;
      % Update filename where results are stored
      rplc_fl_nm = 1 ;
      end
    end
  end
  
  %% If they are not fields, create them
    if ~isfield( opt, 'r_source_mas' ) 
    opt.r_source_mas = sqrt( opt.x_source_mas^2 + opt.y_source_mas^2 ) ;
    end
    %%% Recall matlab column major convention and atan2(y,x)
    if ~isfield( opt, 'psi_source_deg' )
      if opt.r_source_mas % if the radius is zero, the angle does not matter.
      opt.psi_source_deg = atan2( opt.y_source_mas, opt.x_source_mas ) * 180 / pi ; % deg
      else
      opt.psi_source_deg = 0 ;
      end
    end

% For ease identification of the pixel on the image
opt.x_source_mas = opt.r_source_mas * cos( opt.psi_source_deg * pi / 180 ) ;
opt.y_source_mas = opt.r_source_mas * sin( opt.psi_source_deg * pi / 180 ) ;

% Hot pixels in the electric fields (by default, do not erase them)
  if ~isfield( opt, 'erase_hot_pixels' )
  opt.erase_hot_pixels = 0 ;
  end

% Replace by a file in FITS format with an Nx x Nx array if you want a specific pupil
  if ~isfield( opt, 'pupil_file' )
  opt.pupil_file = [ opt.path_occulter 'pupil_D1Kpix_2048.fits' ] ;
  end

% NB: Name of the filename to store the results at the end of the code
% Saving all the output results and images (0=No, 1=Yes)
  if ~isfield( opt, 'save_all' )
  opt.save_all = 0 ;
  end

% Saving only output results
  if ~isfield( opt, 'save' )
  opt.save = 0 ;
  end
% Saving figures
  if ~isfield( opt, 'save_fig' )
  opt.save_fig = 0 ;
  end

% If all saved, then:
  if opt.save_all
  opt.save = 1 ;
  opt.save_fig = 1 ;
  end

% paths to save the results
  if ~isfield( opt, 'save_path' )
  opt.save_path = './out' ;
    if ( opt.developer )
    opt.save_path = './out_dev/' ;
    end
  end

% Saving the figures
  if ~isfield( opt, 'save_path_fig' )
  opt.save_path_fig = './fig' ;
    if ( opt.developer )
    opt.save_path_fig = './fig_dev/' ;
    end
  end

% For the interpolation analysis
  if ~isfield( opt, 'star' )
  opt.star = 0 ;
  end

  if ~isfield( opt, 'planet' )
  opt.planet = 0 ;
  end

  if ~isfield( opt, 'polar' )
  opt.polar = 0 ;
  end

  if ~isfield( opt, 'super_resolution' )
  opt.super_resolution.res = 1 ;
  opt.super_resolution.interp_method = 'linear' ;
  end

  if ~isfield( opt.super_resolution, 'res' )
  opt.super_resolution.res = 1 ;
  end

  if ~isfield( opt.super_resolution, 'interp_method' )
  opt.super_resolution.interp_method = 'linear' ;
  end

  if ~isfield( opt, 'low_resolution' )
  opt.low_resolution.res = 2 ;
  opt.low_resolution.interp_method = 'linear' ;
  end

  if ~isfield( opt.low_resolution, 'res' )
  opt.low_resolution.res = 2 ;
  end

  if ~isfield( opt.low_resolution, 'interp_method' )
  opt.low_resolution.interp_method = 'linear' ;
  end

  % Number of exact simulations to derive interpolation results: for instance, 2*N+1, -N, -N+1, ..., -1, 0, 1, ..., N-1, N
  if ~isfield( opt, 'n_basis_interpolation' )
  opt.n_basis_interpolation = 5 ;
  end

  % Step in mas for a series of simulations
  if ~isfield( opt, 'step_mas' )
  opt.step_mas = 5 ;
  end
 
% For the simulation work

  if ~isfield( opt, 'input_image' )
  opt.input_image.exozodi = 1 ;
  end
  
  % Some input images
  if ~isfield( opt.input_image, 'exozodi' )
  opt.input_image.exozodi = 1 ;
  end

% Rotating Starshade
  if ~isfield( opt, 'starshade_rotation_rad' )
  opt.starshade_rotation_rad = 0 ;
  end

%%%%%%%%% Filename where to store the results %%%%%%%
  if ~isfield( opt, 'save_filename' ) || ( rplc_fl_nm )
  opt.save_filename = sprintf( 'starshade_UtotL_Nx_%i_pix_%04i_%04i_dl_%inm_dr_%3.1f_mas_psi_%3.1f_deg', ...
  opt.Nx_pupil_pix, opt.lambda_1_nm, opt.lambda_2_nm, opt.delta_lambda_nm, opt.r_source_mas, opt.psi_source_deg ) ;
    if strcmp( opt.pupil_file, '0' ) == 1, opt.save_filename = sprintf( '%s_ideal', opt.save_filename ) ; end
    if opt.diam_img_mas ~= 2000, opt.save_filename = sprintf( '%s_diam_%04i', opt.save_filename, opt.diam_img_mas ) ; end
  % Rotating Starshade
    if ( opt.starshade_rotation_rad ~= 0 )
    str_alph = strrep( sprintf( '%1.4f', opt.starshade_rotation_rad ), '.', 'p' ) ;
    opt.save_filename = sprintf( '%s_alpha_%s', opt.save_filename, str_alph ) ; 
    end
  % Different occulter file than NI2
    if ~strcmp( opt.occulter_name, 'NI2' )
    opt.save_filename = sprintf( '%s_%s', opt.save_filename, opt.occulter_name ) ;
    end
% Not used for now.
%  if opt.Nx_img ~= 400, opt.save_filename = sprintf( '%s_Nx_img_%04i', opt.save_filename, opt.Nx_img ) ; end
  end
