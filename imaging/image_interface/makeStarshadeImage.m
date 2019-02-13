function [ UTotL lambdaIn vecPetalArray pupil telescopeDiameter ] = makeStarshadeImage( opt_in )

% makeStarshadeImage
% A sample program which creates a locus of edge point from an apodization
% function and propagates the resulting field to a telescope focal plane.
% History:
% 4/25/17: first complete version, Eric Cady (JPL/Caltech)
% 4/29/17: modification of the interface, pupil handling, ..., Sergi Hildebrandt (JPL/Caltech)
% 12/12/17: further interface modifications (lambda range), UTotL computed separately, etc... (Sergi Hildebrandt)
% TBD:
% 1) Include RGB colors

% Check some options are being passed
  if ~exist( 'opt_in', 'var' )
  disp( '(makeStarshadeImage) Provide some options: opt.x=y ; makeStarshadeImage( opt ) ; returning ...' )
  return
  end

% Get default options (look inside the function for specific definitions)
opt = get_default_options( opt_in ) ;

% Developer version? Here is just this
  if ( opt.developer ), units_image ; else, units ; end

% Main parameters
Nx = opt.Nx_pupil_pix ;
dlt_lmbd = opt.delta_lambda_nm ;
r_src = opt.r_source_mas ;
psi_src = opt.psi_source_deg ;
ppl_fl = opt.pupil_file ;

% Settings for saving fields
useSave = opt.save ; % 1 save, 0 don't
savePath = opt.save_path ;
  if ~isdir( savePath ), mkdir( savePath ) ; end
% Skipping the simulation if it is saved and does not need to be re-done
if useSave == 1
  if ( ~opt.redo ) && ( exist( [ savePath '/' opt.save_filename '.mat' ] ) == 2 )
  disp( sprintf( '(makeStarshadeImage) Simulation %s exists. Skipping.', opt.save_filename ) )
  load( [ savePath '/' opt.save_filename '.mat' ] )
  return
  end
  if ( opt.redo ) && ( exist( [ savePath '/' opt.save_filename '.mat' ] ) == 2 )
  disp( sprintf( '(makeStarshadeImage) Simulation %s exists, but re-doing it.', opt.save_filename ) )
  end
end

%---------------------------
% Step 1: Load up starshade
%---------------------------
  if ( opt.developer )
  units_image
  else
  units
  end

% Load occulter definition (change this for your file structure)
occulterName = [ opt.path_occulter opt.occulter_name '.mat' ] ;
load( occulterName ); % Load up the comparison occulter
disp( sprintf( '(makeStarshadeImage) Using: %s', occulterName ) )
% Some perturbed files for the shape of the petals do not come with Z
% NI2
NI2_prt = 0 ;
  if numel( strfind( opt.occulter_name, 'NI2' ) ) && ( numel( opt.occulter_name ) > 3 )
  NI2_prt = 1 ;
  end
  if ( NI2_prt )
  Z = 3.724225668350351e+07 ; % m (425-552 nm)
  telescopeDiameter = 2.4 ;
  end
% NW2
NW2_prt = 0 ;
  if numel( strfind( opt.occulter_name, 'NW2' ) ) && ( numel( opt.occulter_name ) > 3 )
  NW2_prt = 1 ;
  end
  if ( NW2_prt )
  Z = 1.197666616918624e+08 ; % m (1000 nm)
  telescopeDiameter = 4 ;
  end

% Adapting the distance of the starshade depending on the band
Z = set_starshade_distance( Z, opt ) ;

% Range of wavelengths to be produced
lambdaIn = nm * ( opt.lambda_1_nm + dlt_lmbd * ( 0 : 1 : ( opt.lambda_2_nm - opt.lambda_1_nm ) / dlt_lmbd ) ) ;
  if lambdaIn( end ) ~= ( nm * opt.lambda_2_nm ) 
  lambdaIn( end + 1 ) = lambdaIn( end ) + dlt_lmbd * nm ;
  end
n_lmbd = numel( lambdaIn ) ;
disp( sprintf( 'Considering %i wavelengths', numel( lambdaIn ) ) ) 
% Build/load
if strcmp( ppl_fl, '0' )
    secondarySize = opt.secondary_size ; % If you're fine with a generic circular pupil with a secondary (but no struts), set this to the fraction of the radius the secondard covers (e.g. 0.2)
    if ( opt.developer )
    pupil = makePupil_image( Nx, Nx, 1, secondarySize, 0, 0 ) ;
    else
    pupil = makePupil(Nx, Nx, 1, secondarySize, 0, 0);
    end
disp( sprintf( 'Considering a circular pupil with a secondary blocking a fraction of the radius=%f', secondarySize ) ) ;
else
    if strfind( ppl_fl, '.fits' )
    pupil = fitsread( ppl_fl );
    disp( [ 'File for the pupil read: ' ppl_fl ] ) 
    else
    load( ppl_fl ) ;
    % Test from NG
    pupil = pupil_mask ;
    disp( [ 'File for the pupil read: ' ppl_fl ] )
    end
    n_ppl = sqrt( numel( pupil ) ) ;
    % Check
    disp( sprintf( 'Reference pupil size is %ix%i pixels.', n_ppl, n_ppl ) ) 
    % Reducing the size of the pupil grid (fast<1s)
      if Nx ~= n_ppl
      % finding where the pupil starts
      q_ppl_1 = min( find( pupil ~= 0 ) ) ;
      q_ppl_2 = max( find( pupil ~= 0 ) ) ;
      % Column where the pupil data starts and ends
      clmn_1 = floor( q_ppl_1 / n_ppl ) ;
      clmn_2 = floor( q_ppl_2 / n_ppl ) ;
      % Actual resize (method bicubic, some low pass filter is applied first with 11 pix box over the 2048, fine)
      %% Asssuming symmetric pupil when cutting it out
      pupil = resizem( pupil( clmn_1 : clmn_2, clmn_1 : clmn_2 ), Nx / ( clmn_2 - clmn_1 + 1 ), 'bilinear' ) ; 
      % Sharping the edges
      pupil( find( pupil < .7 ) ) = 0 ;
      pupil( find( pupil ~= 0 ) ) = 1 ;
      disp( sprintf( 'New pupil size is reduced to %ix%i pixels.', Nx, Nx ) )
      end
% Exception for checks
%pupil = fitsread( 'pupil_D1Kpix_64.fits' );
end
% Lateral offsets in meters
deltaX = 0;
deltaY = 0;

% -------------------------------
% Step 2: Build edge
% -------------------------------

% Build edge; this function is overkill for an unaberrated edge but will do
% the job
tic
  % Case without perturbations
  vecPetalArray = NaN ;
  % The number of petals is set by default to 24, but it can be changed
  if ~( NI2_prt ) && ~( NW2_prt )
    if ( opt.developer )
    vecPetalArray = createErrorProfileV1_image(r, Profile, occulterDiameter, petalLength, opt.n_ptl, {});
    else
    vecPetalArray = createErrorProfileV1(r, Profile, occulterDiameter, petalLength, opt.n_ptl, {});
    end
  t = toc ;
    if ( opt.developer )
    createErr_lbl = 'createErrorProfileV1_image' ;
    else
    createErr_lbl = 'createErrorProfileV1' ;
    end
  disp( sprintf( '%s took %3.2f seconds', createErr_lbl, t ) )

  tic
  tmpxVals = [];
  tmpyVals = [];
  tmpzVals = [];%**
  for j = 1 : opt.n_ptl
      tmpxVals = [tmpxVals vecPetalArray{j}{1}(1, :)]; %#ok<AGROW>
      tmpyVals = [tmpyVals vecPetalArray{j}{1}(2, :)]; %#ok<AGROW>
      tmpzVals = [tmpzVals vecPetalArray{j}{1}(3, :)]; %#ok<AGROW> %**
  end
  xVals = [tmpxVals tmpxVals(1)];
  yVals = [tmpyVals tmpyVals(1)];
  zVals = [tmpzVals tmpzVals(1)];
  else
  % New xVals, yVals
  load( [ opt.path_occulter opt.occulter_name ] )
    if ~exist( 'zVals', 'var' )
    zVals = 0 * xVals ;
    end
  disp( sprintf( 'Loaded new xVals and yVals from %s', opt.occulter_name ) )
  % Sometimes the x/y/zVals are given as #x1 arrays instead of 1x#
    if size( xVals, 1 ) ~= 1
    xVals = xVals' ;
    end
    if size( yVals, 1 ) ~= 1
    yVals = yVals' ;
    end
    if size( zVals, 1 ) ~= 1
    zVals = zVals' ;
    end
  end

  % Closing the polygon, if it's open
  if (xVals(1) ~= xVals(end)) || (yVals(1) ~= yVals(end))
  xVals = [xVals(:); xVals(1)];
  yVals = [yVals(:); yVals(1)];
  end


% Simple rotation in the XY plane
  if ( opt.starshade_rotation_rad )
  disp( sprintf( '*** Rotating the Starshade an angle %d deg', opt.starshade_rotation_rad * 180 / pi ) ) ;
  alph_rt = opt.starshade_rotation_rad ;
  xVals_nw = xVals * cos( alph_rt ) + yVals * sin( alph_rt ) ;
  yVals_nw = -xVals * sin( alph_rt ) + yVals * cos( alph_rt ) ;
  xVals = xVals_nw ;
  yVals = yVals_nw ;
  clear xVals_nw yVals_nw
  end

% xVals, yVals, zVals give the 3D edge locus
% Under most circumstances zVals will be all zeros
t = toc ;
disp( sprintf( 'xyzVals took %3.2f seconds', t ) )

%--------------------------
% Step 3: Compute field at 
%--------------------------

% Propagate to telescope aperture plane with line integral
tic
  if ( opt.developer )
  UTotL = bdwf_image(xVals, yVals, zVals, Z, lambdaIn, telescopeDiameter/Nx, Nx, r_src * mas, psi_src * degree, deltaX, deltaY);
  else
  UTotL = bdwf(xVals, yVals, zVals, Z, lambdaIn, telescopeDiameter/Nx, Nx, r_src * mas, psi_src * degree, deltaX, deltaY);
  end
t = toc;
  if ( opt.developer )
  bdwf_lbl = 'bdwf_image' ;
  else
  bdwf_lbl = 'bdwf' ;
  end

disp( sprintf('%s took %3.2f seconds, %3.2f per wavelength bin', bdwf_lbl, t, t / n_lmbd ) )

  if useSave == 1
  % These fields do not belong to this computation and may vary later on
    if isfield( opt, 'Nx_image_pix' ), opt = rmfield( opt, 'Nx_image_pix' ) ; end
    if isfield( opt, 'Nx_img' ), opt = rmfield( opt, 'Nx_img' ) ; end
    if isfield( opt, 'diam_img_mas' ), opt = rmfield( opt, 'diam_img_mas' ) ; end
  % Avoiding issues with the opt
  opt_make = opt ;
  save([savePath '/' opt.save_filename '.mat' ], 'lambdaIn', 'opt_make', 'pupil', 'telescopeDiameter', 'UTotL' )
  end

dbstop if error
%make_an_error

function Z_new = set_starshade_distance( Z, opt )
% Default do nothing
Z_new = Z ;
Z_0 = Z ;
% WFIRST
  if strfind( opt.occulter_name, 'NI2' )
  % Control there's a change if NI2 is used
  Z_new = 0 ;
  % There are three bands
  % For the default band: 425-552 nm, Z is 3.724225668350351e+07 m, IWA=72 mas
    if ( opt.lambda_1_nm >= 425 ) && ( opt.lambda_2_nm <= 552 ) && ( opt.lambda_1_nm <= opt.lambda_2_nm )
    Z_new = Z ;
    end
  % Longer wavelength, closer distance Z=2.612037672633376e+07 m, greater IWA, IWA=102.6571 mas
    if ( opt.lambda_1_nm >= 606 ) && ( opt.lambda_2_nm <= 787 ) && ( opt.lambda_1_nm <= opt.lambda_2_nm )
    Z_new = Z * ( 425 + 552 ) / ( 606 + 787 ) ; 
    end
  % Z= 2.119142969119565e+07, IWA=126.5343 mas
    if ( opt.lambda_1_nm >= 747 ) && ( opt.lambda_2_nm <= 970 ) && ( opt.lambda_1_nm <= opt.lambda_2_nm )
    Z_new = Z * ( 425 + 552 ) / ( 747 + 970 ) ;
    end

  % Final message if successful
    if Z_new
      if Z_new ~= Z_0
      disp( sprintf( 'New distance between the telescope and the Starshade changed from %f km to %f km. New IWA=%3.2f', Z / 1000, Z_new / 1000, 72 * Z / Z_new ) ) ;
      else
      disp( sprintf( 'Distance between the telescope and the Starshade is %f km', Z / 1000 ) ) ;
      end
    end
  
  dbstop if error
  % Checking there was some change
    if ~Z_new 
    disp( sprintf( 'Walength is not within one of the three expected bands: 425-552 nm, 606-787 nm, 747-970 nm. The minimum wavelength set is %04f nm, and the maximum wavelength is %04f nm', opt.lambda_1_nm, opt.lambda_2_nm ) )
    make_an_error
    end
  else
  disp( sprintf( 'Distance between the telescope and the Starshade is %f km', Z / 1000 ) ) ;
  end


