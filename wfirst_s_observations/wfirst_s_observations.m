% Function to plot the range of possible observations with WFIRST+Starshade = WFIRST-S
function opt = wfirst_s_observations( opt )

% History:
% May 2016: First version by Sergi Hildebrandt
% June 2016: Addition of ExoCat, Kevin Tang
% June 2016: Addition of background disambiguation, Sergi Hildebrandt & Kevin Tang


% Examples:
% PS: notice that the variable 'opt' is a structure with information about the mission and the targets that can be combined to produce other plots/summaries.
% 1) Basic example (default options)
% clear ; opt.target.label = 'number'; opt.plot_pm=1; wfirst_s_observations( opt ) ;
% 2) Plotting the results in Ecliptic coordinates:
% clear ; opt.target.label = 'number'; opt.plot_pm=1; opt.coord='ecl'  ; opt = wfirst_s_observations( opt ) ;
% 3) Saving the image and shifting the mission day from 0 to 134:
% clear ; opt.target.label = 'number'; opt.plot_pm=1; opt.save_image= 1 ; opt.day_of_mission = 134 ; opt.coord = 'radec' ; wfirst_s_observations( opt ) ;
% Slightly less coverage: some objects become un-observed
% 4) clear ; opt.target.label = 'number'; opt.plot_pm=1; opt.coord='ecl'; opt.beta_1 = 60 ; opt.beta_2 = 75 ;opt=wfirst_s_observations(opt) ;
% Almost complete coverage: same observational time per year
% 5) clear ; opt.target.label = 'number'; opt.plot_pm=1; opt.coord='ecl'; opt.beta_1 = 0.1 ; opt.beta_2 = 89.9 ;opt=wfirst_s_observations(opt) ;
% Very narrow coverage
% 6) clear ; opt.target.label = 'number'; opt.plot_pm=1; opt.coord='ecl'; opt.beta_1 = 65 ; opt.beta_2 = 70 ;opt=wfirst_s_observations(opt) ;
% Only plotting parallax motion (no proper motion)
% clear ; opt.target.label = 'number'; opt.plot_pm=1; opt.coord='ecl'; opt.pax_only = 1 ;opt=wfirst_s_observations(opt) ;
% 7) Get a figure of merit based on background disambiguation across all ExoCat Catalog
% clear ;  opt.target.label = 'number'; opt.plot_pm=1; opt.scan_exocat = 1 ; wfirst_s_observations( opt ) ;

% The reference coordinates are J2000. For instance, a final position is the result of its value in J2000 plus the proper motion that would tell us where
% to find its new position in J2000 coordinates, and with parallax added, where it is in that time of the mission again in terms of J2000 coordinates. 
% Obviously it will be in a different J2000 position that the official star's position in J2000 since the time of observation is not 2000.
% Motion of the starshade increasing Ecliptic longitude, like Earth.
% Q1: The orbit of WFIRST-S is assumed to be circular (e=0.016, OK) but the effect of other planets on the precise location of the Earth from 2000, assumes
% constant orbital angular velocity (might not be accurate for some parallaxes?)
% Optional improvement: use the coordinate transformations tools from:
% https://webhome.weizmann.ac.il/home/eofek/matlab/fun_content.html
% Nota Bene: Matlab cannot read binary Excel sheets unless in Windows OS. So, this is ExoCat.xlsb, then saved as xlsx.

% Quick summary:
% Circular orbit at L2, locked with Earth. No particular orbit around L2.
% It assumes WFIRST is orbiting in the direction of *increasing* Ecliptic longitude in time
% It assumes two limiting solar avoidance angles:
% opt.beta_1=54 ; degrees (default value 54 degrees)
% opt.beta_2=83 ; degrees (default value 83 degrees)
% It assumes coordinates of objects are given in J2000. Later projected into any date supllied by the user, e.g., 2028
% Quick help
% The initial position of WFIRST-S in the orbit can be set with the angle alpha_0 measured in ecliptical cooordinates:
% opt.alpha_0 = 100.41670324 ; (default value 100.41670324 degrees) The Heliocentric Ecliptic coordinates for the Earth for 2028 1st January at 0h0m0s is 100 deg 25' 07.9'' (http://cosinekitty.com/solar_system.html). NB: Notice that this code assumes constant orbital speed (with Earth's period computed below) and does not consider leap-years and regular years.
% The starting time of the observation (for coordinate corrections)
% opt.first_light_date = 2028 ; (default value year 2028.000. PS: it can be decimal, to include any date, e.g,  2028.123)
% Days into the mission (day to be added to the first light date)
% opt.day_of_mission = 0 ; days
% opt.end_light_date = 2030 ; (default 2030.000. PS: it can be decimal)
% Steps for the grid of coordinates
% opt.grid = 0.5 ; degrees (default, 0.5 degrees)
% Coordinate system to display the results
% opt.coord = 'radec' ; (default value, radec. Other cases: 'gal', 'ecl')
% In case, the observation window is to be computed for a particular object
% opt.inv = 0 ; (default, 0. If 1, it will do it, that is invert the pointing solution)
% Object to be considered (list of objects):
% opt.target = { 'tau_ceti' } ; (default, 'tau_ceti'. A string, see the sub-function)
% The overall score in the figure of merit for the target (list of objects)
% opt.target.fom = [ 1 ] % default 1
% Alternatively, the *Equatorial* coordinates of the object can be given (list of objects):
% opt.target.radec = { [ a, b ] } ; degrees (coordinate system set by opt.coord) 
% Proper motion (list of objects)
% opt.target.pm = { [ ma, mb ] } ; mas/yr (in RA*cos(DEC), DEC). Notice the convention in SIMBAD to have projected the RA motion.
% Parallax (list of objects)
% opt.target.pax = [ p ] ; mas
% The days that WFIRST-S can observe a target given the first light date and the day_of the mission.
% opt.target.day_obs_string = '' ;
% The total number of days per orbital period (year) that a target can be observed.
% opt.target.n_days_obs = 0 ;
% The number of times that the object is observable by WFIRST-S: none (0), 1 (single crossing with WFIRST-S field of view), 2 (two periods: annular circle)
% opt.target.n_window = 0 ;
% opt.target.days_list = { [ NaN, NaN, NaN, NaN ] } ;
% Pixel size camera (mas)
% opt.camera_pixel = 60*1000/1000 ; Assuming a 1 arcmin in 1000 pixels
% Optional lapse of time for plot with multi-observations
% opt.days_interval = 50 ; days
% Plot options
% opt.target.label = '' ; '': none, 'name', 'number'
% opt.pax_only=0/1 ; plot only parallax motion without adding proper motion (No/Yes=0/1)
% Plot quiver vectors
% opt.plot_pm = 0 ; 0/1: No/Yes
% Factor to enlarge proper motion vectors
% opt.pm_fact = 1e5 ;
% opt.save_image = 0 ; 0/1 (No/Yes)
% opt.path_image = 'img' ; path to store the image (if it does not exist, it will be created)

% Default options
opt = get_opt( opt ) ;
if ( 1 )
% % % % % % % % % % % % % % % % % % % % % % % % %
% Examples of some possible calls and figures:  %
% % % % % % % % % % % % % % % % % % % % % % % % % 

% 1.-Creating the grid of observed coordinates
[ LON_GRID, LAT_GRID ECL_LAT_GRID ] = observed_grid( opt ) ;

% 2.- Observation map for a given date
clf
plot_generic( LON_GRID, LAT_GRID, ECL_LAT_GRID, opt ) ;
  if ( opt.save_image )
    if ~isdir( opt.path_image ), system( [ 'mkdir -p ' opt.path_image ] ) ; end
  saveas( gca, [ opt.path_image 'wfirst-s_' opt.coord '_single_day' ], 'png' ) ;
  end

% 3.- Getting the positions at a given date for an object
% For this example, create a list of objects but it should normally be a user-supplied list of object(s)
opt = list_of_targets( opt ) ;
  for i_trgt = 1 : numel( opt.target.radec )
  [ opt.target.radec_obs{ i_trgt }( 1 ), opt.target.radec_obs{ i_trgt }( 2 ) ] = get_radec_pos_target( opt, i_trgt ) ;
  end
  
% 4.- Getting the days a target will be observable
i_trgt;
  for i_trgt = 1 : numel( opt.target.radec )
  [ Q_DFF opt.target.day_obs_string{ i_trgt } opt.target.n_days_obs( i_trgt ) opt.target.n_window( i_trgt ) opt.target.days_list{ i_trgt } ] = ...
  get_day_of_obs_string( opt.target.radec_obs{ i_trgt }( 1 ), opt.target.radec_obs{ i_trgt }( 2 ), opt ) ;
  end

% 5.- Plot of targets' positions in a period of time
% 5.1.- Individual plots
  for i_trgt = 1 : numel( opt.target.radec )
  plot_single_target_in_time( opt, i_trgt ) ;
    if ~isdir( opt.path_image ), system( [ 'mkdir -p ' opt.path_image ] ) ; end
  pth_nm_img = [ opt.path_image strrep( opt.target.name{ i_trgt }, ' ', '_' ) '_observability_by_wfirst_s' ] ;
  saveas( gca, pth_nm_img, 'png' ) ;
  disp( [ '(WFIRST_S_OBSERVATIONS) Stored image: ' pth_nm_img '.png' ] ) ;
  end
% 5.2.- Altogether
plot_targets_in_time( opt ) ;

end

% 6.- Overplotting the expected trajectory on the sky over a catalog image
% One example for now (needs Aladin Lite communication TBD)
i_tau_ceti = 1 ;
plot_target_on_catalog( opt, i_tau_ceti ); 

% 7.- Run through all targets in ExoCat and collect figure of merit
  if (opt.scan_exocat)
  get_back_disambiguation_exocat( opt ) ;
  % Calculating figure of merit
  histogram_back_disambiguation( opt ) ;
  end

% End of the main program

% Sub-functions
function opt_out = get_opt( opt_in )
opt_out = opt_in ;
% First limit angle on the orbit plane
  if ~isfield( opt_in, 'beta_1' )
  opt_out.beta_1 = 54 ; % degrees
  end
  if ~isfield( opt_in, 'beta_2' )
  opt_out.beta_2 = 83 ; % degrees
  end
  if ~isfield( opt_in, 'alpha_0' )
  opt_out.alpha_0 = 100.41670324 ; % degrees (see comments at the beginning)
  end
  if ~isfield( opt_in, 'first_light_date' )
  opt_out.first_light_date = 2028 ; % year date, can be fractional
  end
  if ~isfield( opt_in, 'end_light_date' )
  opt_out.end_light_date = 2030 ; % year date, can be fractional
  end
  if ~isfield( opt_in, 'day_of_mission' )
  opt_out.day_of_mission = 0 ; % year date, can be fractional
  end
  if ~isfield( opt_in, 'grid' )
  opt_out.grid = 0.5 ; % degrees
  end
  if ~isfield( opt_in, 'coord' )
  opt_out.coord = 'radec' ; %
  end
  if ~isfield( opt_in, 'inv' )
  opt_out.inv = 0 ; %
  end
  if ~isfield( opt_in.target, 'name' )
  opt_out.target.name = { 'tau_ceti' } ; % string with the name of the star or system
  end
  if ~isfield( opt_in.target, 'fom' )
  opt_out.target.fom = [ 1 ] ; % It can be a list of objects
  end
  if ~isfield( opt_in.target, 'radec' )
  opt_out.target.radec = { [ 0, 0 ] } ; % degrees. It can be a list of objects
  end
  if ~isfield( opt_in.target, 'pm' )
  opt_out.target.pm = { [ 0, 0 ] } ; % degrees. It can be a list of objects
  end
  if ~isfield( opt_in.target, 'pax' )
  opt_out.target.pax = [ 0 ] ; % degrees. It can be a list of objects
  end
  if ~isfield( opt_in.target, 'day_obs_string' )
  opt_out.target.day_obs_string = '' ; 
  end
  if ~isfield( opt_in.target, 'n_days_obs' )
  opt_out.target.n_days_obs = 0 ; 
  end
  if ~isfield( opt_in.target, 'n_window' )
  opt_out.target.n_window = 0 ; 
  end
  if ~isfield( opt_in.target, 'days_list' )
  opt_out.target.days_list = { [ NaN, NaN, NaN, NaN ] } ;
  end
  if ~isfield( opt_in, 'camera_pixel' )
  opt_out.camera_pixel = 60 * 1000 / 1000 ; % mas. Default 1 arcmin in 1000 pixels
  end
  if ~isfield( opt_in, 'days_interval' )
  opt_out.days_interval = 90 ; % days
  end
  if ~isfield( opt_in.target, 'label' )
  opt_out.target.label = '' ; % '': none, 'name', 'number'
  end
  if ~isfield( opt_in, 'pax_only' )
  opt_out.pax_only = 0 ; % 0/1
  end
  if ~isfield( opt_in, 'plot_pm' )
  opt_out.plot_pm = 0 ; %
  end
  if ~isfield( opt_in, 'pm_fact' )
  opt_out.pm_fact = 1e5 ; %
  end
  if ~isfield( opt_in, 'save_image' )
  opt_out.save_image = 0 ; %
  end
  if ~isfield( opt_in, 'path_image' )
  opt_out.path_image = 'img/' ; % 
  end
  if ~isfield (opt_in, 'scan_exocat')
  opt_out.scan_exocat = 0; %
  end 
  if ~isfield( opt_in, 'lambda' )
  opt_out.lambda = 750 ; %nanometers
  end
  if ~isfield( opt_in, 'diam' )
  opt_out.diam = 2.4 ; %meters
  end
  if ~isfield( opt_in, 'snr' )
  opt_out.snr = 5 ; %
  end
  if ~isfield(opt_in, 'snr_time')
  opt_out.snr_time = 1; %days to get to snr
  end
  if ~isfield(opt_in, 't_obs')
  opt_out.t_obs = 5; %min days needed to observe target
  end
% Some check
  if ~strcmp( opt_out.coord, 'gal' ) && ~strcmp( opt_out.coord, 'ecl' ) && ~strcmp( opt_out.coord, 'radec' )
  disp( sprintf( 'WARNING: please, choose opt.coord=%s, %s or %s, or do not set it (default are Ecliptic). You have set opt.coord=%s. Entering debug mode with keyboard ...', 'radec', 'ecl', 'gal', opt_out.coord ) ) ; 
  keyboard
  end

function [ LON_GRID, LAT_GRID, ECL_LAT_GRID ] = observed_grid( opt )

% Creating the grid of coordinate points. The region that can be observed with Starshade and WFIRST is a circular sector on the sphere. 
% 1.-Considering the axis Sun-WFIRST first:
PHI =  -0.1 : opt.grid : 360.1 ;
THETA = (90 - opt.beta_2 ) : opt.grid : ( 90 - opt.beta_1 ) ;
[ PHI_GRID THETA_GRID ] = meshgrid( PHI, THETA ) ;
PHI_GRID = PHI_GRID * pi / 180 ;
THETA_GRID = THETA_GRID * pi / 180 ;

% 2.-Rotating the coordinate system by +90 degrees to have it in Ecliptic coordinates 
COLAT_GRID = acos( -cos( THETA_GRID ) .* cos( PHI_GRID ) ) ;
LON_GRID = atan2( cos( THETA_GRID ) .* sin( PHI_GRID ), sin( THETA_GRID ) ) ;
LAT_GRID = 90 - COLAT_GRID * ( 180 / pi ) ;
LON_GRID = LON_GRID * ( 180 / pi ) ;

% If a day of the mission is provided, then it means that a new starting ecliptic longitude has been chosen
% Assuming Starshade is orbiting in the direction of increasing ecliptic longitude
% The addition of 180 degrees reflects the fact that the region where the starshade can be located is opposite to the vector Sun-Telescope (WFIRST), which coincides with Earth's ecliptic longitude, in orer to have an L2 orbit.
LON_GRID = LON_GRID + opt.alpha_0 + 180 ;
  if ( opt.day_of_mission )
  disp( sprintf( 'One day of the mission has been chosen: %d. Changing the initial ecliptic longitude.', opt.day_of_mission ) )
  LON_GRID = LON_GRID + 360 / l2_period() * opt.day_of_mission ;
  end

% For plotting purposes only:
ECL_LAT_GRID = LAT_GRID ;
% Matlab is column-major order
ECL_LAT_GRID = ECL_LAT_GRID' ;

% Resizing the mesh for the change of coordinates
LON_GRID = reshape( LON_GRID, 1, numel( LON_GRID ) ) ;
LAT_GRID = reshape( LAT_GRID, 1, numel( LAT_GRID ) ) ;

% Change to RADEC if requested
  if strcmp( opt.coord, 'radec' ) || strcmp( opt.coord, 'gal' )
  [ LON_GRID, LAT_GRID ] = eclip2equat( LON_GRID * pi / 180, LAT_GRID * pi / 180 ) ;
  LON_GRID = LON_GRID * ( 180 / pi ) ;
  LAT_GRID = LAT_GRID * ( 180 / pi ) ;
  end

% Change to Galactic coordinates if requested
  if strcmp( opt.coord, 'gal' )
  [ LON_GRID, LAT_GRID ] = RA_Dec_to_Gal( LON_GRID, LAT_GRID ) ;
  end

% Again as a mesh
LON_GRID = reshape( LON_GRID, numel( THETA ), numel( PHI ) ) ;
LAT_GRID = reshape( LAT_GRID, numel( THETA ), numel( PHI ) ) ;

% Transpose to have 'bins in longitude'x'bins in latitude'
LON_GRID = LON_GRID' ;
LAT_GRID = LAT_GRID' ;

% Coordinate changes 
% Ecliptic to Equatorial
% Dmitry Savransky
% From https://www.mathworks.com/matlabcentral/fileexchange/23285-conversion-between-equatorial-and-ecliptic-coordinates
function [alpha,delta] = eclip2equat(lambda,beta)
% ECLIP2EQUAT - Convert from ecliptic to equatorial celestial coordinates.
%   [alpha,delta] = eclip2equat(lambda,beta) returns the right ascension 
%   and declination in the equatorial coordinate system corresponding
%   to the longitudinal and latitudinal coordinates in the ecliptic
%   coordinate system (lambda and beta)
% 
%   All inputs and outputs are in units of radians.  alpha will be in the
%   range of [-pi,pi] while delta will be in the range [-pi/2, pi/2]. 
%   Inputs must be vectors of equal size;
% 
%   Example:
%       [alpha, delta] = eclip2equat([ -0.2295;-0.5555],[0.1517;0.2792]);
 
%   Written by Dmitry Savransky, 3 March 2009

%error checking
if nargin ~= 2
    error('You must input both Declination and Right Ascension.');
end
lambda = lambda(:);
beta = beta(:);
if length(lambda) ~= length(beta)
    error('Inputs must be the same length.');
end

%define Earth's axial tilt
epsilon = (23 + 26/60 + 21.448/3600)*pi/180; %radians

%calculate trigonometric combinations of coordinates
sd = sin(epsilon)*sin(lambda).*cos(beta) + cos(epsilon)*sin(beta);
cacd = cos(lambda).*cos(beta);
sacd = cos(epsilon)*sin(lambda).*cos(beta) - sin(epsilon).*sin(beta);

%calculate coordinates
alpha = atan2(sacd,cacd);
r = sqrt(cacd.^2 + sacd.^2);
delta = atan2(sd,r);
r2 = sqrt(sd.^2+r.^2);
% SRH: making the angle belong to [ 0, 2pi )!
alpha = mod( alpha + 4*pi, 2*pi ) ;
%sanity check: r2 should be 1
if sum(abs(r2 - 1)) > 1e-11 ; % 1e-12 originally, a bit exaggerated
    warning('equat2eclip:radiusError',...
        'Latitude conversion radius is not uniformly 1. ');
end

% Equatorial to Ecliptic
%  Dmitry Savransky
% From https://www.mathworks.com/matlabcentral/fileexchange/23285-conversion-between-equatorial-and-ecliptic-coordinate
function [lambda,beta] = equat2eclip(alpha,delta)
% EQUAT2ECLIP - Convert from equatorial to ecliptic celestial coordinates.
%   [lambda,beta] = equat2eclip(alpha,delta) returns the longitudinal and
%   latitudinal coordinates in the ecliptic coordinate system corresponding
%   to coordinates of right ascension and declination (alpha and delta) in
%   the equatorial coordinate system.
% 
%   All inputs and outputs are in units of raidans.  lambda will be in the
%   range of [-pi,pi] while beta will be in the range [-pi/2, pi/2]. Inputs
%   must be vectors of equal size;
% 
%   Example:
%       [lambda,beta] = equat2eclip([-0.2700;-0.6132],[0.0492;0.0512]);
 
%   Written by Dmitry Savransky, 3 March 2009

%error checking
if nargin ~= 2
    error('You must input both Declination and Right Ascension.');
end
delta = delta(:);
alpha = alpha(:);
if length(delta) ~= length(alpha)
    error('Inputs must be the same length.');
end

%define Earth's axial tilt
epsilon = (23 + 26/60 + 21.448/3600)*pi/180; %radians

%calculate trigonometric combinations of coordinates
sb = sin(delta)*cos(epsilon) - cos(delta)*sin(epsilon).*sin(alpha);
cbcl = cos(delta).*cos(alpha);
cbsl = sin(delta)*sin(epsilon) + cos(delta)*cos(epsilon).*sin(alpha);

%calculate coordinates
lambda = atan2(cbsl,cbcl);
r = sqrt(cbsl.^2 + cbcl.^2);
beta = atan2(sb,r);
r2 = sqrt(sb.^2+r.^2);
% SRH: making the angle be in [ 0, 360 )!
lambda = mod( lambda + 4*pi, 2*pi ) ;

%sanity check: r2 should be 1
if sum(abs(r2 - 1)) > 1e-12
    warning('equat2eclip:radiusError',...
        'Latitude conversion radius is not uniformly 1. ');
end


% Equatorial to Galactic
function [l,b]=RA_Dec_to_Gal(RA,Dec)
% RA_Dec_to_Gal: Convert RA and Dec to galactic coordinates
%
% [l,b]=RA_Dec_to_Gal(RA,Dec)
%
% ARGUMENTS
%  RA    The Right Ascension [degrees]
%  Dec    The Declination [degrees]
%
% RETURNS
%  l    Galactic longitude [degrees]
%  b    Galactic latitude [degrees]
%
% NOTES
%  Presuming J2000.
%
%  Uses the approximate formulae from "Allen's Astrophysical
%  Quantities", Aurthur N. Cox, Ed. 2000.
%  C.A. Murray, 1988, A&Ap 218, 325 supposedly has a more precise formula.
%
% See also: PRECESS

% AUTHOR: Eric Tittley
%
% HISTORY:
%  00 11 17 First version
%  07 11 05 Modified Comments
%
% COMPATIBILITY: Matlab, Octave
%
% LICENSE
%  Copyright Eric Tittley 2000
%  See the Gnu GPL license.
%  Essentially, free to use.  Free to modify.  But cannot be re-created in
%  whole or in part in anything sold or traded for something of worth.

Dec=Dec/180*pi;
RA=RA/180*pi;

Coef1=62.87/180*pi;
Coef2=282.86/180*pi;
Coef3=32.9319186/180*pi;
Coef4=2*pi+Coef3;

b=asin(sin(Dec)*cos(Coef1)-cos(Dec).*sin(RA-Coef2)*sin(Coef1));
l=0*b;

% Catch where cos(b)==0 and set the values by hand.  These correspond to the
% Galactic poles.
b_eq_2pi = find(cos(b)==0);
b(b_eq_2pi) = b(b_eq_2pi)*0+pi/2;
l(b_eq_2pi) = l(b_eq_2pi)*0;

% This is the break that separates into hemispheres the domains for which we
% have to get the signs right.  We also catch cos(b)==0, which would give a
% divide by zero error.
Dec_le_break = find( (Dec <= (Coef1-pi/2)*sin(RA-Coef2)) & cos(b)~=0 );
Dec_gt_break = find( (Dec  > (Coef1-pi/2)*sin(RA-Coef2)) & cos(b)~=0 );

l(Dec_gt_break)=Coef3 + ...
 acos( cos(Dec(Dec_gt_break)).*cos(RA(Dec_gt_break)-Coef2)./cos(b(Dec_gt_break)));

l(Dec_le_break)=Coef4 - ...
 acos( cos(Dec(Dec_le_break)).*cos(RA(Dec_le_break)-Coef2)./cos(b(Dec_le_break)));

l_ge_2pi = find(l>=2*pi);
l(l_ge_2pi)=l(l_ge_2pi)-2*pi;

b=b*180/pi;
l=l*180/pi;

% Some generic plot
function h_plt = plot_generic( LON_GRID, LAT_GRID, ECL_LAT_GRID, opt ) 
figure( 1 ) ; clf
setwinsize(gcf, 980,340);
hold all
% Dealing with the circular symmetry. 
LON_GRID = LON_GRID  + 360 - floor( max( LON_GRID( : ) / 360 ) ) * 360 ;
q = find( LON_GRID > 360 ) ;
  if numel( q ) ~= 0
  LON_GRID_1 = LON_GRID * nan ;
  LAT_GRID_1 = LAT_GRID * nan ;
  ECL_LAT_GRID_1 = ECL_LAT_GRID * nan ; 
  LON_GRID_1( q ) = LON_GRID( q ) - 360 ;
  LAT_GRID_1( q ) = LAT_GRID( q ) ;
  ECL_LAT_GRID_1( q ) = ECL_LAT_GRID( q ) ;
  % Setting off warning message about nan
  warning( 'off', 'MATLAB:contour:NonFiniteData' )
  % Spacing between angles is 5 degrees
  contour( LON_GRID_1, LAT_GRID_1, ECL_LAT_GRID_1, 10 ) ;
  xlim( [ 0, 360 ] ) ;
  ylim( [ -90, 90 ] ) ;
  % Borders (the equal longitude limits are unnecessary. Commented out.
  %contour( LON_GRID_1(1:2, : ), LAT_GRID_1(1:2, : ), ECL_LAT_GRID_1( 1:2, : ), 600, 'k' ) ;
  %contour( LON_GRID_1(end-1:end, : ), LAT_GRID_1(end-1:end, : ), ECL_LAT_GRID_1( end-1:end, : ), 600, 'k' ) ;
  contour( LON_GRID_1(:, 1:2 ), LAT_GRID_1( :, 1:2 ), ECL_LAT_GRID_1( :, 1:2 ), 600, 'k' ) ;
  contour( LON_GRID_1(:, end-1:end ), LAT_GRID_1( :, end-1:end ), ECL_LAT_GRID_1( :, end-1:end ), 600, 'k' ) ;
  clear q LON_GRID_1 LAT_GRID_1 ECL_LAT_GRID_1
  q = find( LON_GRID < 360 ) ;
    if numel( q ) ~= 0
    LON_GRID_1 = LON_GRID * nan ;
    LAT_GRID_1 = LAT_GRID * nan ;
    ECL_LAT_GRID_1 = ECL_LAT_GRID * nan ;
    LON_GRID_1( q ) = LON_GRID( q ) ;
    LAT_GRID_1( q ) = LAT_GRID( q ) ;
    ECL_LAT_GRID_1( q ) = ECL_LAT_GRID( q ) ;
    % Spacing between angles is 5 degrees
    contour( LON_GRID_1, LAT_GRID_1, ECL_LAT_GRID_1, 10 ) ;
    xlim( [ 0, 360 ] ) ;
    ylim( [ -90, 90 ] ) ;
    % Borders (the equal longitude limits are unnecessary. Commented out.
    %contour( LON_GRID_1(1:2, : ), LAT_GRID_1(1:2, : ), ECL_LAT_GRID_1( 1:2, : ), 600, 'k' ) ;
    %contour( LON_GRID_1(end-1:end, : ), LAT_GRID_1(end-1:end, : ), ECL_LAT_GRID_1( end-1:end, : ), 600, 'k' ) ;
    contour( LON_GRID_1(:, 1:2 ), LAT_GRID_1( :, 1:2 ), ECL_LAT_GRID_1( :, 1:2 ), 600, 'k' ) ;
    contour( LON_GRID_1(:, end-1:end ), LAT_GRID_1( :, end-1:end ), ECL_LAT_GRID_1( :, end-1:end ), 600, 'k' ) ;
    clear q LON_GRID_1 LAT_GRID_1 ECL_LAT_GRID_1
    end
  end

xlim( [ 0, 360 ] ) ;
ylim( [ -90, 90 ] ) ;
LBL_COORD = { 'ECLIPTIC LON', 'ECLIPTIC LAT' } ;
  if strcmp( opt.coord, 'radec' ), LBL_COORD = { 'RA', 'DEC' } ; end
  if strcmp( opt.coord, 'gal' ), LBL_COORD = { 'GALACTIC LON', 'GALACTIC LAT' } ; end
xlabel( [ LBL_COORD{ 1 }, ' (Degrees) ' ], 'FontSize', 15 ) ;
ylabel( [ LBL_COORD{ 2 }, ' (Degrees) ' ], 'FontSize', 15 ) ;
grid
  if ~strcmp( opt.coord, 'ecl' )
  h = colorbar ;
  ylabel(h, { 'ECLIPTIC LATITUDE WFIRST-S', '(Degrees) ' }, 'FontSize', 15) ;
  end

ttl = sprintf( 'MISSION FIRST LIGHT %04d. DAY OF MISSION %3.2f', opt.first_light_date, opt.day_of_mission )  ;
  if ( opt.plot_pm ), ttl = { ttl, sprintf( 'PROPER MOTIONS IN %d YEARS', opt.pm_fact ) } ; end
h_ttl = title( ttl ) ;
set( h_ttl, 'FontSize', 15 ) ;
set(gca,'FontSize',15) ;

% Adding some systems
opt = list_of_targets( opt ) ;
oplot_targets( opt ) ;
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',@myupdatefcn)
box on

% A list of targets for testing purposes
function opt_out = list_of_targets( opt )
opt_out = opt ;

 opt_out.target.name = { 'Tau Ceti', 'Trappist-1', 'Beta Cen', 'Ups And', 'Q01 Eri', 'HD 13931 b',  'Eps Eridani', 'HD 30562', 'GJ 179b', ...
                    'HD 33636 b', 'Pi. Men', '47 UMa', 'Mu Ara','16 Cyg B', 'HR 8734' } ;
 % Value in the figure of merit
 opt_out.target.fom = [ 15, 15, 6, 7, 5, 8, 10, 7, 12, 9, 10.5, 9.5, 6, 5, 6] ;
 % RADEC coordinates
 opt_out.target.radec = { [ 26.02140, -15.93955 ], [ 346.62233, -5.04144 ], [ 210.95586, -60.37304 ], ...
                          [ 24.1993, 41.40546 ], [ 25.62215, -53.74083 ], [ 34.19741, 43.7730 ], [ 53.23269, -09.45826 ], ...
                          [ 72.15160, -5.67404 ], [ 73.02388, 6.47657 ], [ 77.94353, 4.40353 ], [ 84.29122, -80.46912 ], ...
                          [ 164.86656, 40.43026 ], [ 266.03625, -51.83405 ], [ 295.46655, 50.51752 ], [ 344.56475, -02.39538 ] } ;
 % Proper motions (mas/yr in RA, DEC)
 opt_out.target.pm = { [ -1721.05 , 854.16 ], [ 922.1, -471.9 ], [ -33.27, -23.16 ], [ -173.33, -381.80 ], [ 166.32, -106.52 ], ...
                       [ 99.03, -183.19 ], [ -975.17, 19.49 ], [ 311.04, -249.44 ], [ 142.98, -309.39 ], [ 179.69, -138.40 ], ...
                       [ 312.01, 1050.38 ], [ -317.01, 54.64 ], [ -16.85, -190.60 ], [ -135.11, -163.78 ], [ -6.35, -15.80 ]  } ;
 % Parallax (mas!)
 opt_out.target.pax = [ 273.96, 82.58, 8.32, 74.12, 57.36, 22.61, 310.94, 37.85,  81.38, 35.25, 54.60, 71.11, 64.47, 47.14, 50.36] ;
 
% List of targets from Exocat
function opt_out = list_of_targets_exocat( opt )

% Nota Bene: Matlab cannot read binary Excel sheets unless in Windows OS. So, this is ExoCat.xlsb, then saved as xlsx.
[num,txt,raw] = xlsread('in/ExoCat1.xlsx'); %import data from ExoCat
% IDs must be common names or Hipparcos ID #s
% Examples: 'HIP 19849', 'HIP 30920'
opt_out.target.name = lower(opt_out.target.name)
txt = lower(txt);
L = length(opt_out.target.name);
s1 = 'HIP';

for x=1:L
    target=opt_out.target.name(x);
    %find data for HIP ID
    isHIP = strncmpi(s1,target,3);
    if isHIP == 1
        [HIP,HIPID] = strtok(target);
        id_num = str2double(HIPID);
        ind = find(num(:,1)==id_num);
        notfound = isempty(ind);
        if notfound == 1
            error('Error. %s was not found in ExoCat. Try a different identifier.', target)
        end
        [row, col] = ind2sub(size(num),ind);
        ra(x) = num((row),12);
        dec(x) = num((row),13); 
        pmra(x) = num((row),14);
        pmdec(x) = num((row),15);
        opt_out.target.pax(x) = 1000*(1/num((row),18)); 
    end 
       
    %find data for common names
    if isHIP ==0 
        ind = find(txt==target);
        notfound = isempty(ind);
        if notfound == 1
            error('Error. %s was not found in ExoCat. Try a different identifier.', target)
        end
        [row, col] = ind2sub(size(txt),ind);
        ra(x) = num((row)-1,12);
        dec(x) = num((row)-1,13); 
        pmra(x) = num((row)-1,14);
        pmdec(x) = num((row)-1,15);
        opt_out.target.pax(x) = 1000*(1/num((row)-1,18));
    end 
    %compile data into cell pairs
    opt_out.target.fom = [ 15, 10];
    opt_out.target.radec{x} = [ra(x), dec(x)];
    opt_out.target.pm{x} = [pmra(x), pmdec(x)];
end 

% Overplotting some targets
function oplot_targets( opt )

  if numel( opt.target.radec ) ~= numel( opt.target.fom )
  disp( sprintf( 'WARNING: The number of targets is %d, but there are %d values with figure of merit (fom).', numel( opt.target.radec ), numel( opt.target.fom ) ) ) ;
  return
  end
set(gcf,'Units','normalized');
hold all
  for i_trgt = 1 : numel( opt.target.fom )
    % if necessary, change the coordinates
    if ~( strcmp( opt.coord, 'radec' ) )
      if strcmp( opt.coord, 'ecl' )
      [ X_COORD, Y_COORD ] = equat2eclip( opt.target.radec{ i_trgt }( 1 ) / 180 * pi, opt.target.radec{ i_trgt }( 2 ) / 180 * pi ) ; 
      X_COORD = 180 / pi * X_COORD ;
      Y_COORD = 180 / pi * Y_COORD ;
      end
      if strcmp( opt.coord, 'gal' ), [ X_COORD, Y_COORD ] = RA_Dec_to_Gal( opt.target.radec{ i_trgt }( 1 ), opt.target.radec{ i_trgt }( 2 ) ) ; end
    else
    X_COORD = opt.target.radec{ i_trgt }( 1 ) ;
    Y_COORD = opt.target.radec{ i_trgt }( 2 ) ;
    end
  % Longitudinal coordinate beween 0 and 360 degrees
  X_COORD = mod( X_COORD + 720, 360 ) ;
  plot( X_COORD, Y_COORD, 'ok', 'MarkerSize', opt.target.fom( i_trgt ),'MarkerFaceColor','k' ) % red circles
  plot( X_COORD, Y_COORD, 'or', 'MarkerSize', opt.target.fom( i_trgt ) / 1.3,'MarkerFaceColor','r' ) % black inner circles
  % Plotting the proper motion
    if ( opt.plot_pm )
    [ PM_X PM_Y ] = pm_vect( i_trgt, opt ) ;
    % mas/yr to degrees/yr -> fct = 1/1000/3600. Then x(opt.pm_fact) for visualization purposes
    quiver( X_COORD, Y_COORD, PM_X/3.6e6, PM_Y/3.6e6, opt.pm_fact, 'LineWidth', 2 ) ;
    end
  % Small number on top.
    if strcmp( opt.target.label, 'number' )
    text( X_COORD-7, Y_COORD+5, num2str( i_trgt ) ) ;
    end
    if strcmp( opt.target.label, 'name' )
    text( X_COORD-7, Y_COORD+5, opt.target.name( i_trgt ) ) ; 
    end
  end % i_trgt
 
hold off

% Earth year (circular orbit)
function EARTH_PERIOD = earth_period()
EARTH_PERIOD = 2 * pi / sqrt( 6.67e-11 * 1.98855e30 ) * sqrt( ( 149.6e9 )^3 ) / 86400 ; % days
% L2 period: for now, considering WFIRST-S orbit locked with Earth's orbit and tehrefore same period.
% If at some point, an orbit around L2 is considered, at the very least, the places where l2_orbit() is used should be updated
function L2_PERIOD = l2_period()
L2_PERIOD = earth_period();

function plot_full_period_generic( LON_GRID, LAT_GRID, ECL_LAT_GRID, opt )
clf
% Number of intervals
% L2 orbital period
L2_PERIOD = l2_period() ;
n_interval = round( L2_PERIOD / opt.days_interval ) ;

opt_tmp = opt ;
opt_tmp.day_of_mission = 0 ;
  for i_int = 1 : n_interval
  % Creating again the grid of observed coordinates
    if i_int > 1, opt_tmp.alpha_0 = mod( opt_tmp.alpha_0 + ( 360 / n_interval ), 360 ) ; end
  [ LON_GRID, LAT_GRID ECL_LAT_GRID ] = observed_grid( opt_tmp ) ;
  % Some generic plot
  plot_generic( LON_GRID, LAT_GRID, ECL_LAT_GRID, opt_tmp ) ;
  end
% Adding some systems
opt_tmp = list_of_targets( opt_tmp ) ;
oplot_targets( opt_tmp ) ;
  if ~mod( n_interval, 2 ), grid ; end
ttl = sprintf( 'MISSION FIRST LIGHT %04d. TIME LAPSE %d DAYS', opt_tmp.first_light_date, opt_tmp.days_interval ) ;
  if ( opt_tmp.plot_pm ), ttl = { ttl, sprintf( 'PROPER MOTIONS IN %d YEARS', opt_tmp.pm_fact ) } ; end
h_ttl = title( ttl ) ;
set( h_ttl, 'FontSize', 15 ) ;
  if ( opt_tmp.save_image )
    if ~isdir( opt_tmp.path_image ), system( [ 'mkdir -p ' opt_tmp.path_image ] ) ; end
  saveas( gca, [ opt_tmp.path_image 'wfirst-s_' opt_tmp.coord '_multiple_days' ], 'png' ) ;
  end

function setwinsize(winhandle,x,y)
% setwinsize(handle,x,y)
%
% Set the size of a window without moving the top left corner
%
% eg: setwinsize(gcf,1000,500)

set(winhandle,'Units','pixels');
p=get(winhandle,'Position');
top=p(2)+p(4);
p(3)=x;
p(4)=y;
p(2)=top-p(4);
set(winhandle,'Position',p);

% function to derive the tangent vector of the proper motion in different coordinate systems
function [ w_x, w_y ] = pm_vect( i_target, opt )
% If other coordinates are requested
RA_1 = opt.target.radec{ i_target }( 1 ) ;
DEC_1 = opt.target.radec{ i_target }( 2 ) ;
% Coordinates of the proper motion are given in SIMBAD convention:
% pm-ra : mu-ra*cos(dec) (expressed in the ICRS system in mas/yr)
% pm-dec : mu-dec (expressed in the ICRS system in mas/yr)
w_x = opt.target.pm{ i_target }( 1 ) / cos( DEC_1 * pi /180 ) ;
w_y = opt.target.pm{ i_target }( 2 ) ;
% mas per year (1 year)=mas -> deg
RA_2 = RA_1 + w_x / 1e3 / 3600 ;
DEC_2 = DEC_1  + w_y / 1e3 / 3600 ;
  if strcmp( opt.coord, 'ecl' )
  [ X_1, Y_1 ] = equat2eclip( RA_1 / 180 * pi, DEC_1 / 180 * pi ) ;
  [ X_2, Y_2 ] = equat2eclip( RA_2 / 180 * pi, DEC_2 / 180 * pi ) ;
  w_x = ( X_2 - X_1 ) * 180 / pi ;
  w_y = ( Y_2 - Y_1 ) * 180 / pi ;
  % Back to mas
  w_x = w_x * 1e3 * 3600 ;
  w_y = w_y * 1e3 * 3600 ;
  end
  if strcmp( opt.coord, 'gal' ) 
  [ X_1, Y_1 ] = RA_Dec_to_Gal( RA_1, DEC_1 ) ; 
  [ X_2, Y_2 ] = RA_Dec_to_Gal( RA_2, DEC_2 ) ;  
  w_x = X_2 - X_1 ;
  w_y = Y_2 - Y_1 ;
  % Back to mas
  w_x = w_x * 1e3 * 3600 ;
  w_y = w_y * 1e3 * 3600 ;
  end

% Function to get the position of the target at some given time in Equatorial coordinates
% Includes parallax 
function [ ra_3 dec_3 ra_2 dec_2 ra_4 dec_4 ] = get_radec_pos_target( opt, i_trgt )
% RA, DEC of the object
ra_1 = opt.target.radec{ i_trgt }( 1 ) ; 
dec_1 = opt.target.radec{ i_trgt }( 2 ) ; 
% Get proper motion and translate it into RA, DEC 
% SIMBAD convention:
% pm-ra : mu-ra*cos(dec) (expressed in the ICRS system in mas/yr)
% pm-dec : mu-dec (expressed in the ICRS system in mas/yr) 
pm_ra = opt.target.pm{ i_trgt }( 1 ) ;
pm_dec = opt.target.pm{ i_trgt }( 2 ) ;
pm_ra = pm_ra / cos( dec_1 * pi /180 ) ;
% Add the subtended angle from J2000 to the date of observation
% mas per year (1 year)=mas -> deg
% Date of observation in fraction of a year
dt_obs = opt.first_light_date + opt.day_of_mission / earth_period() ;
ra_2 = ra_1 + pm_ra / 1e3 / 3600 * ( dt_obs - 2000 ) ;
dec_2 = dec_1  + pm_dec / 1e3 / 3600  * ( dt_obs - 2000 );
% Parallax
[ ra_3, dec_3 ] = add_parallax( ra_2, dec_2, opt, i_trgt ) ;
% Only parallax (for visual purposes)
[ ra_4 dec_4 ] = add_parallax( ra_1, dec_1, opt, i_trgt ) ;

% Function to get the days of observation for a given target
% Notice that the input is the (RA,DEC) coordinates at the *time of observation*, not J2000. Use the function get_radec_obs()
% Assuming WFIRST-S orbits in a direction of increasing Ecliptic longitude in time
function [ Q_DFF, day_of_observation, n_days_of_observation, n_window, days_list ] = get_day_of_obs_string( RA_OBS, DEC_OBS, opt )
n_window = 0 ;
n_days_of_observation = 0 ;
day_of_observation = '' ;
% The maximum four days of crossing with WFIRST-S field of view
days_list = [ NaN, NaN, NaN, NaN ] ;
% Create the grid of observation in Ecliptic coordinates
opt_2 = opt ;
% Ecliptic
opt_2.coord = 'ecl' ;
[ LON_GRID, LAT_GRID, DUMMY ] = observed_grid( opt_2 ) ;
% Transforming the position of the target
% Needs to be updated to Ecliptical coordinates at the time of observation, not J2000
[ ECL_LON, ECL_LAT ] = equat2eclip( RA_OBS / 180 * pi, DEC_OBS / 180 * pi ) ;
ECL_LON = ECL_LON * 180 / pi ;
% Shifting the longitude grid to have the target at the origin and assuming WFIRST-S orbit in the direction of
% increasing Ecliptic longitude as time passes, like Earth about the Sun.
LON_TRGT = ECL_LON - LON_GRID ;
ECL_LAT = ECL_LAT * 180 / pi ;
% Finding the latitude cut with WFIRST-S observational strategy (precision up to opt.grid degrees)
DFF_LAT = ( ECL_LAT - LAT_GRID ) ;
Q_DFF = find( abs( DFF_LAT( : ) ) < opt_2.grid ) ;
% Whether it is observable
  if numel( Q_DFF ) == 0
  disp( '(WFIRST_S_OBSERVATIONS/GET_DAY_OF_OBSERVATION) The target is not observable with WFIRST-S' ) ;
  day_of_observation = 'UNOBSERVABLE' ;
  return
  end
LON_OBS = sort( LON_TRGT( Q_DFF ) ) ; 
% Finding the boundaries of the observed longitude
D_LON = LON_OBS - circshift( LON_OBS, 1 ) ;
Q_LON = find( abs( D_LON ) > opt.grid * 10 ) ;
N_Q_LON = numel( Q_LON ) ;
  if N_Q_LON > 2
  disp( '(WFIRST_S_OBSERVATIONS/GET_DAY_OF_OBSERVATION) The latitude of the target cuts the WFIRST-S field of view more than twice. Odd. returning.' ) ;
  keyboard
  end
  if N_Q_LON == 0
  disp( 'Target close to tangential observation by Starshade field of view. Define a finer grid with opt.grid to make sure whether it is observable and how much. Now, not considering it like a target.');
  day_of_observation = 'TANGENTIAL';
  return 
  end
  if N_Q_LON == 1
  % Single interval of days to observe the target
  n_window = 1 ;
  % WFIRST-S is orbiting in the direction of increasing Ecliptic longitude, like the Earth
  DAY_1 = LON_OBS( 1 ) / 360 * l2_period() - mod( opt_2.first_light_date, 1 ) * earth_period() ;
  DAY_2 = LON_OBS( end ) / 360 * l2_period() - mod( opt_2.first_light_date, 1 ) * earth_period() ;
  n_days_of_observation = DAY_2 - DAY_1 ;

  % Different possibilities
    if ( DAY_1 < 0 ) && ( DAY_2 < 0 )
    day_of_observation = sprintf( 'The target is observable in the second orbit and further from days %3.2f to %3.2f', DAY_1 + l2_period(), DAY_2 + l2_period ) ;
    end
    if ( DAY_1 < 0 ) && ( DAY_2 >=0 )
    day_of_observation = sprintf( 'The target is observable from days 0 to %3.2f in each orbit. And in the second and further orbits from days %3.2f to %3.2f (end of orbital period)', DAY_2, DAY_1 + l2_period(), l2_period() ) ;
    end
    if ( DAY_1 >= 0 ) && ( DAY_2 >= 0 )
    day_of_observation = sprintf( 'The target is observable from days % to % in each orbit', DAY_1, DAY_2 ) ;
    end
    if ( DAY_1 >= 0 ) && ( DAY_2 < 0 )
    disp( '(WFIRST_S_OBSERVATIONS/GET_DAY_OF_OBSERVATION) This case should not happen. Look at it. Returning' ) ;
    day_of_observation = 'ODD CASE. Look at get_day_of_observation_string.' ;	
    end
  days_list = [ DAY_1, DAY_2, NaN, NaN ] ;
  return
  else
  % Two cuts over the WFIRST-S field of view
  n_window = 2 ;
  DAY_1 = LON_OBS( 1 ) / 360 * l2_period() - mod( opt_2.first_light_date, 1 ) * earth_period() ;
  DAY_2 = LON_OBS( Q_LON( 2 ) - 1 ) / 360 * l2_period() - mod( opt_2.first_light_date, 1 ) * earth_period() ;
  DAY_3 = LON_OBS( Q_LON( 2 ) ) / 360 * l2_period() - mod( opt_2.first_light_date, 1 ) * earth_period() ;
  DAY_4 = LON_OBS( end ) / 360 * l2_period() - mod( opt_2.first_light_date, 1 ) * earth_period() ;
  n_days_of_observation = ( DAY_2 - DAY_1 ) + ( DAY_4 - DAY_3 ) ;
  % Several cases
  % DAY_4 < 0 (all <0)
    if ( DAY_4 < 0 )
    day_of_observation = sprintf( 'The target is observable in the second and further orbits. From days %3.2f to %3.2f and from %3.2f to %3.2f', DAY_1 + l2_period(), ...
    DAY_2 + l2_period(), DAY_3 + l2_period(), DAY_4 + l2_period() ) ;
    end
  % DAY_4 >= 0, but DAY_3 < 0 (DAY_1, DAY_2 must also be <0) 
    if ( ( DAY_4 >= 0 ) && ( DAY_3 < 0 ) )
    day_of_observation = sprintf( 'The target is observable in each orbit from days 0 to %3.2f. And in the second and further orbits from days %3.2f to %3.2f and from %3.2f to %3.2f (end of orbital period)', DAY_4, DAY_1 + l2_period(), DAY_2 + l2_period(), DAY_3 + l2_period(), l2_period() ) ; 
    end
  % DAY_3>=0, but DAY_2<0 (DAY_1 is also<0)
    if ( ( DAY_3 >= 0 ) && ( DAY_2 < 0 ) )
    day_of_observation = sprintf( 'The target is observable in each orbit from days %3.2f to %3.2f. And in the second and further orbits from days %3.2f to %3.2f.', DAY_3, DAY_4, DAY_1 + l2_period(), DAY_2 + l2_period() ) ;
    end
  % DAY_2>=0, but DAY_1<0 (DAY_3 and DAY_4 must be >=0)
    if ( ( DAY_2 >= 0 ) && ( DAY_1 < 0 ) )
    day_of_observation = sprintf( 'The target is observable in each orbit from days 0 to %3.2f and from days %3.2f to %3.2f. And in the second and further orbits from day %3.2f to %3.2f (end of orbital period)', DAY_2, DAY_3, DAY_4, DAY_1 + l2_period(), l2_period() ) ;
    end
  % DAY_1>=0, and DAY_4>0 (all must be >=0)
    if ( ( DAY_1 >= 0 ) && ( DAY_4 >= 0 ) )
    day_of_observation = sprintf( 'The target is observable in each orbit from days %3.2f to %3.2f and from %3.2f to %3.2f', DAY_1, DAY_2, DAY_3, DAY_4 ) ;
    end
  % DAY_1>=0, but DAY_4<0 and DAY_3 >= 0
    if ( ( DAY_1 >= 0 ) && ( DAY_3 >= 0 ) )
    day_of_observation = sprintf( 'The target is observable in each orbit from days %3.2f to %3.2f and from %3.2f to %3.2f (end of orbital period). And in the second and further orbits from 0 to %3.2f', DAY_1, DAY_2, DAY_3, l2_period(), DAY_4 + l2_period() ) ;
    end
  % DAY_1>=0, but DAY_4<0 and DAY_3<0
    if ( ( DAY_1 >= 0 ) && ( DAY_3 < 0 ) && ( DAY_4 < 0 ) )
    day_of_observation = sprintf( 'The target is observable in each orbit from days %3.2f to %3.2f. And from days %3.2f to %3.2f', DAY_1, DAY_2, DAY_3 + l2_period(), DAY_4 + l2_period() ) ;
    end
  days_list = [ DAY_1, DAY_2, DAY_3, DAY_4 ] ;
  return
  end

% Function to get 1/0 if the target is observed/Not observed
% Notice that the input is the (RA,DEC) coordinates at the *time of observation*, not J2000. Use the function get_radec_obs()
% Assuming WFIRST-S orbits in a direction of increasing Ecliptic longitude in time
function target_observed = is_target_observed( RA_OBS, DEC_OBS, opt )
% Create the grid of observation in Ecliptic coordinates
opt_2 = opt ;
% Ecliptic
opt_2.coord = 'ecl' ;
[ LON_GRID, LAT_GRID, DUMMY ] = observed_grid( opt_2 ) ;
% Transforming the position of the target
% Needs to be updated to Ecliptical coordinates at the time of observation, not J2000
[ ECL_LON, ECL_LAT ] = equat2eclip( RA_OBS / 180 * pi, DEC_OBS / 180 * pi ) ;
ECL_LON = ECL_LON * 180 / pi ;
ECL_LAT = ECL_LAT * 180 / pi ;
% Shifting the longitude grid to have the target at the origin and assuming WFIRST-S orbit in the direction of
% increasing Ecliptic longitude as time passes, like Earth about the Sun.
  for i = 1 : numel( ECL_LON )
  LON_TRGT = ECL_LON - LON_GRID( i ) ;
  % Finding the latitude cut with WFIRST-S observational strategy (precision up to opt.grid degrees)
  DFF_LAT = ( ECL_LAT( i ) - LAT_GRID ) ;
  Q_DFF = find( abs( DFF_LAT( : ) ) < opt_2.grid ) ;
  % Whether it is observable
  target_observed( i ) = 0 ;
    if numel( Q_DFF ) ~= 0, target_observed( i ) = 1 ; end
  end

% Click on object to select data point
function txt = myupdatefcn(empt,event_obj)
pos = get(event_obj,'Position');
txt = {['Time: ',num2str(pos(1))],...
       ['Amplitude: ',num2str(pos(2))], ...
       ['Total: ', num2str( pos(1)+pos(2) )]};

% This sub-function is not used. Delete it if it is not recycled in some time (March 2017)
% Plotting individual objects
function plot_individual_target( opt )
n_trgt = numel( opt.target ) ;
setwinsize( gcf, 800, 350 ) ;
  for i_target = 1 : n_trgt
  clf
  hold off
  h_ttl = title( [ 'TARGET: ' opt.target{ i_target } ] ) ;
  set( h_ttl, 'FontSize', 15 ) ;
  subplot(121)
  % Calculating the range of coordinates where the target may be observed:
  % Derive its Ecliptic coordinates. Find the times where the ecliptic longitude will be visible. Derive the time vector.
  % Recompute the coordinates in the chosen system. Plot.
  [ ECL_LON, ECL_LAT ] = equat2eclip( opt.target.radec{ i_target }( 1 ) / 180 * pi, opt.target.radec{ i_target }( 2 ) / 180 * pi ) ;
  ECL_LON = ECL_LON * 180 / pi ;
  ECL_LAT = ECL_LAT * 180 / pi ;
  % Looking at one limiting crossing. Remember, assuming WFIRST_S orbits in the direction of increasing Ecliptic longitude
  [ LON_GRID, LAT_GRID ECL_LAT_GRID ] = observed_grid( opt ) ;
  LON_ARRAY = ( min( [ D_LON_1, D_LON_2 ] ) : opt.grid : max( [ D_LON_1, D_LON_2 ] ) ) ;
  LAT_ARRAY = repmat( ECL_LAT, 1, numel( LON_ARRAY ) ) ;
  LAT_WFIRST_S_UP = repmat( min( ECL_LAT( : ) ), 1, numel( LON_ARRAY ) ) ;
  LAT_WFIRST_S_DOWN = repmat( max( ECL_LAT( : ) ), 1, numel( LON_ARRAY ) ) ;
  % The times
  TIME_ARRAY = LON_ARRAY / 360 * l2_period() ;
  % Negative times get an orbital period added (maybe only for the ticks)
  %Q = find( TIME_ARRAY < 0 ) ;
  %TIME_ARRAY( Q ) = TIME_ARRAY( Q ) + l2_period ;
  % The latitude in other coordinates (Gal->Ecl to RADEC+RADEC->gal)
  subplot(121)
  hold all
  plot( TIME_ARRAY, LAT_ARRAY, 'b', 'LineWidth', 1.5 )
  text(  mean( TIME_ARRAY ), mean( LAT_ARRAY ) + 7, opt.target( i_target ) ) ;
  plot( TIME_ARRAY, LAT_WFIRST_S_UP, 'k', 'LineWidth', 1.5 )
  text( mean( TIME_ARRAY ), min( [ mean( LAT_WFIRST_S_UP ) + 10, 85 ] ), 'MAX LAT WFIRST-S' ) ;
  plot( TIME_ARRAY, LAT_WFIRST_S_DOWN, 'k', 'LineWidth', 1.5 )
  text( mean( TIME_ARRAY ), max( [ mean( LAT_WFIRST_S_DOWN ) - 10, -85 ] ), 'MIN LAT WFIRST-S' ) ;
  ylim( [ -90, 90 ] )
  xlim( [ min( TIME_ARRAY ), max( TIME_ARRAY ) ] ) ;
  ttl = [ 'MISSION FIRST LIGHT: ' num2str( opt.first_light_date ) ] ;
  h_ttl = title( ttl ) ;
  set( h_ttl, 'FontSize', 15 ) ;
  xlabel( 'DAYS AFTER FIRST LIGHT', 'FontSize', 15 ) ;
  ylbl = 'ECLIPTIC LATITUDE (DEGREES) ' ;
  ylabel( ylbl, 'FontSize', 15 ) ;
  grid
  set(gca,'FontSize',15) ;
  % second plot (proper motion and parallax)
  subplot(122)
  hold all
  SZ_SQR = 3 ; % arcsec
  % Camera coordinates: Ecliptical
  opt_camera = opt ;
  opt_camera.coord = 'ecl' ;
  [ PM_X PM_Y ] = pm_vect( i_target, opt_camera ) ;
  % mas after a year in camera pixels
  fct_mission = l2_period() / earth_period()  ;
  % Remember camera's pixel is in mas
  quiver( 0, 0, PM_X/opt.camera_pixel, PM_Y/opt.camera_pixel, fct_mission, 'LineWidth', 2 ) ;
  xlim( [ -SZ_SQR*1e3/opt.camera_pixel, SZ_SQR*1e3/opt.camera_pixel ] ) ;
  ylim( [ -SZ_SQR*1e3/opt.camera_pixel, SZ_SQR*1e3/opt.camera_pixel ] ) ;
  xlabel( 'CAMERA PIXELS', 'FontSize', 15 ) ;
  ylabel( 'CAMERA PIXELS', 'FontSize', 15 ) ;
  ttl = { 'TARGET MOTION AFTER ONE ORBIT', sprintf( 'FOV=%dx%d arcsec^2. 1 PIX=%d mas', SZ_SQR, SZ_SQR, opt.camera_pixel ) } ;
  h_ttl = title( ttl ) ;
  set( h_ttl, 'FontSize', 15 ) ;
  grid
    if ( opt.save_image )
      if ~isdir( opt.path_image ), system( [ 'mkdir -p ' opt.path_image ] ) ; end
    saveas( gca, [ opt.path_image 'wfirst-s_' strrep( opt.target{ i_target }, ' ', '_' ) '_observability' ], 'png' ) ;
    end
  end % i_target

% Deriving the parallax in RADEC
% Simplifications:
% 1: WFIRST-S Orbit and Earth orbit same plane (Ecliptic latitude WFIRST-S = 0)
% 2: Distance WFIRST-S to Sun same as Earth to Sun: 1% error, parallax error less than a hundredth of arcsec (negligible)
% 3: WFIRST-S orbit locked with Earth (exact L2) in a circular (0.0167 -> ~another 1, negligible) orbit of 1 year period (this may not be good enough, since the influence of other planets on Earth's orbit might shift the position of stars by a similar amount as the parallax from J2000 to 2028
% The heliocentric position vector of a star is the same. From WFIRST-S perspective, the star moves along its orbit -> parallax (pax)
function [ ra_wpax_obs dec_wpax_obs ] = add_parallax( ra_obs, dec_obs, opt, i_trgt )
% Orbital phase (recall #3). NB: subtracting J2000 is just to remember that the reference is J2000, with the circular orbit assumption and constant angular velocity, it just adds a multiple of 2 pi.
phi_orb = 2 * pi * mod( ( opt.first_light_date - 2000 + opt.day_of_mission / earth_period() ), 1 ) ; 
% Vector position of WFIRST-S with respect the Sun amd J2000
d_ES = 1.496e11 ; % m
x_wfirst_s = d_ES * cos( phi_orb ) ;
y_wfirst_s = d_ES * sin( phi_orb ) ;
z_wfirst_s = 0 ; % #1
% Vector position of the star from WFIRST-S
% Converting parallax to distance (opt.target.pax is in mas)
d_str = 1000 / opt.target.pax( i_trgt ) * 3.0857e16 ; % 1 pc = 3.0857e16 m
% Converting ra_obs, dec_obs to Ecliptic coordinates
[ ecl_lon_obs, ecl_lat_obs ] = equat2eclip( ra_obs * pi / 180, dec_obs * pi / 180 ) ;
% for orbital phase equal to 0, WFIRST-S aligned with the reference time that defines J2000
x_str_wfirst_s = d_str * cos( ecl_lat_obs ) * cos( ecl_lon_obs ) - x_wfirst_s ;
y_str_wfirst_s = d_str * cos( ecl_lat_obs ) * sin( ecl_lon_obs ) - y_wfirst_s ;
z_str_wfirst_s = d_str * sin( ecl_lat_obs ) ;
% Normalized position vector from WFIRST-S:
d_wfirst_s = sqrt( x_str_wfirst_s * x_str_wfirst_s + y_str_wfirst_s * y_str_wfirst_s + z_str_wfirst_s * z_str_wfirst_s ) ;
x_nrm = x_str_wfirst_s / d_wfirst_s ;
y_nrm = y_str_wfirst_s / d_wfirst_s ;
z_nrm = z_str_wfirst_s / d_wfirst_s ;
% Deriving the new Ecliptic coordinates
ecl_lat_wpax_obs = asin( z_nrm ) ;
% atan2(Y,X)
ecl_lon_wpax_obs = atan2( y_nrm, x_nrm ) ;
% Ecliptic to Equatorial
[ ra_wpax_obs dec_wpax_obs ] = eclip2equat( ecl_lon_wpax_obs, ecl_lat_wpax_obs ) ;
% Back to degrees
ra_wpax_obs = ra_wpax_obs * 180 / pi ;
dec_wpax_obs = dec_wpax_obs * 180 / pi ;

% Target's motion on the sky for a given period of time
function [ ra_obs_intrp dec_obs_intrp ra_obs_wfirst_s dec_obs_wfirst_s ra_pm_intrp dec_pm_intrp dy_obs_wfirst_s ra_pax_intrp dec_pax_intrp ] = get_radec_obs_arry( opt, i_trgt ) 
ra_obs_intrp = NaN ;
dec_obs_intrp = NaN ;
ra_obs_wfirst_s = NaN ;
dec_obs_wfirst_s = NaN ;
ra_pm = NaN ;
dec_pm = NaN ;
dy_obs_wfirst_s = NaN ;
% period of time to be plotted
dy_1 = 0 ;
dy_2 = earth_period() * ( opt.end_light_date - opt.first_light_date - opt.day_of_mission / earth_period() ) ;
% time lapse between evaluations (proper motion and parallax are small)
% fraction of a year (format to compute the orbital phase in get_radec_pos_target)
dlt_dy = 10 ; % days
dy_arry = dy_1 : dlt_dy : dy_2 ;
n_dy = numel( dy_arry ) ;

% Getting the motion in arrays
opt2= opt ;
  for i_dy = 1 : n_dy
  opt2.day_of_mission = dy_arry( i_dy ) ;
  [ ra_obs( i_dy ) dec_obs( i_dy ) ra_pm( i_dy ) dec_pm( i_dy ) ra_pax( i_dy ) dec_pax( i_dy ) ] = get_radec_pos_target( opt2, i_trgt ) ;
  end
% Interpolating the results to a finer time grid
dlt_dy_intrp = 0.1 ; % days
dy_arry_intrp = dy_1 : dlt_dy_intrp : dy_2 ;
ra_obs_intrp = spline( dy_arry, ra_obs, dy_arry_intrp ) ;
dec_obs_intrp = spline( dy_arry, dec_obs, dy_arry_intrp ) ;
ra_pm_intrp = spline( dy_arry, ra_pm, dy_arry_intrp ) ;
dec_pm_intrp = spline( dy_arry, dec_pm, dy_arry_intrp ) ;
ra_pax_intrp = spline( dy_arry, ra_pax, dy_arry_intrp ) ;
dec_pax_intrp = spline( dy_arry, dec_pax, dy_arry_intrp ) ;

% Finding the days the target is observable
% Looping over the four possible crossings with the WFIRST-S field of view during the period of time between the first and end light.
clear dy_crss
i_elem = 1 ;
i_dy_wfirst_s = 1 ;
  for i = 0 : ceil( opt.end_light_date - opt.first_light_date ) + 1
    for j = 1 : 2
    dy_tmp_1 = opt.target.days_list{ i_trgt }( 2 * j - 1 ) ;
    dy_tmp_2 = opt.target.days_list{ i_trgt }( 2 * j ) ;
      % if the first day is NaN, the corresponding second day is also a NaN: skip
      if isnan( dy_tmp_1 ), continue ; end
    arry_tmp = abs( dy_arry_intrp - dy_tmp_1 - earth_period() * i ) ;
    % If there is a crossing, it has to be within the time resolution
    q_tmp_1 = find( arry_tmp <= 0.5 * dlt_dy_intrp ) ;
    arry_tmp = abs( dy_arry_intrp - dy_tmp_2 - earth_period() * i ) ;
    q_tmp_2 = find( arry_tmp <= 0.5 * dlt_dy_intrp ) ;
      % If there are no crossings, continue
      if ( numel( q_tmp_2 ) == 0 ) && ( numel( q_tmp_1 ) == 0 ), continue ; end
      % If the second day does not cross the time interval, either the second day is beond the array limits, or day_1 neither can cross
      if ( numel( q_tmp_2 ) == 0 ) && ( numel( q_tmp_1 ) > 0 )
        if ( dy_tmp_2 + earth_period() * i > dy_arry_intrp( end ) )
        q_tmp_2 = [ numel( dy_arry_intrp ) ] ; 
        else 
        continue 
        end
      end
      % If the second day crosses the time interval, still the first might have no crossing: set the first crossing to time=0
      if numel( q_tmp_1 ) == 0, q_tmp_1 = [ 1 ] ; end
    % feed the observed arrays. 
    ra_obs_wfirst_s( i_elem : i_elem + q_tmp_2( 1 ) - q_tmp_1( 1 ) ) = ra_obs_intrp( q_tmp_1( 1 ) : q_tmp_2( 1 ) ) ;
    dec_obs_wfirst_s( i_elem : i_elem + q_tmp_2( 1 ) - q_tmp_1( 1 ) ) = dec_obs_intrp( q_tmp_1( 1 ) : q_tmp_2( 1 ) ) ;  
    i_elem = i_elem + q_tmp_2( 1 ) - q_tmp_1( 1 ) + 1 ;
    dy_obs_wfirst_s( i_dy_wfirst_s ) = dy_arry_intrp( q_tmp_1( 1 ) ) ;
    dy_obs_wfirst_s( i_dy_wfirst_s + 1 ) = dy_arry_intrp( q_tmp_2( 1 ) ) ;
    i_dy_wfirst_s = i_dy_wfirst_s + 2 ;

    end
  end

% Plot of a single target's position in a period of time
function plot_single_target_in_time( opt, i_trgt )
figure ; clf ; setwinsize( gcf, 500, 550 ) ;
[ ra_trgt dec_trgt ra_obs_wfirst_s dec_obs_wfirst_s ra_trgt_pm dec_trgt_pm dy_obs_wfirst_s ra_trgt_pax dec_trgt_pax ] = get_radec_obs_arry( opt, i_trgt ) ;
% Conversion from degrees to mas
deg2mas = 60 * 60 * 1e3 ;
avg_ra = mean( ra_trgt ) ;
avg_dec = mean( dec_trgt ) ;
avg_ra_pm = mean( ra_trgt_pm ) ;
avg_dec_pm = mean( dec_trgt_pm ) ;
avg_ra_pax = mean( ra_trgt_pax ) ;
avg_dec_pax = mean( dec_trgt_pax ) ;
alpha = deg2mas * ( ra_trgt - avg_ra ) ;
delta = deg2mas * ( dec_trgt - avg_dec ) ;
alpha_pm = deg2mas * ( ra_trgt_pm - avg_ra_pm ) ;
delta_pm = deg2mas * ( dec_trgt_pm - avg_dec_pm ) ;
alpha_obs = deg2mas * ( ra_obs_wfirst_s - avg_ra ) ;
delta_obs = deg2mas * ( dec_obs_wfirst_s - avg_dec ) ;
alpha_pax = deg2mas * ( ra_trgt_pax - avg_ra_pax ) ;
delta_pax = deg2mas * ( dec_trgt_pax - avg_dec_pax ) ;

% PS: plotting pax only is for visual purposes only
  if ~( opt.pax_only )
  plot( alpha, delta, 'k' ) ;
  else
  plot( alpha_pax, delta_pax, 'k' ) ;
  end
hold all
plot( alpha_pm, delta_pm, 'k--' ) ;
  if ~( opt.pax_only ), h_obs = plot( alpha_obs, delta_obs, '.', 'Color',  'Blue' ) ; end
  if isnan( alpha_obs ), text( 0, max( delta_pm ), 'UNOBSERVED' ) ; end
grid
% Keeping the same scale for the option with parallax only (to make it easier to visualize the relative effect)
xlim( [ min( alpha ) * 1.15, max( alpha ) * 1.15 ] ) ;
ylim( [ min( delta ) * 1.15, max( delta ) * 1.15 ] ) ;
title( sprintf( '%s (%3.2f days/year)', regexprep(opt.target.name{ i_trgt },'(\<[a-z])','${upper($1)}'), opt.target.n_days_obs( i_trgt ) ) ) ;
ylabel( 'mas', 'FontSize', 15 ) ;
xlabel( 'mas', 'FontSize', 15 ) ;
n_dy = numel( dy_obs_wfirst_s ) ;
% Is an even number (pairs)
dy_obs_wfirst_s_str{ 1 } = 'Days observed:' ;
  for i_dy = 1 : n_dy / 2
  dy_obs_wfirst_s_str{ i_dy + 1 } = [ num2str( dy_obs_wfirst_s( 2 * i_dy - 1 ) ) ' to ' num2str( dy_obs_wfirst_s( 2 * i_dy ) ) ] ;
  end
text( 0, max( delta )/1.2, dy_obs_wfirst_s_str, 'FontSize', 14 ) ;
ttl_1 = sprintf( 'WFIRST-S. FROM: %3.0f, DAY %3.2f TO %3.0f, DAY %3.2f ', floor( opt.first_light_date ), mod( opt.first_light_date, 1 ) * earth_period() + opt.day_of_mission, floor( opt.end_light_date ), mod( opt.end_light_date, 1 ) * earth_period() ) ;
  if ~( opt.pax_only ), ttl_2 = 'Black: total motion; {\color{blue}Blue}: observable by WFIRST-S; Dotted: proper motion' ; else, ttl_2 = 'Black: parallax motion; Dotted: proper motion' ; end
suptitle( { ttl_1, ttl_2 } ) ;

% Plot of targets' positions in a period of time 
function plot_targets_in_time( opt )
figure( 2 ) ; clf ; setwinsize( gcf, 1250, 650 ) ;
  for i_trgt = 1 : numel( opt.target.radec )
  [ ra_trgt dec_trgt ra_obs_wfirst_s dec_obs_wfirst_s ra_trgt_pm dec_trgt_pm dy_obs_wfirst_s ra_trgt_pax dec_trgt_pax ] = get_radec_obs_arry( opt, i_trgt ) ;
  subplot( 3, 5, i_trgt )
   % Conversion from degrees to mas
  deg2mas = 60 * 60 * 1e3 ;
  avg_ra = mean( ra_trgt ) ;
  avg_dec = mean( dec_trgt ) ;
  avg_ra_pm = mean( ra_trgt_pm ) ;
  avg_dec_pm = mean( dec_trgt_pm ) ;
  avg_ra_pax = mean( ra_trgt_pax ) ;
  avg_dec_pax = mean( dec_trgt_pax ) ;
  alpha = deg2mas * ( ra_trgt - avg_ra ) ;
  delta = deg2mas * ( dec_trgt - avg_dec ) ;
  alpha_pm = deg2mas * ( ra_trgt_pm - avg_ra_pm ) ;
  delta_pm = deg2mas * ( dec_trgt_pm - avg_dec_pm ) ;
  alpha_obs = deg2mas * ( ra_obs_wfirst_s - avg_ra ) ;
  delta_obs = deg2mas * ( dec_obs_wfirst_s - avg_dec ) ;
  alpha_pax = deg2mas * ( ra_trgt_pax - avg_ra_pax ) ;
  delta_pax = deg2mas * ( dec_trgt_pax - avg_dec_pax ) ;
  % PS: plotting pax only is for visual purposes only
    if ~( opt.pax_only )
    plot( alpha, delta, 'k' ) ;
    else
    plot( alpha_pax, delta_pax, 'k' ) ;
    end
  hold all
  plot( alpha_pm, delta_pm, 'k--' ) ;
    if ~( opt.pax_only ), h_obs = plot( alpha_obs, delta_obs, '.', 'Color',  'Blue' ) ; end
    if isnan( alpha_obs ), text( 0, max( delta_pm ), 'UNOBSERVED' ) ; end
  grid
  % Keeping the same scale for the option with parallax only (to make it easier to visualize the relative effect)
  xlim( [ min( alpha ) * 1.1, max( alpha ) * 1.1 ] ) ;
  ylim( [ min( delta ) * 1.1, max( delta ) * 1.1 ] ) ;
  title( sprintf( '%s (%3.2f days/year)', regexprep(opt.target.name{ i_trgt },'(\<[a-z])','${upper($1)}'), opt.target.n_days_obs( i_trgt ) ) ) ;
    if mod( i_trgt - 1, 5 ) == 0, ylabel( 'mas', 'FontSize', 15 ) ; end
    if ( i_trgt / 5 > 2 ), xlabel( 'mas', 'FontSize', 15 ) ; end
  % One legend as an example
    if ~( opt.pax_only ) && ( i_trgt <= 1 )
    n_dy = numel( dy_obs_wfirst_s ) ;
    % Is an even number (pairs)
    dy_obs_wfirst_s_str{ 1 } = 'Days observed:' ;
      for i_dy = 1 : n_dy / 2
      dy_obs_wfirst_s_str{ i_dy + 1 } = [ num2str( dy_obs_wfirst_s( 2 * i_dy - 1 ) ) ' to ' num2str( dy_obs_wfirst_s( 2 * i_dy ) ) ] ; 
      end
    text( 0, max( delta )/1.8, dy_obs_wfirst_s_str ) ;
    end
  end
ttl_1 = sprintf( 'WFIRST-S. FROM: %3.0f, DAY %3.2f TO %3.0f, DAY %3.2f ', floor( opt.first_light_date ), mod( opt.first_light_date, 1 ) * earth_period() + opt.day_of_mission, floor( opt.end_light_date ), mod( opt.end_light_date, 1 ) * earth_period() ) ;
  if ~( opt.pax_only ), ttl_2 = 'Black: total motion; {\color{blue}Blue}: observable by WFIRST-S; Dotted: proper motion;' ; else, ttl_2 = 'Black: parallax motion; Dotted: proper motion;' ; end
 suptitle( { ttl_1, ttl_2 } ) ;

% Function to overplot the trajectory of a target over an image form a catalog
% Needs communication with Aladin Lite, TBD

function plot_target_on_catalog( opt, i_trgt )
figure ; setwinsize( gcf, 475, 350 ) ;
% Notice this is a place holder for now. The script should get the image from a catalog on the web, centered at the coordinates of the object
% at the first light of WFIRST-S, including identified background objects
disp( 'Not sure whether within one arcsec an overplot will be indicative. The halo covers larger areas.' )
% From: http://www.peteryu.ca/tutorials/matlab/plot_over_image_background
% replace with an image of your choice
img = imread('img/Tau_ceti.png');
% set the range of the axes
% The image will be stretched to this.
min_x = 0;
max_x = 8;
min_y = 0;
max_y = 6;
% Flip the image upside down before showing it
imagesc([min_x max_x], [min_y max_y], flipud(img(:,:,1)));
% NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
hold on;
% Test
% make data to plot - just a line.
x = min_x:max_x;
y = (6/8)*x;
plot(x,y,'k-*','linewidth',1.5);
% set the y-axis back to normal.
set(gca,'ydir','normal');

% Function to get the times an object could be disambiguated from the background
function [bad_window, max_delta_L, N_sigma_obs, N_sigma_per_day_obs] = get_sigma_obs(opt,i_trgt)

lambda = opt.lambda/10^9;
idx_snr = (opt.snr_time *10) + 1; %1 entry in array = .1 days, need 11 indicies to get to snr
numel(opt.target.radec);
  [ ra_trgt dec_trgt ra_obs_wfirst_s dec_obs_wfirst_s ra_trgt_pm dec_trgt_pm dy_obs_wfirst_s ra_trgt_pax dec_trgt_pax ] = get_radec_obs_arry( opt, i_trgt ) ;

deg2mas = 60 * 60 * 1e3 ;
avg_ra = mean( ra_trgt ) ;
avg_dec = mean( dec_trgt ) ;
alpha_obs = deg2mas * ( ra_obs_wfirst_s - avg_ra ) ;
delta_obs = deg2mas * ( dec_obs_wfirst_s - avg_dec ) ;
delta_ra = circshift(alpha_obs, 1)-alpha_obs;
delta_dec = circshift(delta_obs, 1)-delta_obs;
dec_obs_wfirst_s = dec_obs_wfirst_s * pi/180;
%skip jumps from last index to first index
delta_ra(1) = 0;
delta_ra(end) = 0;
delta_dec(1) = 0;
delta_dec(end) = 0;

%special case in which target is almost perfectly tangental to FOV
if numel(alpha_obs) == 1 || numel(delta_obs) == 1
    N_sigma_obs = -1;
    N_sigma_per_day_obs = -1;
    max_delta_L = 0;
    bad_window = 0;
    
    return
end 
r=sort(abs(delta_ra),'descend');
t=sort(abs(delta_dec),'descend');
idx= find(abs(delta_ra)>r(6)+.01);
idx_2 = find(abs(delta_dec)>t(6)+.01);

if isempty(idx)
    N_sigma_obs = 0;
    N_sigma_per_day_obs = 0;
    max_delta_L = 0;
    bad_window = 0;
    return
end 

if (numel(idx) > numel(idx_2))
    idx_2 = idx;
elseif numel(idx_2) > numel(idx)
    idx = idx_2;
end 

if (numel(idx)>5) || (numel(idx_2)>5)
    msg = no;
    error(msg)
end

if idx(end) ~= numel(delta_ra)-1
    idx(end+1) = numel(delta_ra) - 1;
end
n_idx = numel(idx);
idx_science = idx+idx_snr;
i_o = idx_snr;

all_delta_L = sqrt(((delta_ra.*cos(dec_obs_wfirst_s)).^2) + ((delta_dec).^2));
for i=1:(n_idx)
    if idx(i) <= idx_snr %case where very small obs window at beginning
        L(i) = 0;
        diff_1 = 0;
        diff_2 = 0;
    continue
    elseif i == n_idx
        if idx(i)-idx(i-1) < idx_snr %case where very small obs window at end
        L(i) = 0;
        diff_1 = 0;
        diff_2 = 0;
        continue 
        end
        idx(i) = idx(i) +1;
    else 
        if i~=1 && idx(i-1) == idx_snr;
            i_o = i_o+1;
        end 
    end 
    delta_L = sqrt(((delta_ra(i_o:idx(i)-1).*cos(dec_obs_wfirst_s(i_o:idx(i)-1))).^2) + ((delta_dec(i_o:idx(i)-1)).^2));
    L(i) = sum(delta_L);
    max_delta_L(i) = max(delta_L);
    if max_delta_L(i) > mean(delta_L) + 2
        max_delta_L(i) = max(delta_L(delta_L~=max(delta_L)));
    end 
    idx_max(i) = find(all_delta_L == max_delta_L(i));
    if i == 1
        diff_1(i) = abs(idx_max(i) - idx_snr);
    else
        diff_1(i)=abs(idx_max(i) - idx_science(i-1));
    end 
    diff_2(i) = abs(idx(i)-1 - idx_max(i));
    i_o=idx(i)+idx_snr;
    end  
%determine if obs window is long enough
days_of_science = (diff_1 + diff_2)/10;
bad_window = days_of_science > opt.t_obs;
FWHM = (lambda/opt.diam) * 206264806.247;
sigma = 5* FWHM / opt.snr; % distance to disambiguate centroid from background: (Agreed as of 07/13/17 to consider 1 full PSF when SNR=5)

N_sigma_obs = L/sigma;%amount of sigmas per observation window
dy_obs_wfirst_s;
difference = abs(circshift(dy_obs_wfirst_s,1) - dy_obs_wfirst_s);
difference(1) = 0;
dy_period_obs = difference(2:2:end);
N_sigma_per_day_obs = sum(N_sigma_obs)/(sum(dy_period_obs)-(numel(dy_period_obs)+1));

% Collecting all the information in Exocat
function opt_out = all_targets(opt)
opt_out = opt ;
data = xlsread('in/ExoCat1.xlsx'); %import data from ExoCat
HIPID_all = xlsread('in/ExoCat1.xlsx', 'A:A'); 
L=length(HIPID_all);

%put data into arrays
for x=1:L 
    ra(x) = data((x),12);
    dec(x) = data((x),13);
    pmra(x) = data((x),14);
    pmdec(x) = data((x),15);
    opt_out.target.pax(x) = 1000*(1/data((x),18)); %finding parallax from distance
    opt_out.target.radec{x} = [ra(x), dec(x)];
    opt_out.target.pm{x} = [pmra(x), pmdec(x)];
end 

% Creating a figure of merit based on background disambiguation run for all targets in ExoCat
function get_back_disambiguation_exocat( opt )

opt = all_targets(opt);

  for i_trgt = 1 : numel( opt.target.radec )
  [ opt.target.radec_obs{ i_trgt }( 1 ), opt.target.radec_obs{ i_trgt }( 2 ) ] = get_radec_pos_target( opt, i_trgt ) ;
  end

  for i_trgt = 1 : numel( opt.target.radec )
  [ Q_DFF opt.target.day_obs_string{ i_trgt } opt.target.n_days_obs( i_trgt ) opt.target.n_window( i_trgt ) opt.target.days_list{ i_trgt } ] = ...
  get_day_of_obs_string( opt.target.radec_obs{ i_trgt }( 1 ), opt.target.radec_obs{ i_trgt }( 2 ), opt ) ;

  [bad_window_tmp, max_delta_L_tmp, N_sigma_obs_tmp, N_sigma_per_day_obs_tmp] = get_sigma_obs(opt,i_trgt);
  bad_window(i_trgt, 1:numel(bad_window_tmp)) = bad_window_tmp;
  max_delta_L(i_trgt) = max(max_delta_L_tmp);
  N_sigma_obs( i_trgt, 1:numel(N_sigma_obs_tmp) ) = N_sigma_obs_tmp;
  N_sigma_per_day_obs(i_trgt,1:numel(N_sigma_per_day_obs_tmp) )= N_sigma_per_day_obs_tmp;
  end

obs_data = N_sigma_per_day_obs(N_sigma_per_day_obs > 0);
days_for_N_sigma = 1./obs_data;
N_dummy = N_sigma_per_day_obs;
num = [1:numel(N_dummy)];
max_delta_L = max_delta_L';
idx_no = find(all(bad_window==0,2) == 1);
N_dummy(idx_no) = [];
max_delta_L(idx_no) = [];
num(idx_no) = [];

f = find(N_dummy >= .2);
max_delta_L_new = max_delta_L(f);
num = num(f);
max_delta_L_new=sort(max_delta_L_new,'descend');

last = horzcat((num+1)', max_delta_L_new);

%create histograms of data
figure
ratio = numel(obs_data) / numel(N_sigma_per_day_obs)
s = std(obs_data);
edges = [min(obs_data) : s/8 : max(obs_data)+s/8];
h = histogram(obs_data,edges)
xlim([min(obs_data(:)), max(obs_data(:))])
ylim = ([-10,400]);
[N,edges] = histcounts(obs_data,edges)
title('ExoCat Histogram');
xlabel('N Sigma per Day Observed');
ylabel('Frequency');

days_for_N_sigma = 1./obs_data;
figure
s_2 = std(days_for_N_sigma);
edges = [min(days_for_N_sigma):s_2/4:max(days_for_N_sigma)+s/4];
h4 = histogram(days_for_N_sigma,edges)
title('ExoCat Histogram');
xlabel('Days of Observation to N_(sigma)');
ylabel('Frequency');

figure
subplot(2,1,1);
s = std(obs_data);
edges = [min(obs_data) : s/8 : max(obs_data)+s/8];
a = prctile(obs_data, 0);
b = prctile(obs_data, 90);
h2 = histogram(obs_data,edges)
xlim([a,b]);
title('Adjusted ExoCat Histogram');
xlabel('N_(sigma) per Day Observed');
ylabel('Frequency');

subplot(2,1,2);
cla reset
s = std(obs_data);
edges = [min(obs_data) : s/8 : max(obs_data)+s/8];
h3 = histogram(obs_data, edges)
set(gca,'YScale','log')
title('Adjusted ExoCat Histogram');
xlabel('N_(sigma) per Day Observed');
ylabel('Frequency');

% Creating a histogram of the background disambiguation
function histogram_back_disambiguation( opt )
  for i_trgt = 1 : numel( opt.target.radec )
  [bad_window_tmp, max_delta_L_tmp, N_sigma_obs_tmp, N_sigma_per_day_obs_tmp] = get_sigma_obs(opt,i_trgt);
  bad_window(i_trgt, 1:numel(bad_window_tmp)) = bad_window_tmp;
  max_delta_L(i_trgt) = max(max_delta_L_tmp);
  N_sigma_obs( i_trgt, 1:numel(N_sigma_obs_tmp) ) = N_sigma_obs_tmp;
  N_sigma_per_day_obs(i_trgt,1:numel(N_sigma_per_day_obs_tmp) )= N_sigma_per_day_obs_tmp;
  end
idx_no = find(all(bad_window==0,2) == 1);
max_delta_L = max_delta_L'
max_delta_L(idx_no) = []
numel(max_delta_L)
max_delta_L=sort(max_delta_L,'descend')

figure
%make histogram
N_sigma_per_day_obs;
s = std(N_sigma_per_day_obs);
edges = [min(N_sigma_per_day_obs):s/4:max(N_sigma_per_day_obs)+s/4];
h = histogram(N_sigma_per_day_obs, edges)
xlim([min(N_sigma_per_day_obs), max(N_sigma_per_day_obs)]);
[N,edges] = histcounts(N_sigma_per_day_obs,edges)

days_for_N_sigma = 1./N_sigma_per_day_obs

figure
s = std(days_for_N_sigma);
edges = [min(days_for_N_sigma(:)):s/3:max(days_for_N_sigma(:))+s/4];
h4 = histogram(days_for_N_sigma,edges)
[N,edges] = histcounts(days_for_N_sigma,edges)
%make table of data
T = array2table(N_sigma_obs);
Ti = array2table((opt.target.name).');
 Tnew = [Ti,T];
 Tnew.Properties.VariableNames{1} = 'Target'
 Tnew;


