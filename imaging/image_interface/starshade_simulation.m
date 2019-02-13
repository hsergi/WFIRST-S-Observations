% NB: include imagescnan in the distribution

function starshade_simulation( opt_in )

% Simple interface to test some simulation setup
% opt has the ingredients for the simulation
% opt.r_star position of the star in mas, usually on-axis (=0)
% opt.psi_star angle position of the star with respect the on-axis in degrees
% opt.r_planet is an array of planet orbital radius in mas
% opt.phase_planet is an array with the orbital phases of the planets in degrees
% opt.contrast_planet is an array of contrast values for each planet
% History:
% 05/01/17: first version. Sergi Hildebrandt (JPL/Caltech)
% 05/14/17: adapted to use opt. Sergi Hildebrandt (JPL/Caltech)

% get default values for the options
opt = get_default_options( opt_in ) ;

  if isfield( opt, 'developer' )
  units_image
  else
  units
  end

n_plnt = numel( opt.r_planet_mas ) ;
disp( sprintf( '(starshade_simulation) Simulating a star and %i planets', n_plnt ) )

% Simulation
%% 1) Star
opt_str = opt ;
opt_str.r_source_mas = opt.r_star_mas ;
opt_str.psi_source_deg = opt.psi_star_deg ;
[ efDefectImg_str lambdaIn  vecPetalArray ] = makeStarshadeImage( opt_str ) ;
IntDefectImg_sim = abs( efDefectImg_str ).^2 ; % This is correct, it does compute the intensity for each wavelength

%% 2) Accumulating the results for each planet
opt_plnt = opt ;
  for i_plnt = 1 : n_plnt
  opt_plnt.r_source_mas = opt.r_planet_mas( i_plnt ) ;
  opt_plnt.psi_source_deg = opt.phase_planet_deg( i_plnt ) ;
  [ efDefectImg_tmp ] = makeStarshadeImage( opt_plnt ) ;
    if ( 0 )
    % A figure
    % Non-linear transformation of the original image
    colormap( bone )
    i_lmbd = 1 ;
    int_img = squeeze( IntDefectImg_sim( :, :, i_lmbd ) + abs( efDefectImg_tmp( :, :, i_lmbd) ).^2 ) ;
    fct_scl = 0. ;
    nw_img = log10( int_img - fct_scl * min( int_img( : ) ) ) ;
    %nw_img = int_img ;
    set( 0, 'defaultlinelinewidth', 1.0 ) ;
    set( 0, 'DefaultAxesFontSize', 14 ) ;
    figure( 1 )
    clf ;
    x_ax = ( -200 : 200 ) * 5 ;
    y_ax = ( -200 : 200 ) * 5 ;
    x_ax = x_ax( 1 : end - 1 ) ;
    y_ax = y_ax( 1 : end - 1 ) ;
    setwinsize( gcf, 600, 550 )
    imagescnan( x_ax, y_ax, nw_img, [ -5, 0 ] ) ;
    xlabel( 'mas' )
    ylabel( 'mas' )
    h = colorbar ;
    ylabel( h, 'Log10(Intensity)', 'FontSize', 16 )
    suptitle( sprintf( 'Starshade simulation at %3.0f nm', lambdaIn( i_lmbd ) / nm ) )
    saveas( gcf, sprintf( 'fig_dev/img_%03i', i_plnt ), 'jpeg' )
    setwinsize( gcf, 600, 420 )
    imagescnan( x_ax, y_ax, nw_img, [ -5, 0 ] ) ;
    xlabel( 'mas' )
    ylabel( 'mas' )
    h = colorbar ;
    ylabel( h, 'Log10(Intensity)', 'FontSize', 16 )
    suptitle( sprintf( 'Starshade simulation at %3.0f nm', lambdaIn( i_lmbd ) / nm ) )
    xlim( [ -400, 400 ] ) ; ylim( [ -500, 500 ] ) ;
    saveas( gcf, sprintf( 'fig_dev/img_zoomed_%03i', i_plnt ), 'jpeg' )
dbstop if error
%make_an_error
    end
  IntDefectImg_sim = IntDefectImg_sim + opt.contrast_planet( i_plnt ) * abs( efDefectImg_tmp ).^2 ;
  end

% Temporary ... FIXME
saveFilename = 'starshade_star_planet_sim' ; 
  if opt.save == 1
  pth_fl_sv = [ opt.save_path '/' saveFilename '.mat' ] ;
  save( pth_fl_sv, 'IntDefectImg_sim', 'lambdaIn', 'vecPetalArray')
  disp( sprintf( '(starshade_simulation) Output array of the simulation stored in: %s', pth_fl_sv ) )
  end

% A figure
% Non-linear transformation of the original image
colormap( bone )
i_lmbd = 1 ;
abs_img = squeeze( IntDefectImg_sim( :, :, 1 ) ) ;
fct_scl = 0. ;
nw_img = log10( abs_img - fct_scl * min( abs_img( : ) ) ) ;
%nw_img = abs_img ;
set(0,'defaultlinelinewidth',1.0);
set(0,'DefaultAxesFontSize',14);
figure( 1 )
clf ;
setwinsize(gcf,600,550)
x_ax = ( 1 : 400 ) * 5 ;
y_ax = ( 1 : 400 ) * 5 ;
imagescnan( x_ax, y_ax, nw_img', [ log10( min( opt.contrast_planet ) ) - 2.5, log10( max( opt.contrast_planet ) ) ] ) ;
%imagescnan( nw_img, [ ( min( opt.contrast_planet ) ) / 1000, 10*( max( opt.contrast_planet ) ) ] ) ;
xlabel( 'mas' ) ;
ylabel( 'mas' ) ;
h = colorbar ;
ylabel( h, 'Log10(Intensity)', 'FontSize', 16 )
suptitle( sprintf( 'Starshade simulation at %3.0f nm (Rayleigh scattering)', lambdaIn( i_lmbd ) / nm ) )
title( sprintf( 'Contrasts: %1.1e, %1.1e, %1.1e, %1.1e', opt.contrast_planet( 1 ), opt.contrast_planet( 2 ), opt.contrast_planet( 3 ), opt.contrast_planet( 4 ) ), 'FontSize', 16 )  
title( sprintf( 'Contrasts: %1.1e, %1.1e, %1.1e, %1.1e', opt.contrast_planet( 1 ), opt.contrast_planet( 4 ), opt.contrast_planet( 7 ), opt.contrast_planet( 10 ) ), 'FontSize', 16 )
text( 200*5, 200*5, 'A' )
text( 222*5, 200*5, 'B' )
text( 186*5, 198*5, 'C' )
text( 158*5, 203*5, 'D' )



dbstop if error
make_an_error

