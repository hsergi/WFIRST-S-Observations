function plot_starshade( alpha )
% Plot the Starshade and optionally overplot some lines showing some directions

% From makeStarshdeImage
load 'in/NI2'
vecPetalArray = createErrorProfileV1_image(r, Profile, occulterDiameter, petalLength, numPetals, {});
tic
tmpxVals = [];
tmpyVals = [];
tmpzVals = [];%**
for j = 1:numPetals
    tmpxVals = [tmpxVals vecPetalArray{j}{1}(1, :)]; %#ok<AGROW>
    tmpyVals = [tmpyVals vecPetalArray{j}{1}(2, :)]; %#ok<AGROW>
    tmpzVals = [tmpzVals vecPetalArray{j}{1}(3, :)]; %#ok<AGROW> %**
end
xVals = [tmpxVals tmpxVals(1)];
yVals = [tmpyVals tmpyVals(1)];
zVals = [tmpzVals tmpzVals(1)];
% xVals, yVals, zVals give the 3D edge locus
% Under most circumstances zVals will be all zeros
t = toc ;
disp( sprintf( 'xyzVals took %3.2f seconds', t ) )
close all
figure( 1 )
setwinsize( gcf, 500, 500 )
plot( xVals, yVals )
title( 'Starshade for occulter NI2', 'FontSize', 16 )

  if exist( 'alpha', 'var' )
  n_alph = numel( alpha ) ;
  hold all
    for i_alph = 1 : n_alph
    % Using the same colors as in other programs (e.g. starshade_radial_interpolation)
    rgb_tmp = get_rgb_colors( i_alph ) ;
      if alpha( i_alph ) ~= 90
      h( i_alph ) = plot( 0 : 0.01 : 15, ( 0 : 0.01 : 15 ) * tan( alpha( i_alph ) * pi / 180. ) ) ; 
      else
      h( i_alph ) = plot( 0 * ( 0 : 0.01 : 15 ), 0 : 0.01 : 15 ) ;
      end
    alph_str{ i_alph } = sprintf( '%2.1f deg', alpha( i_alph ) ) ;
    end
  legend( h, alph_str, 'Location', 'NorthWest' ) ;
  end
whitebg()

% Storing the image:
img = getframe( gcf ) ;
  if exist( 'alpha', 'var' )
  imwrite( img.cdata, sprintf( 'fig_dev/starshade_NI2_n%i.%s', numel( alpha ), 'png' ) ) ;
  else
  imwrite( img.cdata, sprintf( 'fig_dev/starshade_NI2.%s', 'png' ) ) ;
  end
tic
