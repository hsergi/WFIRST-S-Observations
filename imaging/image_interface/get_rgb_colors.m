% Getting a list of colors
function rgb = get_rgb_colors( idx )

  if idx <=0, disp( '(GET_RGB_COLORS) Color index must be positive. Defaulting to idx=1.' ) ; idx = 1 ; end

  if idx > 34, disp( '(GET_RGB_COLORS) Maximum number of colors reached!' ) ; end
idx = mod( idx, 35 ) ;
  if idx == 0, idx = 1 ; end
% From https://www.mathworks.com/matlabcentral/newsreader/view_thread/2465
rg_arry = [ 0 1 0; 1 0 1; 0 1 1; 1 0 0; .2 .6 1;
0.7 0.3 0.4; 1 .4 1; 0 0 1; 1 .2 .6; .2 1 .6;
.6 1 .2; .6 .2 1; 1 1 0; 0 .6 0; .6 0 .6;
0 .6 .6; .6 .6 0; .7 .7 .7; .6 0 0; .2 .2 .7;
.5 .5 .5; .7 .2 .2; .2 .7 .2; 0 0 .6; .3 .3 .3;
0 .9 .4; 0 .4 .9; .9 .4 0; .4 .9 0; .9 0 .4;
.4 0 .9; .8 .5 .5; .5 .8 .5; .5 .5 .8; ] ;

rgb = rg_arry( idx, : ) ;

