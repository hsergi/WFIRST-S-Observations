% Function to translate the UTot field from makeStarshadeImage to the image plane

function efDefectImg = UTot2ImagePlane(  lambdaIn, opt, pupil, telescopeDiameter, UTotL, mssg )

tic
units_image()
Nx_img = opt.Nx_img ;
masPerPixel = opt.diam_img_mas / Nx_img ;
  if ( opt.developer )
  [ ~, ~, impeak ] = mft_image( 1, Nx_img, pupil ) ; % Use 1 L/D square
  else
  [ ~, ~, impeak ] = mft( 1, Nx_img, pupil ) ; % Use 1 L/D square
  end
peak = max(abs(impeak(:)));

efDefectImg = zeros(Nx_img, Nx_img, length(lambdaIn));
for ll = 1:length(lambdaIn)
    UTot = UTotL(:,:,ll);
    lambda = lambdaIn(ll);
    disp(['Wavelength: ' num2str(lambda*1e9) 'nm'])

    D = telescopeDiameter;
    imagePlaneDiameterInMAS = Nx_img*masPerPixel;
    imagePlaneDiameterInLambdaOverD = imagePlaneDiameterInMAS*mas*D/lambda;

      if ( opt.developer )
      [Xout, Yout, imagePlane] = mft_image(imagePlaneDiameterInLambdaOverD, Nx_img, UTot.*pupil);
      else
      [Xout, Yout, imagePlane] = mft(imagePlaneDiameterInLambdaOverD, Nx_img, UTot.*pupil);
      end
    imagePlane = imagePlane/peak; % Normalize to unocculted on-axis peak = 1

    efDefectImg(:,:,ll) = imagePlane;
end
t = toc ;
  if exist( 'mssg', 'var' )
    if ( mssg )
    disp( sprintf( 'efDefectImg took %3.2f seconds', t ) )
    end
  end
