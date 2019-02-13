% Boundary diffraction wave, fast
% Eric Cady, Caltech/JPL   1/24/12
%  8/2/13: Switched to use outPoly on the clip tests in place of
%   xVals/yVals
%
% Following Miyamoto and Wolf (I and II) 1962
% Version #4: takes edges in directly, can do off-axis and out-of-plane,
% can do lateral errors, can do multiple wavelengths at once, very fast if
% you're not overlapping the geometric outline of the occulter.
%
% Note: Still assumes an underlying grid.

function E = bdwf_image(xVals, yVals, zVals, Z, lambda, dxO, nO, psi1, psi2, deltaX, deltaY)

% Set up output grid
Nx = nO;
Ny = nO;

width = nO*dxO;
xO = -width/2+dxO/2:dxO:width/2-dxO/2;
yO = xO - deltaY;
xO = xO - deltaX;

% Prep flags
if isempty(zVals) || isequal(zVals, zeros(size(zVals)))
    flagZ = false;
else
    flagZ = true;
end

if psi1 ~= 0
    flagP1 = true;
else
    flagP1 = false;
end

inOccFlag = 0;
outPoly = [xVals(:), yVals(:)]; %polyclip(xVals, yVals, min(xO)-dxO/2, max(xO)+dxO/2, min(yO)-dxO/2, max(yO)+dxO/2, 0);
% outPoly = polyclip(xVals, yVals, min(xO)-dxO/2, max(xO)+dxO/2, min(yO)-dxO/2, max(yO)+dxO/2, 0);
% if isequal(size(outPoly),[5,2]) && max(outPoly(:,1)) == max(xO)+dxO/2 && min(outPoly(:,1)) == min(xO)-dxO/2 ...
%         && max(outPoly(:,2)) == max(yO)+dxO/2 && min(outPoly(:,2)) == min(yO)-dxO/2
%     inOccFlag = 1;
% else
%     inOccFlag = 0;
% end

% Assign variables prior to loop
p2l = 2*pi./lambda;
pil = pi./lambda;
pilz = pi./(lambda*Z);
nLambda = length(lambda);
E = zeros(Nx, Ny, nLambda);
wind = zeros(Nx, Ny);
vt = outPoly;

s1 = sin(psi1);
c1 = cos(psi1);
s2 = sin(psi2);
c2 = cos(psi2);

% Use midpoints for integral, endpoints to get vector ell for the dot
% product and dl.
xm = (xVals(2:end) + xVals(1:end-1))/2;
ym = (yVals(2:end) + yVals(1:end-1))/2;

xl = xVals(2:end) - xVals(1:end-1);
yl = yVals(2:end) - yVals(1:end-1);

% Do edge integral
if flagZ
    % Out-of-plane errors present
    zm = (zVals(2:end) + zVals(1:end-1))/2;
    zl = zVals(2:end) - zVals(1:end-1);
    
    if flagP1
        % Off-axis source present, as well as out-of-plane errors
        tilt = zeros(nO, nO);
        
        for jj = 1:Nx
            for kk = 1:Ny
                
                dx = (xm - xO(jj));
                dy = (ym - yO(kk));
                dz = (zm - Z);
                
                f = -s1*dz + c1*c2*dx + c1*s2*dy;
                g =        -    s2*dx +    c2*dy;
                h = -c1*dz - s1*c2*dx - s1*s2*dy;
                
                fSquarePlusGSquare = f.*f + g.*g;
                sHatCrossPDotLdl = xl.*(f*s2 + g*c1*c2) + yl.*(-f*c2 + g*c1*s2) + zl.*(-g*s1);
                
                for qq = 1:nLambda
                    E(jj, kk, qq) = sum(exp(1i*pil(qq)*(fSquarePlusGSquare./h))./(fSquarePlusGSquare).*sHatCrossPDotLdl);
                end
                
                wind(jj, kk) = polywindFlag_image(vt, [(xO(jj) - Z*s1*c2) (yO(kk) - Z*s1*s2)], inOccFlag);
                tilt(jj, kk) = xO(jj)*c2 + yO(kk)*s2;
            end
        end
        
        for qq = 1:nLambda
            eikz = exp(1i*p2l(qq)*Z*c1)*exp(1i*p2l(qq)*s1*tilt);
            E(:,:,qq) = eikz./(2*pi).*E(:,:,qq);
            E(:,:,qq) = eikz.*(wind == 0) - E(:,:,qq);
        end
    else
        % Out-of-plane but on-axis
        for jj = 1:Nx
            for kk = 1:Ny
                h = (Z - zm);
                
                fSquarePlusGSquare = (xm - xO(jj)).^2 + (ym - yO(kk)).^2;
                sHatCrossPDotLdl = xl.*(ym - yO(kk)) - yl.*(xm - xO(jj));
                for qq = 1:nLambda
                    E(jj, kk, qq) = sum(exp(1i*pil(qq)*(fSquarePlusGSquare./h))./(fSquarePlusGSquare).*sHatCrossPDotLdl);
                end
                
                wind(jj, kk) = polywindFlag_image(vt, [xO(jj) yO(kk)], inOccFlag);
            end
        end
        
        for qq = 1:nLambda
            eikz = exp(1i*p2l(qq)*Z);
            E(:,:,qq) = eikz./(2*pi).*E(:,:,qq);
            E(:,:,qq) = eikz.*(wind == 0) - E(:,:,qq);
        end
    end
    
else
    % Entirely in-plane
    if flagP1
tic
% SRH
%Nx_0 = Nx ;
%Ny_0 = Ny ;
%Nx = 32 ;
%Ny = 32 ;
%sprintf( 'Nx, Ny=%i,%i changed to %i,%i', Nx_0, Ny_0, Nx, Ny )

        % Off-axis source present but in-plane
        tilt = zeros(nO, nO);
        
        for jj = 1:Nx
            for kk = 1:Ny
                dx = (xm - xO(jj));
                dy = (ym - yO(kk));
% 6% running time for f, g and h                
                f = s1*Z + c1*c2*dx + c1*s2*dy;
                g =      -    s2*dx +    c2*dy;
                h = c1*Z - s1*c2*dx - s1*s2*dy;
                
% 3% running time spent here:
                fSquarePlusGSquare = f.*f + g.*g;
%                sHatCrossPDotLdl = xl.*(s2*s1*Z + c1*dy) + yl.*(-c2*s1*Z - c1*dx);
% SRH
% 7% running time spent with tmp_1 and tmp_2
                tmp_1 = fSquarePlusGSquare./h ;
                tmp_2 = ( xl.*(s2*s1*Z + c1*dy) + yl.*(-c2*s1*Z - c1*dx) ) ./ fSquarePlusGSquare ;
                for qq = 1:nLambda
%                     E(jj, kk, qq) = sum(exp(1i*pil(qq)*(fSquarePlusGSquare./h))./(fSquarePlusGSquare).*sHatCrossPDotLdl);
% SRH
% 32% of the running time spent here
                    E(jj, kk, qq) = sum(exp(1i*pil(qq)*tmp_1).*tmp_2) ; 
                end
% 48% running time spent here                
                wind(jj, kk) = polywindFlag_image(vt, [(xO(jj) - Z*s1*c2) (yO(kk) - Z*s1*s2)], inOccFlag);
% Negligible time spent here
                tilt(jj, kk) = xO(jj)*c2 + yO(kk)*s2;
            end
        end
toc        
%  Negligible time spent here
        for qq = 1:nLambda
            eikz = exp(1i*p2l(qq)*Z*c1)*exp(1i*p2l(qq)*s1*tilt);
            E(:,:,qq) = eikz./(2*pi).*E(:,:,qq);
            E(:,:,qq) = eikz.*(wind == 0) - E(:,:,qq);
        end
    else
        % On-axis source & flat occulter
% SRH:
tic
%Nx_0 = Nx ;
%Ny_0 = Ny ;
%Nx = 4 ; 
%Ny = 4 ;
%sprintf( 'Nx, Ny=%i,%i changed to %i,%i', Nx_0, Ny_0, Nx, Ny ) 
        for jj = 1:Nx
            for kk = 1:Ny
                fSquarePlusGSquare = (xm - xO(jj)).^2 + (ym - yO(kk)).^2;
                sHatCrossPDotLdl = xl.*(ym - yO(kk)) - yl.*(xm - xO(jj));
% SRH
                tmp = sHatCrossPDotLdl ./ fSquarePlusGSquare ;                
                for qq = 1:nLambda
% SRH
                    E(jj,kk,qq) = sum(exp(1i*fSquarePlusGSquare*pilz(qq)).*tmp);
%                    E(jj,kk,qq) = sum(exp(1i*pilz(qq)*(fSquarePlusGSquare))./(fSquarePlusGSquare).*sHatCrossPDotLdl);
%                    D(jj,kk,qq) = E2(jj,kk,qq) - E(jj,kk,qq) ;
                end
                
                wind(jj, kk) = polywindFlag_image(vt, [xO(jj) yO(kk)], inOccFlag);
            end
        end
toc
dbstop if error
%make_an_error
        
        for qq = 1:nLambda
            eikz = exp(1i*p2l(qq)*Z);
            E(:,:,qq) = eikz./(2*pi).*E(:,:,qq);
            E(:,:,qq) = eikz.*(wind == 0) - E(:,:,qq);
        end
    end
end
%make_an_error
