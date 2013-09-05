classdef angleOfArrivalStats
   
    methods (Static)
        
        function  out = variance(atm,d,geometry)
            %% VARIANCE variance of the angle of arrival
            %
            % out = variance(atm,d,geometry) computes the variance of the
            % angle of arrival for an atmosphere atm and a subaperture size
            % d. The subaperture is either a square or a circle
            
            switch lower(geometry)
                case 'square_cov'
                    funInt = @(x,y) angleOfArrivalStats.pointCovariance(hypot(x,y),atan2(y,x),atm).*...
                        (1-abs(x/d)).*(1-abs(y/d));
                    out = integral2( funInt, -d, d, -d, d);
                    out = out./d.^2;
                case 'square_psd'
                    funInt = @(fx,fy) fx.^2.*phaseStats.spectrum(hypot(fx,fy),atm).*tools.sinc(d.*fx).^2.*tools.sinc(d.*fy).^2;
                    lim = inf;
                    out = integral2(funInt,-lim,lim,-lim,lim);
                    out = atm.wavelength.^2.*out;
                case 'circle'
                    funInt = @(f) pi.*f.^3.*phaseStats.spectrum(f,atm).*(2.*tools.sombrero(1,pi.*d.*f)).^2;
                    out = integral(funInt,0,inf);
                    out = atm.wavelength.^2.*out;
                case 'series'
                    c = gamma(11/6)^2*(24*gamma(6/5)/5)^(5/6)/(2*pi^(11/3));
                    red = pi*d/atm.L0;
                    out = 4*pi*c.*atm.r0^(-5/3)*(pi.*d)^(-4)*atm.L0^(11/3).*...
                        tools.oneParameterExample3(2,1,2,11/6,red,10);
                    out = atm.wavelength.^2.*out;
            end
        end
        
        function out = pointCovariance(r,o,atm)
            r(r==0) = eps;
            cst = gamma(11/6)*(24*gamma(6/5)/5)^(5/6)/(2^(5/6)*pi^(8/3));
            red = 2*pi.*r/atm.L0;
            out = (2*cos(o).^2/3 + sin(o).^2).*red.^(-1/6).*besselk(1/6,red) - ...
                cos(o).^2.*red.^(5/6).*besselk(5/6,red);
            out = cst.*atm.wavelength.^2.*atm.r0.^(-5/3).*atm.L0.^(-1/3).*out;
        end
        
        function out = covariance(r,o,atm,d,dir,method)
            %% COVARIANCE covariance of the angle of arrival
            %
            % out = covariance(r,o,atm,d,dir) computes the covariance of
            % the angle of arrival at a baseline defined by the polar
            % coordinates (r,o), for an atmosphere atm and a subaperture
            % size d. dir is 0 for the x-axis angle of arrival and pi/2 for
            % the y-axis.
            
            if nargin<6 || isempty(method)
                
                funInt = @(f,r_,o_) pi.*f.^3.*phaseStats.spectrum(f,atm).*(2.*tools.sombrero(1,pi.*d.*f)).^2.*...
                    ( besselj(0,2*pi*f.*r_) - cos(2*(o_-dir)).*besselj(2,2*pi*f.*r_) );
                out = arrayfun( @(r__,o__) integral(@(x) funInt(x,r__,o__) ,0,inf) , r, o);
                out = atm.wavelength.^2.*out;
               
            else
                
                index = r~=0;
                out = ones(size(r)).*angleOfArrivalStats.variance(atm,d,'series');
                r = r(index);
                c = gamma(11/6)^2*(24*gamma(6/5)/5)^(5/6)/(2*pi^(11/3));
                red = 2*pi*r;
                c = pi*c.*atm.r0^(-5/3)*(red).^(-4)*atm.L0^(11/3);
                red = red/atm.L0;
                a = tools.oneParameterExample2(4,0,2,11/6,red,10);
                b = tools.oneParameterExample2(4,2,2,11/6,red,10);
                out(index) = atm.wavelength.^2.*c.*(a - cos(2*(o(index)-dir)).*b);
                
            end
        end
        
        function out = xyCovariance(r,o,atm,d,~,method)
            %% XYCOVARIANCE covariance between the x and y angle of arrival
            %
            % out = xyCovariance(r,o,atm,d) computes the covariance between
            % the x and y components of the angle of arrival at a baseline
            % defined by the polar coordinates (r,o) for an atmosphere atm
            % and a subaperture size d.
            
            if nargin<6 || isempty(method)
                
                funInt = @(f,r_,o_) pi.*f.^3.*phaseStats.spectrum(f,atm).*(2.*tools.sombrero(1,pi.*d.*f)).^2.*...
                    sin(2*(o_)).*besselj(2,2*pi*f.*r_);
                out = arrayfun( @(r__,o__) integral(@(x) funInt(x,r__,o__) ,0,inf) , r, o);
                
            else
                
                index = r~=0;
                out = zeros(size(r));
                r = r(index);
                c = gamma(11/6)^2*(24*gamma(6/5)/5)^(5/6)/(2*pi^(11/3));
                red = 2*pi*r;
                c = pi*c.*atm.r0^(-5/3)*(red).^(-4)*atm.L0^(11/3);
                red = red/atm.L0;
                a = tools.oneParameterExample2(4,2,2,11/6,red,10);
                out(index) = c.*sin(2*o(index)).*a;
                
            end
            
            out = -atm.wavelength.^2.*out;
        end
        
        function out = phaseCovariance(r,o,atm,d,dir,method)
            %% PHASECOVARIANCE covariance between the wavefront and the angle of arrival
            %
            % out = phaseCovariance(r,o,atm,d,dir) computes the covariance
            % between the wavefront and the angle of arrival at a baseline
            % defined by the polar coordinates (r,o), for an atmosphere atm
            % and a subaperture size d. dir is 0 for the x-axis angle of
            % arrival and pi/2 for the y-axis.
            
            if nargin<6 || isempty(method)
                
                funInt = @(f,r_,o_) 2*pi.*f.^2.*phaseStats.spectrum(f,atm).*(2.*tools.sombrero(1,pi.*d.*f)).^2.*cos(o_-dir).*besselj(1,2*pi*f.*r_);
                out = arrayfun( @(r__,o__) integral(@(x) funInt(x,r__,o__) ,0,inf) , r, o);
                
            else
                
                index = r~=0;
                out = zeros(size(r));
                r = r(index);
                c = gamma(11/6)^2*(24*gamma(6/5)/5)^(5/6)/(2*pi^(11/3));
                red = 2*pi*r;
                c = 2*pi*c.*atm.r0^(-5/3)*(red).^(-3)*atm.L0^(11/3);
                red = red/atm.L0;
                a = tools.oneParameterExample2(3,1,2,11/6,red,10);
                out(index) = c.*cos(o(index)-dir).*a;
                
            end
            
            out = atm.wavelength.*out;
        end
        
        function out = tipTiltCovarianceZ(r,o,atm,d1,g1,d2,g2,dir)
            d1 = g1.*d1;
            if nargin>7
                funInt = @(f,r_,o_) f.^2.*phaseStats.spectrum(f,atm).*...
                    tools.sombrero(1,pi.*d1.*f).*tools.sombrero(2,pi.*d2.*f).*...
                    ( besselj(0,2*pi*f.*r_)-cos(2*(o_-dir)).*besselj(2,2*pi*f.*r_));
            else
                funInt = @(f,r_,o_) f.^2.*phaseStats.spectrum(f,atm).*...
                    tools.sombrero(1,pi.*d1.*f).*tools.sombrero(2,pi.*d2.*f).*...
                    sin(2*o_).*besselj(2,2*pi*f.*r_);
            end
            out = g2.*arrayfun( @(r__,o__) integral(@(x) funInt(x,r__,o__) ,0,inf) , r, o);
            out = 16*atm.wavelength^2*out/d2;
        end
         
        function out = tipTiltPhaseCovarianceZ(r,o,atm,d1,g1,d2,g2,dir)
            d1 = g1.*d1;
            funInt = @(f,r_,o_) f.*phaseStats.spectrum(f,atm).*...
                tools.sombrero(2,pi.*d2.*f).*cos(o_-dir).*besselj(1,2*pi*f.*r_);
            out = g2.*arrayfun( @(r__,o__) integral(@(x) funInt(x,r__,o__) ,0,inf) , r, o);
            out = 16*atm.wavelength*out/d2;
        end
        
        function out = tipTiltCovarianceMatrix(x1,y1,d1,xa,ya,da,src1,x2,y2,d2,xb,yb,db,src2,atm)
            
            layers = atm.layer;
            x1 = x1(:);
            y1 = y1(:);
            x2 = x2(:);
            y2 = y2(:);
            
            out = cellfun( @(x) 0 , cell(1,7), 'UniformOutput', false);
            
            for kLayer=1:atm.nLayer
                
                theta1 = src1.directionVector.*layers(kLayer).altitude;
                scale1  = 1 - layers(kLayer).altitude./src1.height;
                z1      = complex( scale1*x1 + theta1(1) , scale1*y1 + theta1(2) );
                z1a     = complex( scale1*xa + theta1(1) , scale1*ya + theta1(2) );
                
                theta2 = src2.directionVector.*layers(kLayer).altitude;
                scale2  = 1 - layers(kLayer).altitude./src2.height;
                z2      = complex( scale2*x2 + theta2(1) , scale2*y2 + theta2(2) );
                z2b     = complex( scale2*xb + theta2(1) , scale2*yb + theta2(2) );
                                
                % first row
                z = bsxfun( @minus, z2b , z1.'); 
                r = abs(z);
                o = angle(z);
                
                out{1} = out{1} + angleOfArrivalStats.tipTiltCovarianceZ(r,o,atm,d1,scale1,db,scale2,0);                
                out{2} = out{2} + angleOfArrivalStats.tipTiltCovarianceZ(r,o,atm,d1,scale1,db,scale2,pi/2);
                out{3} = out{3} + angleOfArrivalStats.tipTiltCovarianceZ(r,o,atm,d1,scale1,db,scale2);
                                
                % first row
                z = bsxfun( @minus, z2 , z1a.'); 
                r = abs(z);
                o = angle(z);
                
                out{4} = out{4} + angleOfArrivalStats.tipTiltCovarianceZ(r,o,atm,d2,scale1,da,scale2,0);                
                out{5} = out{5} + angleOfArrivalStats.tipTiltCovarianceZ(r,o,atm,d2,scale1,da,scale2,pi/2);
                out{6} = out{6} + angleOfArrivalStats.tipTiltCovarianceZ(r,o,atm,d2,scale1,da,scale2);
                
            end
            out{7} = zernikeStats.angularCovariance(zernike(2:3,d2),atm,[src1,src2]);
            
        end
        
        function out = tipTiltPhaseCovarianceMatrix(xa,ya,da,src1,x2,y2,d2,src2,atm)
            
            layers = atm.layer;
            x2 = x2(:);
            y2 = y2(:);
            
            out = cellfun( @(x) 0 , cell(1,2), 'UniformOutput', false);
            
            for kLayer=1:atm.nLayer
                
                theta1 = src1.directionVector.*layers(kLayer).altitude;
                scale1  = 1 - layers(kLayer).altitude./src1.height;
                z1a     = complex( scale1*xa + theta1(1) , scale1*ya + theta1(2) );
                
                theta2 = src2.directionVector.*layers(kLayer).altitude;
                scale2  = 1 - layers(kLayer).altitude./src2.height;
                z2      = complex( scale2*x2 + theta2(1) , scale2*y2 + theta2(2) );
                                
                % first row
                z = bsxfun( @minus, z2 , z1a.'); 
                r = abs(z);
                o = angle(z);
                
                out{1} = out{1} + angleOfArrivalStats.tipTiltPhaseCovarianceZ(r,o,atm,d2,scale1,da,scale2,0);                
                out{2} = out{2} + angleOfArrivalStats.tipTiltPhaseCovarianceZ(r,o,atm,d2,scale1,da,scale2,pi/2);  
                
            end
            
        end
       
        function out = angularCovarianceMatrix(x1,y1,src1,x2,y2,src2,atm,d,dir,covFun,method)
            
            if nargin<11
                method = [];
            end
            layers = atm.layer;
            out = 0;
            x1 = x1(:);
            y1 = y1(:);
            x2 = x2(:);
            y2 = y2(:);
            for kLayer=1:atm.nLayer
                
                theta1 = src1.directionVector.*layers(kLayer).altitude;
                scale1  = 1 - layers(kLayer).altitude./src1.height;
                z1      = complex( scale1*x1 + theta1(1) , scale1*y1 + theta1(2) );
                
                theta2 = src2.directionVector.*layers(kLayer).altitude;
                scale2  = 1 - layers(kLayer).altitude./src2.height;
                z2      = complex( scale2*x2 + theta2(1) , scale2*y2 + theta2(2) );
                
                z = bsxfun( @minus, z2 , z1.');
                r = abs(z);
                o = angle(z);
                
                out = out + ...
                    layers(kLayer).fractionnalR0*covFun(r,o,atm,d*scale1,dir,method);
            end
        end
        
        function TBT = angularCovarianceMatrixToeplitz(x1,y1,src1,x2,y2,src2,atm,d,dir,covFun,method)
            
            if nargin<11
                method = [];
            end
            n1 = length(x1);
            n2 = length(x2);
            layers = atm.layer;
            cov = 0;
            x1 = x1(:);
            y1 = y1(:);
            x2 = x2(:);
            y2 = y2(:);
            for kLayer=1:atm.nLayer
                
                theta1 = src1.directionVector.*layers(kLayer).altitude;
                scale1  = 1 - layers(kLayer).altitude./src1.height;
                z1      = complex( scale1*x1 + theta1(1) , scale1*y1 + theta1(2) );
                
                theta2 = src2.directionVector.*layers(kLayer).altitude;
                scale2  = 1 - layers(kLayer).altitude./src2.height;
                z2      = complex( scale2*x2 + theta2(1) , scale2*y2 + theta2(2) );
                                
                % first row
                z = bsxfun( @minus, z2(1) , z1.'); 
                z = reshape(z,n1,n1);
                
                % combined with first column of first block row
                zr = [flipud(z);bsxfun( @minus, z2(2:n2) , z1(1:n1:end).')];
                
                % first column
                z = bsxfun( @minus, z2(n2+1:end) , z1(1) );
                z = reshape(z,n2,n2-1);
                
                % combined with first row of first block column
                zc = [flipud(bsxfun( @minus, z2(n2+1:n2:end) , z1(2:n1).').');z];
                
                z = [fliplr(zr) zc];
                r = abs(z);
                o = angle(z);
                
                cov = cov + ...
                    layers(kLayer).fractionnalR0*covFun(r,o,atm,d*scale1,dir,method);
            end
            TBT = toeplitzBlockToeplitz([n2,n1],[n2,n1],cov);
        end
                
        function MM = fullToeplitx(R)
            
            n = size(R,1);
            nb = 0.5*(n+1);
            
            cidx = (0:nb-1)';
            ridx = nb:-1:1;
            t = cidx(:,ones(nb,1)) + ridx(ones(nb,1),:);  % Toeplitz subscripts

            M = cell(1,size(R,2));
            for k=1:size(R,2)
                x = R(:,k);
                M{k} = x(t);
            end
            
            MM = cell2mat(M(t));
        end
        
        function out = angularCovarianceMatrixBlock(x,y,src1,src2,atm,d,method)
            if nargin<7
                method = [];
            end
            fprintf('@(angleOfArrivalStats)> angle of arrival xx covariance ...\n')
            out{1} = angleOfArrivalStats.angularCovarianceMatrixToeplitz(x,y,src1,x,y,src2,atm,d,0,@angleOfArrivalStats.covariance,method);
            fprintf('@(angleOfArrivalStats)> angle of arrival yy covariance ...\n')
            out{2} = angleOfArrivalStats.angularCovarianceMatrixToeplitz(x,y,src1,x,y,src2,atm,d,pi/2,@angleOfArrivalStats.covariance,method);
            fprintf('@(angleOfArrivalStats)> angle of arrival xy covariance ...\n')
            out{3} = angleOfArrivalStats.angularCovarianceMatrixToeplitz(x,y,src1,x,y,src2,atm,d,0,@angleOfArrivalStats.xyCovariance,method);
        end
        
        function out = fullBlock(B)
            Cxx = full( B{1} );
            Cyy = full( B{2} );
            Cxy = full( B{3} );
            out = [ Cxx , Cxy ; Cxy , Cyy ];
        end
        
        function out = angularCovarianceMatrixMetaBlock(x,y,src,atm,d,method)
            if nargin<6
                method = [];
            end
            nSrc = length(src);
            fprintf('@(angleOfArrivalStats)> AA covariance meta matrix:\n')
            out = cell(nSrc);
            out{1} = angleOfArrivalStats.angularCovarianceMatrixBlock(x,y,src(1),src(1),atm,d,method);
            for kSrc2=1:nSrc-1
                for kSrc1=kSrc2+1:nSrc
                    fprintf(' [%d,%d] ;',kSrc2,kSrc1)
                    out{kSrc2,kSrc1} = angleOfArrivalStats.angularCovarianceMatrixBlock(x,y,src(kSrc1),src(kSrc2),atm,d,method);
                end
                fprintf('\b\n')
            end
        end
        
        function S = fullMetaBlock(Sin)
            S = Sin;
            n = length(S);
            idx = reshape(1:n^2,n,n);
            upDiag = triu(idx,1)>0;
            downDiag = tril(idx,-1)>0;
            S( downDiag ) = S( upDiag );
            S(1:n+1:end) = S(1);
            S = cellfun( @(x) angleOfArrivalStats.fullBlock(x), S , 'UniformOutput', false);
            S( downDiag) = cellfun( @(x) x', S(downDiag) , 'UniformOutput', false);
            S = cell2mat( S );
        end
        
        function out = angularPhaseCovarianceMatrixMetaBlock(x1,y1,src1,x2,y2,src2,atm,d,method)
            if nargin<9
                method = [];
            end
            nSrc1 = length(src1);
            nSrc2 = length(src2);
            out = cell(nSrc2,nSrc1);
            fprintf('@(angleOfArrivalStats)> phase covariance meta matrix:\n')
            for kSrc2=1:nSrc2
                for kSrc1=1:nSrc1
                    fprintf(' [%d,%d] ;',kSrc2,kSrc1)
                    out{kSrc2,kSrc1} = {...
                        angleOfArrivalStats.angularCovarianceMatrixToeplitz(...
                        x1,y1,src1(kSrc1),x2,y2,src2(kSrc2),atm,d,0,...
                        @angleOfArrivalStats.phaseCovariance,method) , ...
                        angleOfArrivalStats.angularCovarianceMatrixToeplitz(...
                        x1,y1,src1(kSrc1),x2,y2,src2(kSrc2),atm,d,pi/2,...
                        @angleOfArrivalStats.phaseCovariance,method) };
                end
                fprintf('\b\n')
            end
        end
        
        function out = fullPhaseMetaBlock(T)
            out = cellfun( @(x) [ full( x{1} ) , full( x{2} ) ] , T , 'uniformOutput', false);
            out = cell2mat(out);
        end
        
        function displayStats
            atm = atmosphere(589e-9,15e-2,30);
            [x,y] = meshgrid( linspace(-5,5,41) );
            [o,r] = cart2pol(x,y);
            Cxx = angleOfArrivalStats.covariance(r,o,atm,1,0);
            Cyy = angleOfArrivalStats.covariance(r,o,atm,1,pi/2);
            Cxy = angleOfArrivalStats.xyCovariance(r,o,atm,1);
            Cpx = angleOfArrivalStats.phaseCovariance(r,o,atm,1,0);
            Cpy = angleOfArrivalStats.phaseCovariance(r,o,atm,1,pi/2);
            figure
            subplot(2,1,1)
            imagesc([Cxx,Cyy])
            axis equal tight
            subplot(2,1,2)
            imagesc([Cpx,Cpy])
            axis equal tight
            figure
            imagesc(Cxy)
            axis equal tight
        end

        function tipTiltStats
%%
    atm = atmosphere(589e-9,15e-2,30);
            nLenslet = 40;
            tel = telescope(10,'resolution',nLenslet*6);
            tel.shape = 'square';
            lenslet = telescope(tel.D/nLenslet,'resolution',nLenslet);
            zern = zernike(lenslet,2:3);
            wfs = shackHartmann(nLenslet,tel.resolution);
            ngs = source('wavelength',photometry.K);
            ngs = ngs.*tel*wfs;
            wfs.INIT
            tel = tel + atm;
            lenslet = lenslet + atm;
            %%
            [x,y] = meshgrid( linspace(-1,1,nLenslet)*tel.R );
            [o,r] = cart2pol(x,y);
            Cxa2 = angleOfArrivalStats.tipTiltCovarianceZ(r,o,atm,1,1,1,1,0);
            Cya3 = angleOfArrivalStats.tipTiltCovarianceZ(r,o,atm,1,1,1,1,pi/2);
            Cxa3 = angleOfArrivalStats.tipTiltCovarianceZ(r,o,atm,1,1,1,1);
            Cpa2 = angleOfArrivalStats.tipTiltPhaseCovarianceZ(r,o,atm,1,1,1,1,0);
            Cpa3 = angleOfArrivalStats.tipTiltPhaseCovarianceZ(r,o,atm,1,1,1,1,pi/2);
            figure
            subplot(2,1,1)
            imagesc([Cxa2,Cya3])
            axis equal tight
            subplot(2,1,2)
            imagesc([Cpa2,Cpa3])
            axis equal tight
            figure
            imagesc(Cxa3)
            axis equal tight
            %%
            C = 0;
            Csa = 0;
            Cpa = 0;
            nIt = 1000;
            h = waitbar(0,'Matrix estimation...!');
            for k=1:nIt
                +tel;
                ngs = ngs.*tel*wfs;
                p = reSample(ngs,nLenslet+1);
                ngs = ngs.*lenslet;
                zern = zern.\ngs;
                C = C + wfs.slopes*wfs.slopes';
                Csa = Csa + wfs.slopes*zern.c';
                Cpa = Cpa + p(:)*zern.c';
                waitbar(k/nIt)
            end
           close(h)
       end
        
        function tomography
            atm = atmosphere(photometry.V,15e-2,30,'altitude',10e3);
            atm.wavelength = photometry.K;
            ngs = source('wavelength',photometry.K);
            lgs = source('asterism',{[3,arcsec(30),0]},'height',90e3,'wavelength',photometry.K);
            nLenslet = 10;
            resolution = nLenslet*6;
            wfs = shackHartmann(nLenslet,resolution);
            tel = telescope(10,'resolution',resolution,'fieldOfViewInArcmin',1);
            ngs = ngs.*tel*wfs;
            wfs.INIT;
            u = 0:nLenslet-1;
            [x,y] = meshgrid(u);
            tel = tel + atm;
            lgs = lgs.*tel*wfs;
            d = 2*(tel.D/nLenslet)/sqrt(pi);
            S = angleOfArrivalStats.angularCovarianceMatrixMetaBlock(x,y,lgs,atm,d);
            Sf = angleOfArrivalStats.fullMetaBlock(S);
            figure,imagesc(Sf),axis equal tight
            C = angleOfArrivalStats.angularPhaseCovarianceMatrixMetaBlock(x,y,ngs,lgs,atm,d);
            Cf = angleOfArrivalStats.fullPhaseMetaBlock(C);
            figure,imagesc(Cf),axis equal tight
            M = Cf/Sf;
            slopes = wfs.slopes;
            slopes(isnan(slopes)) = 0;
            phase = M*slopes(:);
            phase = reshape(phase,nLenslet,nLenslet);
            figure,imagesc(phase),axis square
            ngs = ngs.*tel;
            figure,imagesc(ngs.phase),axis square
        end
        
        function scao
            %%
            atm = atmosphere(photometry.V,15e-2,30,'altitude',10e3);
            atm.wavelength = photometry.K;
            ngs = source('wavelength',photometry.K);
            nLenslet = 20;
            resolution = nLenslet*6;
            wfs = shackHartmann(nLenslet,resolution);
            tel = telescope(10,'resolution',resolution,'fieldOfViewInArcmin',1);
            ngs = ngs.*tel*wfs;
            wfs.INIT;
            u = linspace(-1,1,nLenslet)*0.5*tel.D*(1-1/nLenslet);
            [x,y] = meshgrid(u);
            tel = tel + atm;
            ngs = ngs.*tel*wfs;
            d = 2*(tel.D/nLenslet)/sqrt(pi);
            method = 'series';
            S = angleOfArrivalStats.angularCovarianceMatrixMetaBlock(x,y,ngs,atm,d,method);
%             S_xx = angleOfArrivalStats.angularCovarianceMatrix(...
%                 x,y,ngs,x,y,ngs,atm,d,0,@angleOfArrivalStats.covariance);
%             S_yy = angleOfArrivalStats.angularCovarianceMatrix(...
%                 x,y,ngs,x,y,ngs,atm,d,pi/2,@angleOfArrivalStats.covariance);
%             S_xy = angleOfArrivalStats.angularCovarianceMatrix(...
%                 x,y,ngs,x,y,ngs,atm,d,0,@angleOfArrivalStats.xyCovariance);
            Sf = angleOfArrivalStats.fullMetaBlock(S);
            figure,imagesc(Sf),axis equal tight
            v = linspace(-1,1,nLenslet+1).*tel.R;
%             v = linspace(-1,1,resolution).*tel.R;
            [xPhase,yPhase] = meshgrid(v);
            C = angleOfArrivalStats.angularPhaseCovarianceMatrixMetaBlock(...
                x,y,ngs,xPhase,yPhase,ngs,atm,d,method);
%             C_px = angleOfArrivalStats.angularCovarianceMatrix(...
%                 x,y,ngs,xPhase,yPhase,ngs,atm,d,0,@angleOfArrivalStats.phaseCovariance,'series');
%             C_py = angleOfArrivalStats.angularCovarianceMatrix(...
%                 x,y,ngs,xPhase,yPhase,ngs,atm,d,pi/2,@angleOfArrivalStats.phaseCovariance,'series');
            
            Cf = angleOfArrivalStats.fullPhaseMetaBlock(C);
            figure,imagesc(Cf),axis equal tight
            M = Cf/Sf;
            slopes = wfs.slopes;
            slopes(isnan(slopes)) = 0;
            phase = M*slopes(:);
            phase = reshape(phase,sqrt(length(phase)),[]);
            y = cgs(Sf,slopes(:));
            phaseCg = Cf*y;
            phaseCg = reshape(phaseCg,sqrt(length(phaseCg)),[]);
            figure,imagesc([phase,phaseCg]),axis equal tight
            ngs = ngs.*tel;
            figure,imagesc(ngs.phase),axis square
        end
        
        function out = covTimeVec(G,b)
            
            nGS = length(G);
            n = size(b,1)/nGS;
            
            b = mat2cell( b , n*ones(nGS,1) , 1 );
            b = cellfun( @(x) mat2cell(x,0.5*n*[1;1],1),b,'UniformOutput',false);
            
            out = cell(nGS,1);
            out = cellfun( @(x) cell(2,1) , out ,...
                'UniformOutput',false);
            
            parfor iGS=1:nGS
                
                C = G{1,1}{1};
                out{iGS}{1} =               C*b{iGS}{1};
                C = G{1,1}{3};
                out{iGS}{1} = out{iGS}{1} + C*b{iGS}{2};
                out{iGS}{2} =               C*b{iGS}{1};
                C = G{1,1}{2};
                out{iGS}{2} = out{iGS}{2} + C*b{iGS}{2};
                
                for k=iGS+1:nGS
                    C = G{iGS,k}{1};
                    out{iGS}{1} = out{iGS}{1} + C*b{k}{1};
                    C = G{iGS,k}{3};
                    out{iGS}{1} = out{iGS}{1} + C*b{k}{2};
                    out{iGS}{2} = out{iGS}{2} + C*b{k}{1};
                    C = G{iGS,k}{2};
                    out{iGS}{2} = out{iGS}{2} + C*b{k}{2};
                end
                
                for k=1:iGS-1
                    C = G{k,iGS}{1}.';
                    out{iGS}{1} = out{iGS}{1} + C*b{k}{1};
                    C = G{k,iGS}{3}.';
                    out{iGS}{1} = out{iGS}{1} + C*b{k}{2};
                    out{iGS}{2} = out{iGS}{2} + C*b{k}{1};
                    C = G{k,iGS}{2}.';
                    out{iGS}{2} = out{iGS}{2} + C*b{k}{2};
                end
                
            end
            out = cell2mat(...
                cellfun (@(x)cell2mat(x),out,'uniformOutput',false) ); 
        end
        
        function out = covPhaseTimeVec(G,b)
            
            nGS = length(G);
            n = size(b,1)/nGS;
            m = size(G{1}{1},1);
            
            b = mat2cell( b , n*ones(nGS,1) , 1 );
            b = cellfun( @(x) mat2cell(x,0.5*n*[1;1],1),b,'UniformOutput',false);
            
            out = zeros(m,1);
            
            for kGS=1:nGS
                
                out = out + full(G{kGS}{1})*b{kGS}{1};
                out = out + full(G{kGS}{2})*b{kGS}{2};
                
            end
                        
        end
    end
end