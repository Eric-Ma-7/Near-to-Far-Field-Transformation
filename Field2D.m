classdef Field2D
    properties(Access=public)
        x
        y
        nx
        ny
        kx = []
        ky = []
        fx = []
        fy = []
        freq0
        lambda0
        k0
    end
    
    properties(Access=private)
        c0 = 299792458
        nfgx
        nfgy
        ds
    end
    
    methods(Access=public)
        function obj = Field2D(x_in, y_in, fx_in, fy_in, freq, unit)
            switch unit
                case 'm'
                    g = 1;
                case 'cm'
                    g = 1e-2;
                case 'mm'
                    g = 1e-3;
                case 'um'
                    g = 1e-6;
                case 'nm'
                    g = 1e-9;
                case 'a'
                    g = 1e-10;
            end
            [xlist, ~, indx] = unique(x_in*g,'sorted');
            [ylist, ~, indy] = unique(y_in*g,'sorted');
            [xmesh, ymesh] = meshgrid(xlist,ylist);
            fxmesh = zeros(size(xmesh));
            fymesh = zeros(size(xmesh));
            for i = 1:length(x_in)
                fxmesh(indy(i),indx(i)) = fx_in(i);
                fymesh(indy(i),indx(i)) = fy_in(i);
            end
            obj.x = x_in;
            obj.y = y_in;
            obj.nx = fx_in;
            obj.ny = fy_in;
            obj.freq0 = freq;
            obj.lambda0 = obj.c0/freq;
            obj.k0 = 2*pi/obj.c0*freq;
            obj.nfgx = length(xlist);
            obj.nfgy = length(ylist);
            obj.ds = (xlist(2) - xlist(1))*(ylist(2) - ylist(1));
        end
        
        function obj = getFarFieldCPU(obj, fargrid, varargin)
            narginchk(2,3);
            Field2D.displog('Using CPU to calculate far field')
            poolobj = gcp('nocreate');
            if isempty(poolobj)
                Field2D.displog('No parallel pool founded')
                Field2D.displog('Starting new parallel pool')
                if nargin == 3
                    poolobj = parpool('local',varargin{1});
                else
                    poolobj = parpool('local');
                end
                Field2D.displog('Parallel pool started');
            else
                Field2D.displog('A parallel pool is running');
            end
            Field2D.displog(['Parallel pool size : ' num2str(poolobj.NumWorkers)]);
            Field2D.displog('Calculation of far field Started');
            dS = obj.ds;
            if iscell(fargrid)
                kkx = fargrid{1};
                kky = fargrid{2};
            else
                kkx = linspace(-1,1,fargrid);
                kky = linspace(-1,1,fargrid);
            end
            [Kx,Ky] = meshgrid(kkx,kky);
            Kz = real(sqrt(1-Kx.^2-Ky.^2));
            r0 = 1e6*obj.lambda0;
            % far field coordinate
            xf = r0 * Kx;
            yf = r0 * Ky;
            zf = r0 * Kz;
            inUnitCircle = Kx.^2+Ky.^2<1;
            % far field components
            ffx = zeros(size(Kx));
            ffy = zeros(size(Kx));
            k = obj.k0;
            % near field coordinate
            xn = obj.x;
            yn = obj.y;
            % near field components
            nfx = obj.nx;
            nfy = obj.ny;
            % field integral
            parfor i = 1:(size(Kx,1)*size(Kx,2))
                if inUnitCircle(i)
                    r = sqrt((xf(i)-xn).^2+(yf(i)-yn).^2+zf(i).^2);
                    ir = 1/(2*pi)*exp(1i*k*r).*zf(i)./(r.^2).*(1./r-1i*k);
                    fx_int = ir.*nfx.*dS;
                    fy_int = ir.*nfy.*dS;
                    ffx(i) = sum(fx_int(:));
                    ffy(i) = sum(fy_int(:));
                else
                    ffx(i) = 0;
                    ffy(i) = 0;
                end
            end
            Field2D.displog('Calculation of far field completed');
            obj.fx = ffx;
            obj.fy = ffy;
            obj.kx = Kx;
            obj.ky = Ky;
        end
        
        function obj = getFarFieldGPU(obj,fargrid)
            Field2D.displog('Using GPU to calculate far field');
            narginchk(2,3);
            assert(gpuDeviceCount~=0,'No GPU found');
            gpu = gpuDevice;
            Field2D.displog(['GPU name : ' gpu.Name]);
            if iscell(fargrid)
                kkx = fargrid{1};
                kky = fargrid{2};
            else
                kkx = linspace(-1,1,fargrid);
                kky = linspace(-1,1,fargrid);
            end
            Field2D.displog('Calculation of far field Started');
            dS = obj.ds;
            [Kx,Ky] = meshgrid(kkx,kky);
            Kz = real(sqrt(1-Kx.^2-Ky.^2));
            r0 = 1e6*obj.lambda0;
            % far field coordinate
            xf = r0 * Kx;
            yf = r0 * Ky;
            zf = r0 * Kz;
            inUnitCircle = Kx.^2+Ky.^2<1;
            % far field components
            ffx = zeros(size(Kx),'gpuArray');
            ffy = zeros(size(Kx),'gpuArray');
            k = obj.k0;
            % near field coordinate
            xn = gpuArray(obj.x);
            yn = gpuArray(obj.y);
            % near field components
            nfx = gpuArray(obj.nx);
            nfy = gpuArray(obj.ny);
            % field integral
            for i = 1:(size(Kx,1)*size(Kx,2))
                if inUnitCircle(i)
                    r = sqrt((xf(i)-xn).^2+(yf(i)-yn).^2+zf(i).^2);
                    ir = 1/(2*pi)*exp(1i*k*r).*zf(i)./(r.^2).*(1./r-1i*k);
                    fx_int = ir.*nfx.*dS;
                    fy_int = ir.*nfy.*dS;
                    ffx(i) = sum(fx_int(:));
                    ffy(i) = sum(fy_int(:));
                else
                    ffx(i) = 0;
                    ffy(i) = 0;
                end
            end
            Field2D.displog('Calculation of far field completed');
            Field2D.displog('Transfering array to local workspace');
            ffx = gather(ffx);
            ffy = gather(ffy);
            Field2D.displog('Array transferred');
            obj.fx = ffx;
            obj.fy = ffy;
            obj.kx = Kx;
            obj.ky = Ky;
        end
        
        function varargout = plotNearField(obj,varargin)
            narginchk(1,2);
            if nargin == 3
                clrmap = varargin{1};
            else
                clrmap = 'parula';
            end
            P = abs(obj.nx).^2 + abs(obj.ny).^2;
            s = pcolor(obj.x,obj.y,P);
            set(s,'FaceColor','interp','EdgeAlpha',0);
            axis equal off
            colormap(clrmap);
            if nargout == 1
                varargout{1} = s;
            end
        end
        
        function varargout = plotFarField(obj,quiverSize,density,range)
            P = abs(obj.fx).^2 + abs(obj.fy).^2;
            Pmax = max(P(:));
            th = linspace(0,2*pi,1000);
            ax = abs(obj.fx);
            ay = abs(obj.fy);
            A = max([max(ax(:)) max(ay(:))]);
            ax = ax/A;
            ay = ay/A;
            deltax = angle(obj.fx);
            deltay = angle(obj.fy);
            kkx = unique(obj.kx(:));
            dkx = kkx(2)-kkx(1);
            
            hold on
            s = pcolor(obj.kx,obj.ky,P);
            set(s,'FaceColor','interp','EdgeAlpha',0);
            axis equal off
            colormap hot
            for i = 1:density:size(obj.kx,1)
                for j = 1:density:size(obj.kx,2)
                    if P(i,j)>Pmax*range
                        ex = obj.kx(i,j)+ax(i,j)*cos(th+deltax(i,j))*dkx*quiverSize;
                        ey = obj.ky(i,j)+ay(i,j)*cos(th+deltay(i,j))*dkx*quiverSize;
                        plot(ex,ey,'g')
                    end
                end
            end
            hold off
            if nargout == 1
                varargout{1} = s;
            end
        end
    end
    
    methods(Static, Access=private)
        function displog(str)
            clk = datestr(now,'yyyy/mm/dd HH:MM:SS');
            disp(['[' clk '] ' str]);
        end
    end
    
end
