function ZTD = IGGtropSERA_bl(siteLon, siteLat, sitehgt, doy)
% (c) State Key Laboratory of Geodesy and Earth's Dynamics, Institute of Geodesy and Geophysics, 
% Chinese Academy of Sciences, 2017

% This program calculates ZTD from IGGtropS model, which is established by Li Wei. 
% The coefficients of IGGtropS are derived using ERA-Interim reanalysis pressure level data during 2006.1-2009.12

% References:
% Li W, Yuan YB, Ou JK, He YJ (2018). IGGtrop_SH and IGGtrop_rH: Two improved empirical tropospheric delay models based on vertical 
% reduction functions. IEEE Transactions on Geoscience and Remote Sensing. PP. 1-13. 10.1109/TGRS.2018.2812850.
% Li W, Yuan YB, Ou JK, Li H, Li ZS (2012). A new global zenith tropospheric delay model IGGtrop for GNSS applications. Chin Sci Bull, 57(17): 2132¨C2139

%==========================input============================
% sitelon:longitude [degree]; range: 0 to 360 degree
% sitelat:latitude [degree] ; range: -90 to 90 degree
% siteh:orthometric height [m]; 
% doy: day of year
%===========================================================
%==========================output===========================
% ZTD: zenith tropospheric delay [m]
%===========================================================

% example
% siteLon = 114.0;128.9;
% siteLat = 30.0;-89;
% sitehgt = 30;
% doy = 1:1:365;

    nH = fix(26);   nPara = fix(5);
    dlon = 2.5;   dlat = 2.5;
    nLon = fix(360/dlon);   nLat = fix(180/dlat)+1;
    lon = dlon*(0:1:nLon-1);   lat = 90-dlat*(0:1:nLat-1);
    dH = 1000.0;  normH = dH*(0:1:nH-1);
    
    fid = fopen('IGGtropSERA.bin');
    coeff = fread(fid, 'float');
    fclose(fid);
    coeff = reshape(coeff, [nLon nLat nH nPara]);
    
    ix = (siteLon-lon(1))/(lon(2)-lon(1)) + 1;
    iy = (siteLat-lat(1))/(lat(2)-lat(1)) + 1;
    
    Y_grid = zeros( size(doy,2), 4);
    
    for igrid=1:4 
        a=mod( fix(ix+bitand(igrid-1, 1))-1, nLon) + 1;
        b=fix(iy)+ fix((igrid-1)/2);
        if b>nLat
            b=nLat;
        end
        
        for idoy = 1:size(doy,2)
            t=2*pi/365.25*doy(idoy);
            Y_h = coeff(a,b,:,1) + coeff(a,b,:,2)*cos(t) + coeff(a,b,:,3)*sin(t) + ...
                                 + coeff(a,b,:,4)*cos(2*t) + coeff(a,b,:,5)*sin(2*t);
            Y_h = reshape(Y_h,[1,nH]);                 
            Y_grid(idoy, igrid) = exp((interp1(normH, log(Y_h), sitehgt)));
        end    
    end
    
    p = ix-fix(ix);
    q = iy-fix(iy);
    
    ZTD = (1-p)*(1-q)*Y_grid(:,1) + (1-q)*(p)*Y_grid(:,2) +  (1-p)*(q)*Y_grid(:,3)+ (p)*(q)*Y_grid(:,4);
end






