; svn $Id$
;+
; NAME:
;   MGH_ROMS_ORTHOGONALITY
;
; PURPOSE:
;   Calculate grid orthogonality
;
; POSITIONAL PARAMETERS:
;   x (input, 2-D numeric array)
;     Longitude or Cartesian x position of rho points
;
;   y (input, 2-D numeric array)
;     Latitude or Cartesian y position of rho points
;
; KEYWORD PARAMETERS:
;   LONLAT (input, switch)
;     Specifies whether the input arrays are longitude/latitude or
;     Cartesian coordinates. Default is 1.
;
; RETURN VALUE:
;   The function returns a floating/double array with the same size and
;   shape as the input, containing the stiffness of the grid.
;
; PROCEDURE:
;   Calculate values on a staggered grid, then project back onto original grid
;   by taking maximum neighbouring value,
;
;###########################################################################
;
; This software is provided subject to the following conditions:
;
; 1.  NIWA makes no representations or warranties regarding the
;     accuracy of the software, the use to which the software may
;     be put or the results to be obtained from the use of the
;     software.  Accordingly NIWA accepts no liability for any loss
;     or damage (whether direct of indirect) incurred by any person
;     through the use of or reliance on the software.
;
; 2.  NIWA is to be acknowledged as the original author of the
;     software where the software is used or presented in any form.
;
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2010-10:
;     Written, adpating code from Matlab script orthogonality.m,
;     provided by Julai Moriarty.
;-
function mgh_orthogonality, x, y, LONLAT=lonlat

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   message, 'Not implemented yet!'

;% The program loads grid nodes,
;% computes orthogonality and grid resolution
;% and plots orthogonality OR grid resolution
;% JMM 10/2010
;clear all

;%Gridgen file
;goutput_file='/home/moriarty/pacific/Waipaoa/Modelgrids/MakeGridgen/gridgen_output';

;% Load gridgen-produced grid node locations:
;%For grid produced with a resolution file - check what x_ and y_length are
;%in gridgen_resolution.m
;%Change ny and nx as needed
;nodes=importdata(goutput_file,' ',0);
;ny=118;
;nx=292;
;x_nodes=reshape(nodes(:,1),nx,ny);
;y_nodes=reshape(nodes(:,2),nx,ny);
;clear nodes nx ny

;%For uniform grid: nx and ny are specified in gridgen_param
;% use the following lines instead:
;%output=importdata(goutput_file,' ',1);
;%nodes=output.data;
;%[nx ny]=strread(char(output.textdata),'%*s%n%*s%n');
;%x_nodes=reshape(nodes(:,1),nx,ny);
;%y_nodes=reshape(nodes(:,2),nx,ny);

;% Set map projection, Put gridgen results in lat/lon
;projection='Mercator';
;m_proj(projection)
;Radius=6371229.0;
;[lon_nodes lat_nodes]=m_xy2ll(x_nodes/Radius,y_nodes/Radius);


;% Bathymetry and Coastline files
;bathy_file='/home/moriarty/pacific/Waipaoa/Modelgrids/BATHY/depth_gridded.mat';
;coast_file='/home/moriarty/pacific/Waipaoa/Modelgrids/BATHY/coast/full_coast_nz.mat';

;load(bathy_file)
;lat_bathy=LAT;
;lon_bathy=LON;
;bathy=-1*H;
;clear LAT LON H README
;clear h h_all h_pov_shelf pov_shelf_bathy
;clear lat lat_all lat_pov_shelf lats
;clear lon lon_all lon_pov_shelf lons

;load(coast_file)
;lat_coast=coast(:,1);
;lon_coast=coast(:,2);
;clear coast

;% Compute orthogonality and resolution
;orth=zeros(size(y_nodes,1), size(y_nodes,2));
;resolution=zeros(size(y_nodes,1), size(y_nodes,2));
;for ii=2:size(y_nodes,1)-1;%107,
;    for jj=2:size(y_nodes,2)-1
;       y_up=y_nodes(ii-1,jj);
;       y_do=y_nodes(ii+1,jj);
;       y_ri=y_nodes(ii,jj+1);
;       y_le=y_nodes(ii,jj-1);
;       y_xo=y_nodes(ii,jj);

;       x_up=x_nodes(ii-1,jj);
;       x_do=x_nodes(ii+1,jj);
;       x_ri=x_nodes(ii,jj+1);
;       x_le=x_nodes(ii,jj-1);
;       x_xo=x_nodes(ii,jj);

;       length_up=sqrt(((y_xo-y_up).^2)+((x_xo-x_up).^2));
;       length_do=sqrt(((y_xo-y_do).^2)+((x_xo-x_do).^2));
;       length_ri=sqrt(((y_xo-y_ri).^2)+((x_xo-x_ri).^2));
;       length_le=sqrt(((y_xo-y_le).^2)+((x_xo-x_le).^2));
;       length_ur=sqrt(((y_ri-y_up).^2)+((x_ri-x_up).^2));
;       length_dr=sqrt(((y_ri-y_do).^2)+((x_ri-x_do).^2));
;       length_dl=sqrt(((y_do-y_le).^2)+((x_do-x_le).^2));
;       length_ul=sqrt(((y_up-y_le).^2)+((x_up-x_le).^2));
;%        length_up=spheriq_dist(lon_xo,lat_xo,lon_up,lat_up);
;%        length_do=spheriq_dist(lon_xo,lat_xo,lon_do,lat_do);
;%        length_ri=spheriq_dist(lon_xo,lat_xo,lon_ri,lat_ri);
;%        length_le=spheriq_dist(lon_xo,lat_xo,lon_le,lat_le);
;%        length_ur=spheriq_dist(lon_ri,lat_ri,lon_up,lat_up);
;%        length_dr=spheriq_dist(lon_ri,lat_ri,lon_do,lat_do);
;%        length_dl=spheriq_dist(lon_le,lat_le,lon_do,lat_do);
;%        length_ul=spheriq_dist(lon_le,lat_le,lon_up,lat_up);

;       angle_ur=90-acosd((-1*length_ur.^2+length_up.^2+length_ri.^2)./(2*length_up.*length_ri));
;       angle_dr=90-acosd((-1*length_dr.^2+length_do.^2+length_ri.^2)./(2*length_do.*length_ri));
;       angle_ul=90-acosd((-1*length_ul.^2+length_up.^2+length_le.^2)./(2*length_up.*length_le));
;       angle_dl=90-acosd((-1*length_dl.^2+length_do.^2+length_le.^2)./(2*length_do.*length_le));

;       orth(ii,jj)=max(abs([angle_ur, angle_ul, angle_dr, angle_dl]));

;       resolution(ii,jj)=mean([length_up.*length_ri, length_do.*length_le]);
;    end
;end
;clear *_up *_do *_ri *_le *_xo length*
;clear angle_* ii jj

;% Plot Orthogonality
;figure;clf;
;sethelv
;orth2=orth;orth2(orth2>=3)=3;
;pslice(lon_nodes, lat_nodes, orth2,[min(orth2(:)) max(orth2(:))],'Degrees')
;%pslice(lon_nodes, lat_nodes, resol2,[min(resol2(:)) max(resol2(:))],'M')
;shading flat %
;hold on
;plot3(lat_coast, lon_coast, ones(length(lat_coast),1), 'w')
;xlim([min(lon_nodes(:)) max(lon_nodes(:))])
;ylim([min(lat_nodes(:)) max(lat_nodes(:))])
;contour(lon_bathy,lat_bathy,-1*bathy, [20:20:200])
;dasp(-38.7)

;hold on
;for ii=[1:5:size(lon_nodes,1),size(lon_nodes,1)]
;plot(lon_nodes(ii,:),lat_nodes(ii,:),'k');
;end
;for jj=[1:5:size(lon_nodes,2),size(lon_nodes,2)]
;plot(lon_nodes(:,jj),lat_nodes(:,jj),'k');
;end
;hold on
;clear ii jj

;xlabel('Longitude')
;ylabel('Latitude')
;title('Orthogonality')



end
