;+
; NAME:
;   MGH_ROMS_GRD_NEW
;
; PURPOSE:
;   Generate a new ROMS grid file in lon-lat coordinates.
;
; CALLING SEQUENCE:
;   MGH_ROMS_GRD_NEW, file_grd
;
; POSITIONAL PARAMETERS:
;   file_grd (input, scalar string)
;     The name of a ROMS grid file to be generated (write-only)
;
; SIDE EFFECTS:
;   Writes (if necessary overwrites) the grid file.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2002-12:
;     Written
;   Mark Hadfield, 2003-09:
;     - For a regular grid specified via the BOUNDARIES & GRID_SPACING
;       parameters, the boundaries now specify the limits of the psi
;       grid, rather than the rho grid as previously. This increases
;       the size of the grid by one in each direction.
;     - Function MGH_ROMS_XY_TO_METRIC, called by this routine, now
;       calculates grid parameters using a spherical earth
;       rather than an elliptical one.
;-
pro mgh_roms_grd_new, file_grd, $
     BOUNDARIES=boundaries, GRID_SPACING=grid_spacing, $
     LON_RHO=lon_rho, LAT_RHO=lat_rho

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(file_grd) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'file_grd'

   if n_elements(file_grd) gt 1 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'file_grd'

   message, /INFORM, 'Creating grid file '+strtrim(file_grd,2)

   lonlat = (n_elements(lon_rho) gt 0) && (n_elements(lat_rho) gt 0)

   case lonlat of

      0: begin

         message, /INFORM, 'Longitude & latitude not supplied so creating ' + $
                  'rectangular grid'

         ;; Default is a 1/6 deg NZ-region grid.

         if n_elements(grid_spacing) eq 0 then $
              grid_spacing = 1/6.D0
         if n_elements(grid_spacing) eq 1 then $
              grid_spacing = replicate(grid_spacing, 2)

         if n_elements(boundaries) eq 0 then $
              boundaries = {lon: [155,197], lat: [-60,-25]}

         lon = mgh_range(boundaries.lon[0]-0.5D0*grid_spacing[0], $
                         boundaries.lon[1]+0.5D0*grid_spacing[0], $
                         STRIDE=grid_spacing[0])
         lat = mgh_range(boundaries.lat[0]-0.5D0*grid_spacing[1], $
                         boundaries.lat[1]+0.5D0*grid_spacing[1], $
                         STRIDE=grid_spacing[1])

         n_rho = [n_elements(lon),n_elements(lat)]

         lon_rho = mgh_inflate(n_rho, lon, 1)
         lat_rho = mgh_inflate(n_rho, lat, 2)

         mgh_undefine, lon, lat

      end

      1: begin

         message, /INFORM, 'Longitude & latitude specified explicitly'

         if size(lon_rho, /N_DIMENSIONS) ne 2 then $
              message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'lon_rho'

         if size(lat_rho, /N_DIMENSIONS) ne 2 then $
              message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'lat_rho'

         n_rho = size(lon_rho, /DIMENSIONS)

         if max(size(lat_rho, /DIMENSIONS) ne n_rho) then $
              message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', $
                       'lon_rho', 'lat_rho'

      end

   endcase

   mgh_roms_xy_to_metric, lon_rho, lat_rho, $
        LONLAT=1, ANGLE=angle, PM=pm, PN=pn, DNDX=dndx, DMDE=dmde

   ;; Specify Coriolis force. Omega is angular velocity of Earth's
   ;; rotation, assuming period is 1 sidereal day = 23 hr 56 m 4.1 s, see
   ;;   http://en.wikipedia.org/wiki/Coriolis_parameter
   omega = 7.2921151D-5

   fcor = 2*omega*sin(lat_rho*!dpi/180)

   ogrd = obj_new('MGHncFile', file_grd, /CREATE, /CLOBBER)

   ogrd->AttAdd, /GLOBAL, 'type', "ROMS grid file"

   fmt = '(%"Created by procedure mgh_roms_grd_new at %s")'
   ogrd->AttAdd, /GLOBAL, 'history', $
                 string(FORMAT=fmt, mgh_dt_string(mgh_dt_now()))

   ;; Define dimensions

   ogrd->DimAdd, 'xi_rho', n_rho[0]
   ogrd->DimAdd, 'xi_u', n_rho[0]-1
   ogrd->DimAdd, 'xi_v', n_rho[0]
   ogrd->DimAdd, 'xi_psi', n_rho[0]-1

   ogrd->DimAdd, 'eta_rho', n_rho[1]
   ogrd->DimAdd, 'eta_u', n_rho[1]
   ogrd->DimAdd, 'eta_v', n_rho[1]-1
   ogrd->DimAdd, 'eta_psi', n_rho[1]-1

   ogrd->DimAdd, 'bath'

   ;; Define variables

   ogrd->VarAdd, 'spherical', /CHAR

   ogrd->VarAdd, 'xl'
   ogrd->VarAdd, 'el'

   ogrd->VarAdd, 'pm', ['xi_rho','eta_rho'], /DOUBLE
   ogrd->VarAdd, 'pn', ['xi_rho','eta_rho'], /DOUBLE

   ogrd->VarAdd, 'dmde', ['xi_rho','eta_rho'], /DOUBLE
   ogrd->VarAdd, 'dndx', ['xi_rho','eta_rho'], /DOUBLE

   ogrd->VarAdd, 'lon_rho', ['xi_rho','eta_rho'], /DOUBLE
   ogrd->VarAdd, 'lat_rho', ['xi_rho','eta_rho'], /DOUBLE

   ogrd->VarAdd, 'angle', ['xi_rho','eta_rho'], /DOUBLE

   ogrd->VarAdd, 'f', ['xi_rho','eta_rho'], /DOUBLE

   ;; Store data

   ogrd->VarPut, 'xl', 0.
   ogrd->VarPut, 'el', 0.

   ogrd->VarPut, 'spherical', 'T'

   ogrd->VarPut, 'pm', pm
   ogrd->VarPut, 'pn', pn

   ogrd->VarPut, 'dmde', dmde
   ogrd->VarPut, 'dndx', dndx

   ogrd->VarPut, 'lon_rho', lon_rho
   ogrd->VarPut, 'lat_rho', lat_rho

   ogrd->VarPut, 'angle', angle

   ogrd->VarPut, 'f', fcor

   ;; Close files

   obj_destroy, ogrd

end
