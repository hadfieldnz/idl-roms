;+
; NAME:
;   MGH_ROMS_GRD_FINISH
;
; PURPOSE:
;   Finish off a ROMS grid file by:
;    - Copying the appropriate bathymetry frame into the "h" variable
;    - Ensuring the u, v and psi masks are consistent with the rho
;      mask (if any)
;    - Creating position variables for the u, v and psi points.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   MGH_ROMS_GRD_FINISH, file_grd
;
; POSITIONAL PARAMETERS:
;   file_grd (input, sclar string)
;     The name of a ROMS grid file to be processed
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-01:
;     Written
;   Mark Hadfield, 2006-04:
;     Added staggered masking code.
;   Mark Hadfield, 2013-02:
;     Added staggered position variable code.
;   Mark Hadfield, 2015-12:
;     Added a call to the Sync method on the grid file after
;     variables have been added, to avoid subsequent reads returning
;     spurious data. It's odd that the need for this has not been
;     noticed previously.
;-
pro mgh_roms_grd_finish, file_grd, VERSION=version

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Check file name etc

   if n_elements(file_grd) eq 0 then $
        message, 'Name for grid file not supplied'

   if ~ file_test(file_grd, /READ) then $
        message, 'Grid file cannot be read'

   if ~ file_test(file_grd, /WRITE) then $
        message, 'Grid file cannot be written to'

   ;; Open the grid file

   message, /INFORM, 'Opening grid file '+strtrim(file_grd,2)

   ogrd = obj_new('MGHncFile', file_grd, /MODIFY)

   ;; Specify defaults

   if n_elements(version) eq 0 then version = ogrd->DimInfo('bath', /DIMSIZE) - 1

   ;; Define variables

   ;;; Bathymetry

   if ~ ogrd->HasVar('h') then $
        ogrd->VarAdd, 'h', ['xi_rho','eta_rho'], /DOUBLE

   ;; Masks

   if ogrd->HasVar('mask_rho') then begin
      message, /INFORM, 'Creating staggered mask variables'
      if ~ ogrd->HasVar('mask_u') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'mask_u')
         ogrd->VarAdd, 'mask_u', ['xi_u','eta_u']
      endif
      if ~ ogrd->HasVar('mask_v') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'mask_v')
         ogrd->VarAdd, 'mask_v', ['xi_v','eta_v']
      endif
      if ~ ogrd->HasVar('mask_psi') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'mask_psi')
         ogrd->VarAdd, 'mask_psi', ['xi_psi','eta_psi']
      endif
   endif

   ;; Location variables

   if ogrd->HasVar('lon_rho') then begin
      message, /INFORM, 'Creating staggered longitude data'
      if ~ ogrd->HasVar('lon_u') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'lon_u')
         ogrd->VarAdd, 'lon_u', ['xi_u','eta_u']
      endif
      if ~ ogrd->HasVar('lon_v') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'lon_v')
         ogrd->VarAdd, 'lon_v', ['xi_v','eta_v']
      endif
      if ~ ogrd->HasVar('lon_psi') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'lon_psi')
         ogrd->VarAdd, 'lon_psi', ['xi_psi','eta_psi']
      endif
   endif

   if ogrd->HasVar('lat_rho') then begin
      message, /INFORM, 'Creating staggered latitude data'
      if ~ ogrd->HasVar('lat_u') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'lat_u')
         ogrd->VarAdd, 'lat_u', ['xi_u','eta_u']
      endif
      if ~ ogrd->HasVar('lat_v') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'lat_v')
         ogrd->VarAdd, 'lat_v', ['xi_v','eta_v']
      endif
      if ~ ogrd->HasVar('lat_psi') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'lat_psi')
         ogrd->VarAdd, 'lat_psi', ['xi_psi','eta_psi']
      endif
   endif

   if ogrd->HasVar('x_rho') then begin
      message, /INFORM, 'Creating staggered x data'
      if ~ ogrd->HasVar('x_u') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'x_u')
         ogrd->VarAdd, 'x_u', ['xi_u','eta_u']
      endif
      if ~ ogrd->HasVar('x_v') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'x_v')
         ogrd->VarAdd, 'x_v', ['xi_v','eta_v']
      endif
      if ~ ogrd->HasVar('x_psi') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'x_psi')
         ogrd->VarAdd, 'x_psi', ['xi_psi','eta_psi']
      endif
   endif

   if ogrd->HasVar('y_rho') then begin
      message, /INFORM, 'Creating staggered y data'
      if ~ ogrd->HasVar('y_u') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'y_u')
         ogrd->VarAdd, 'y_u', ['xi_u','eta_u']
      endif
      if ~ ogrd->HasVar('y_v') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'y_v')
         ogrd->VarAdd, 'y_v', ['xi_v','eta_v']
      endif
      if ~ ogrd->HasVar('y_psi') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'y_psi')
         ogrd->VarAdd, 'y_psi', ['xi_psi','eta_psi']
      endif
   endif

   ogrd->Sync

   ;; Write data

   ;;; Bathymetry

   message, /INFORM, 'Copying bathymetry version '+strtrim(version,2)+' to "h" variable'

   hbath = ogrd->VarGet('hraw', OFFSET=[0,0,version], COUNT=[0,0,1])

   ogrd->VarPut, 'h', reform(hbath)

   ;; Masks

   if ogrd->HasVar('mask_rho') then begin

      ;; Get rho mask and calculate others--see ROMS subroutine ana_mask

      mask = ogrd->VarGet('mask_rho')

      mask_u = mask[0:-2,*]*mask[1:-1,*]
      mask_v = mask[*,0:-2]*mask[*,1:-1]

      mask_psi = mask[0:-2,0:-2]*mask[0:-2,1:-1]*mask[1:-1,0:-2]*mask[1:-1,1:-1]

      ogrd->VarPut, 'mask_u', mask_u
      ogrd->VarPut, 'mask_v', mask_v
      ogrd->VarPut, 'mask_psi', mask_psi

   endif

   ;; Locations

   if ogrd->HasVar('lon_rho') then begin

      lon = ogrd->VarGet('lon_rho')

      ogrd->VarPut, 'lon_u', mgh_stagger(lon, DELTA=[-1,0])
      ogrd->VarPut, 'lon_v', mgh_stagger(lon, DELTA=[0,-1])
      ogrd->VarPut, 'lon_psi', mgh_stagger(lon, DELTA=[-1,-1])

   endif

   if ogrd->HasVar('lat_rho') then begin

      lat = ogrd->VarGet('lat_rho')

      ogrd->VarPut, 'lat_u', mgh_stagger(lat, DELTA=[-1,0])
      ogrd->VarPut, 'lat_v', mgh_stagger(lat, DELTA=[0,-1])
      ogrd->VarPut, 'lat_psi', mgh_stagger(lat, DELTA=[-1,-1])

   endif

   if ogrd->HasVar('x_rho') then begin

      x = ogrd->VarGet('x_rho')

      ogrd->VarPut, 'x_u', mgh_stagger(x, DELTA=[-1,0])
      ogrd->VarPut, 'x_v', mgh_stagger(x, DELTA=[0,-1])
      ogrd->VarPut, 'x_psi', mgh_stagger(x, DELTA=[-1,-1])

   endif

   if ogrd->HasVar('y_rho') then begin

      y = ogrd->VarGet('y_rho')

      ogrd->VarPut, 'y_u', mgh_stagger(y, DELTA=[-1,0])
      ogrd->VarPut, 'y_v', mgh_stagger(y, DELTA=[0,-1])
      ogrd->VarPut, 'y_psi', mgh_stagger(y, DELTA=[-1,-1])

   endif

   ;; Close the output file

   obj_destroy, ogrd

end

