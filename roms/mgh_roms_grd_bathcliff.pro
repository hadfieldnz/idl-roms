; svn $Id$
;+
; NAME:
;   MGH_ROMS_GRD_BATHCLIFF
;
; PURPOSE:
;   Given a ROMS grid file, load a "cliffened" version of bathymetry data into it.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   MGH_ROMS_GRD_BATHCLIFF, file_grid, HCLIFF=hcliff, HCLIP=hclip
;
; POSITIONAL PARAMETERS:
;   file_grid (input, scalar string)
;     The name of a ROMS grid file itno which bathymetry data are to
;     be written.
;
; SIDE EFFECTS:
;   Writes (if necessary overwrites) bathymetry data to the hraw
;   variable in the grid file. Unclipped data are written to
;   bathymetry version 0 and (if the HCLIP keyword has been specified)
;   clipped data are written to bathymetry version 1.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-07:
;     Written.
;-

pro MGH_ROMS_GRD_BATHCLIFF, file_grid, $
     DATASET=dataset, ETA_BAND=eta_band, HCLIFF=hcliff, HCLIP=hclip, $
     N_PRE_SMOOTH=n_pre_smooth, XI_BAND=xi_band

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Check file name etc

   if n_elements(file_grid) eq 0 then $
        message, 'Name for grid file not supplied'

   if ~ file_test(file_grid, /READ) then $
        message, 'Grid file cannot be read'

   if ~ file_test(file_grid, /WRITE) then $
        message, 'Grid file cannot be written to'

   ;; Specify defaults

   if n_elements(dataset) eq 0 then dataset = 'sandwell'

   if n_elements(eta_band) eq 0 then eta_band = [2,2]

   if n_elements(hcliff) eq 0 then hcliff = 2000

   ;; This is one of the MOMA standard depths
   if n_elements(hclip) eq 0 then hclip = 4977

   if n_elements(n_pre_smooth) eq 0 then n_pre_smooth = 5

   if n_elements(xi_band) eq 0 then xi_band = [2,2]

   ;; Open the grid file

   message, /INFORM, string(FORMAT='(%"Opening grid file %s")', file_grid)

   ogrd = obj_new('MGHncFile', file_grid, /MODIFY)

   ;; Add "hraw" bathymetry variable if necessary

   if ~ ogrd->HasVar('hraw') then begin
      message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'hraw')
      ogrd->VarAdd, 'hraw', ['xi_rho','eta_rho','bath']
   endif

   ;; Add land-mask variables

   if ~ ogrd->HasVar('mask_rho') then begin
      message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'mask_rho')
      ogrd->VarAdd, 'mask_rho', ['xi_rho','eta_rho']
   endif

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

   ;; Read grid data from file

   grid = {dim_rho: [ogrd->DimInfo('xi_rho', /DIMSIZE), $
                     ogrd->DimInfo('eta_rho', /DIMSIZE)], $
           lon_rho: ogrd->VarGet('lon_rho'), $
           lat_rho: ogrd->VarGet('lat_rho')}

   ;; Specify boundaries over which to get bathymetry data

   bound = {lon: mgh_minmax(grid.lon_rho)+[-2,2], lat: mgh_minmax(grid.lat_rho)+[-2,2]}

   ;; Get bathymetry data

   message, /INFORM, 'Getting bathymetry data'

   bath = mgh_bathymetry(BOUNDARIES=bound, DATASET=dataset)

   ;; Smooth bathymetry on original grid to avoid aliasing

   message, /INFORM, 'Pre-smoothing bathymetry'

   bath.height = smooth(bath.height, n_pre_smooth, /EDGE_TRUNCATE)

   ;; Interpolate bathymetry to the model grid.

   message, /INFORM, 'Interpolating to model grid'

   xx = mgh_locate(bath.lon, XOUT=grid.lon_rho[*])
   yy = mgh_locate(bath.lat, XOUT=grid.lat_rho[*])

   h = reform(interpolate(-bath.height, xx, yy, CUBIC=-0.5), $
              grid.dim_rho[0], grid.dim_rho[1])

   ;; Smooth depths in a band around the edge of the model grid

   if xi_band[0]+xi_band[1]+eta_band[0]+eta_band[1] gt 0 then begin

      ;; Western edge:
      if xi_band[0] gt 0 then begin
         b = xi_band[0]
         for j=0,grid.dim_rho[1]-1 do begin
            jj = (j > eta_band[0]) < (grid.dim_rho[1]-b-1)
            h[0:b-1,j] = h[b,jj]
         endfor
      endif

      ;; Eastern edge:
      if xi_band[1] gt 0 then begin
         b = xi_band[1]
         for j=0,grid.dim_rho[1]-1 do begin
            jj = (j > eta_band[0]) < (grid.dim_rho[1]-b-1)
            h[grid.dim_rho[0]-b:grid.dim_rho[0]-1,j] = h[grid.dim_rho[0]-b-1,jj]
         endfor
      endif

      ;; Southern edge:
      if eta_band[0] gt 0 then begin
         b = eta_band[0]
         for i=0,grid.dim_rho[0]-1 do h[i,0:b-1] = h[i,b]
      endif

      ;; Northern edge:
      if eta_band[1] gt 0 then begin
         b = eta_band[1]
         for i=0,grid.dim_rho[0]-1 do $
              h[i,grid.dim_rho[1]-b:grid.dim_rho[1]-1] = h[i,grid.dim_rho[1]-b-1]
      endif

   endif

   ;; Calculate the mask such that all points shallower than hcliff
   ;; are land, then set all depths hclip

   mask_rho = h gt hcliff

   replicate_inplace, h, hclip

   ;; Calculate other masks--see ROMS subroutine ana_mask

   mask_u = mask_rho[0:grid.dim_rho[0]-2,*] * mask_rho[1:grid.dim_rho[0]-1,*]

   mask_v = mask_rho[*,0:grid.dim_rho[1]-2] * mask_rho[*,1:grid.dim_rho[1]-1]

   mask_psi = mask_rho[0:grid.dim_rho[0]-2,0:grid.dim_rho[1]-2] * $
              mask_rho[0:grid.dim_rho[0]-2,1:grid.dim_rho[1]-1] * $
              mask_rho[1:grid.dim_rho[0]-1,0:grid.dim_rho[1]-2] * $
              mask_rho[1:grid.dim_rho[0]-1,1:grid.dim_rho[1]-1]

   ;; Write bathymetry to file

   bath_version = 0

   message, /INFORM, string(FORMAT='(%"Writing synthetic bathymetric data to ' + $
                            'bathymetry version %d")', bath_version)

   ogrd->VarPut, 'hraw', h, OFFSET=[0,0,bath_version]

   ;; Write masks to grid file

   message, /INFORM, 'Writing land mask to grid file'

   ogrd->VarPut, 'mask_rho', mask_rho
   ogrd->VarPut, 'mask_u', mask_u
   ogrd->VarPut, 'mask_v', mask_v
   ogrd->VarPut, 'mask_psi', mask_psi

   ;; Close the output file

   obj_destroy, ogrd

end

