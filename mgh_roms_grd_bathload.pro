;+
; NAME:
;   MGH_ROMS_GRD_BATHLOAD
;
; PURPOSE:
;   Given a ROMS grid file, load bathymetry data into it.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   MGH_ROMS_GRD_BATHLOAD, file_grid
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
;   Mark Hadfield, 2000-11:
;     Written as MGH_ROMS_LOAD_BATHYMETRY for Chatham Rise work.
;   Mark Hadfield, 2000-12:
;     Added ETA_BAND and XI_BAND keywords.
;   Mark Hadfield, 2001-01:
;     The procedure originally modified a couple of the grid file's
;     attributes ("history" & "CPP-options") but IDL has trouble with
;     moving the file between DEFINE and DATA modes so I have omitted
;     this code.
;   Mark Hadfield, 2002-12:
;     Renamed MGH_ROMS_GRD_BATHLOAD for RANZOM work.
;   Mark Hadfield, 2003-06:
;     Writing of land mask now controlled by a keyword.
;   Mark Hadfield, 2003-07:
;     Pre-smoothing now done by IDL SMOOTH function instead of
;     repeated applications of MGH_SHAPIRO. I have no idea what
;     difference this makes.
;   Mark Hadfield, 2004-11:
;     Now supports optional loading of user-specified bathymetry via
;     BATHYMETRY keyword.
;   Mark Hadfield, 2006-04:
;     Calculation of velocity and stream function masks deleted from
;     this routine, now handled by MGH_ROMS_GRD_FINISH.
;   Mark Hadfield, 2008-02:
;     Name of keyword to control masking changed to CALCULATE_MASK.
;   Mark Hadfield, 2008-10:
;     N_PRE_SMOOTH now defaults to zero. The previous default was 11!
;   Mark Hadfield, 2008-10:
;     N_PRE_SMOOTH now defaults to zero. The previous default was 11!
;-
pro mgh_roms_grd_bathload, file_grid, $
     BATHYMETRY=bathymetry, CALCULATE_MASK=calculate_mask, $
     DATASET=dataset, HCLIP=hclip, $
     N_PRE_SMOOTH=n_pre_smooth, $
     XI_BAND=xi_band, ETA_BAND=eta_band

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

   if n_elements(dataset) eq 0 then dataset = 'ETOPO2'

   if n_elements(eta_band) eq 0 then eta_band = [2,2]

   if n_elements(calculate_mask) eq 0 then calculate_mask = 1B

   if n_elements(n_pre_smooth) eq 0 then n_pre_smooth = 0

   if n_elements(xi_band) eq 0 then xi_band = [2,2]

   ;; Open the grid file

   message, /INFORM, string(FORMAT='(%"Opening grid file %s")', file_grid)

   ogrd = obj_new('MGHncFile', file_grid, /MODIFY)

   ;; Add "hraw" bathymetry variable if necessary

   if ~ ogrd->HasVar('hraw') then begin
      message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'hraw')
      ogrd->VarAdd, 'hraw', ['xi_rho','eta_rho','bath'], /DOUBLE
   endif

   ;; Add land mask if necessary

   if keyword_set(calculate_mask) then begin

      if ~ ogrd->HasVar('mask_rho') then begin
         message, /INFORM, string(FORMAT='(%"Adding variable %s")', 'mask_rho')
         ogrd->VarAdd, 'mask_rho', ['xi_rho','eta_rho']
      endif

   endif

   ;; Read grid data from file

   grid = {dim_rho: [ogrd->DimInfo('xi_rho', /DIMSIZE), $
                     ogrd->DimInfo('eta_rho', /DIMSIZE)], $
           lon_rho: ogrd->VarGet('lon_rho'), $
           lat_rho: ogrd->VarGet('lat_rho')}

   ;;Retrieve bathymetry data if necessary

   case 1B of

      n_elements(bathymetry) gt 0: begin

         message, /INFORM, 'Using user-supplied bathymetry'

         bath = bathymetry

      end

      else: begin

         ;; Get bathymetry data. Fill non-finite values if necessary

         bound = {lon: mgh_minmax(grid.lon_rho)+[-2,2], $
                  lat: mgh_minmax(grid.lat_rho)+[-2,2]}

         message, /INFORM, 'Getting bathymetry data from dataset '+dataset

         bath = mgh_bathymetry(BOUNDARIES=bound, DATASET=dataset)

         bad = where(~ finite(bath.height), n_bad)
         if n_bad gt 0 then bath.height[bad] = 10

      end

   endcase

   ;; What we do now depends on whether the bath variable is a structure
   ;; or an array

   case 1B of

      size(bath, /TYPE) eq 8: begin

         ;; The bath variable is a structure, assumed to contain tags
         ;; lon, lat and height

         ;; Optionally pre-smooth bathymetry

         if n_pre_smooth gt 0 then begin
            message, /INFORM, 'Pre-smoothing bathymetry'
            bath.height = smooth(bath.height, n_pre_smooth, /EDGE_TRUNCATE)
         endif

         ;; Locate the model grid points in the bathymetry grid, then
         ;; interpolate.

         message, /INFORM, 'Interpolating to model grid'

         case size(bath.lon, /N_DIMENSIONS) of
            1: begin
               ix = mgh_locate(bath.lon, XOUT=grid.lon_rho)
               iy = mgh_locate(bath.lat, XOUT=grid.lat_rho)
            end
            2: begin
               loc = mgh_locate2(bath.lon, bath.lat, $
                                 XOUT=grid.lon_rho, YOUT=grid.lat_rho)
               ix = reform(loc[0,*,*])
               iy = reform(loc[1,*,*])
               mgh_undefine, loc
            end
         endcase

         h = interpolate(-bath.height, ix, iy, CUBIC=-0.5)

      end

      else: begin

         ;; The bath variable is not a structure, so assume it is a
         ;; numeric array of the same dimensions as the model grid

         if ~ array_equal(size(bath, /DIMENSIONS), grid.dim_rho) then $
              message, 'Bathymetry data dimensions do not match model grid'

         ;; Optionally pre-smooth bathymetry

         if n_pre_smooth gt 0 then begin
            message, /INFORM, 'Pre-smoothing bathymetry'
            bath = smooth(bath, n_pre_smooth, /EDGE_TRUNCATE)
         endif

         h = -bath

      end

   endcase

   ;; Specify initial bathymetry version

   bath_version = 0

   ;; Write bathymetry to grid file

   message, /INFORM, string(FORMAT='(%"Writing raw bathymetric data to ' + $
                            'bathymetry version %d")', bath_version)

   ogrd->VarPut, 'hraw', h, OFFSET=[0,0,bath_version]

   bath_version += 1

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

      ;; Write bathymetry again

      message, /INFORM, string(FORMAT='(%"Writing edge-smoothed bathymetric data to ' + $
                               'bathymetry version %d")', bath_version)

      ogrd->VarPut, 'hraw', h, OFFSET=[0,0,bath_version]

      bath_version += 1

   endif

   if keyword_set(calculate_mask) then begin

      ;; Calculate mask for RHO points

      mask_rho = h gt 0

      ;; Write masks to grid file

      message, /INFORM, 'Writing land mask to grid file'

      ogrd->VarPut, 'mask_rho', mask_rho

   endif

   ;; Clip bathymetry

   if n_elements(hclip) eq 2 then begin

      h = (h > hclip[0]) < hclip[1]

      ;; Write bathymetry again

      message, /INFORM, 'Writing clipped bathymetric data to bathymetry version '+ $
               strtrim(bath_version,2)

      ogrd->VarPut, 'hraw', h, OFFSET=[0,0,bath_version]

      bath_version += 1

   endif

   ;; Close the output file

   obj_destroy, ogrd

end

