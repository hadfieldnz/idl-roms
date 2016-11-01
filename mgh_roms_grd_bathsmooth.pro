; svn $Id$
;+
; NAME:
;   MGH_ROMS_GRD_BATHSMOOTH
;
; PURPOSE:
;   Smooth the bathymetry data in a ROMS grid file
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   MGH_ROMS_GRD_BATHSMOOTH, file_grid
;
; POSITIONAL PARAMETERS:
;   file_grid (input, sclar string)
;     The name of a ROMS grid file into which bathymetry data are to be written.
;
; SIDE EFFECTS:
;   Data are read from and written to the hraw variable.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2001-01:
;       Written as MGH_ROMS_SMOOTH_BATHYMETRY.
;   Mark Hadfield, 2003-01:
;       Renamed MGH_ROMS_GRD_BATHSMOOTH.
;-
pro mgh_roms_grd_bathsmooth, file_grid, $
     BATHSOAP_KEYWORDS=bathsoap_keywords, BATHSUDS_KEYWORDS=bathsuds_keywords, $
     READ_VERSION=read_version, WRITE_VERSION=write_version

   compile_opt DEFINT32
   compile_opt STRICTARR

   ;; Check file name etc

   if n_elements(file_grid) eq 0 then $
        message, 'Name for grid file not supplied'

   if not file_test(file_grid, /READ) then $
        message, 'Grid file cannot be read'

   if not file_test(file_grid, /WRITE) then $
        message, 'Grid file cannot be written to'

   ;; Open the grid file

   message, /INFORM, 'Opening grid file '+strtrim(file_grid,2)

   ogrid = obj_new('MGHncFile', file_grid, /MODIFY)

   ;; Specify defaults

   if n_elements(read_version) eq 0 then $
        read_version = ogrid->DimInfo('bath', /DIMSIZE) - 1

   if n_elements(write_version) eq 0 then $
        write_version = read_version + 1

   ;; Read grid data from file

   grid = {num_xi_rho: ogrid->DimInfo('xi_rho', /DIMSIZE), $
           num_eta_rho: ogrid->DimInfo('eta_rho', /DIMSIZE)}

   ;; Read bathymetry data from file

   message, /INFORM, string(FORMAT='(%"Reading data from bathymetry ' + $
                            'version %d")', read_version)

   hbath = reform(ogrid->VarGet('hraw', OFFSET=[0,0,read_version], COUNT=[0,0,1]))

   message, /INFORM, string(FORMAT='(%"Min & Max r values are %f %f")', $
                            mgh_minmax(mgh_roms_r_value(hbath)))

   ;; Apply bathsuds smoother

   message, /INFORM, 'Applying the bathsuds smoother'

   hbath = mgh_roms_bathsuds(hbath, N_PASSES=10, THRESHOLD=0.15, $
                             _STRICT_EXTRA=bathsuds_keywords)

   message, /INFORM, string(FORMAT='(%"Min & Max r values are %f %f")', $
                            mgh_minmax(mgh_roms_r_value(hbath)))

   ;; Apply bathsoap smoother.
   message, /INFORM, 'Applying the bathsoap smoother'

   hbath = mgh_roms_bathsuds(hbath, N_PASSES=1, THRESHOLD=0, $
                             _STRICT_EXTRA=bathsoap_keywords)

   message, /INFORM, string(FORMAT='(%"Min & Max r values are %f %f")', $
                            mgh_minmax(mgh_roms_r_value(hbath)))

   ;; Write smoothed bathymetry data

   message, /INFORM, string(FORMAT='(%"Writing smoothed data to bathymetry ' + $
                            'version %d")', write_version)

   ogrid->VarPut, 'hraw', hbath, OFFSET=[0,0,write_version]

   ;; Close the output file

   obj_destroy, ogrid

end

