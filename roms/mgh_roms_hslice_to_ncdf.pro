; svn $Id$
;+
; NAME:
;   MGH_ROMS_HSLICE_TO_NCDF
;
; PURPOSE:
;   Get data from one or more Hslices from a ROMS file and write it to a netCDF file.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   MGH_ROMS_VAR3D_TO_NCDF, history, file_out
;
; POSITIONAL PARAMETERS:
;   history
;     A ROMS history sequence as an MGHromsHistory object or string.
;
;   file_out (input, scalar string)
;     The name of a netCDF file to be created (write-only)
;
; KEYWORD PARAMETERS:
;   DEPTHS (input, numeric scalar or vector)
;     A list of depths to be written
;
;   VARIABLE (input, string scalar)
;     The name of the variable to be extracted and written.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2008-12:
;     Written.
;   Mark Hadfield, 2009-10:
;     Removed calls to widget_event(/NOWAIT).
;-
pro mgh_roms_hslice_to_ncdf, history, file_out, $
     DEPTHS=depths, VARIABLE=variable

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process history argument.

   case size(history, /TNAME) of
      'STRING': begin
         ohis = obj_new('MGHromsHistory', history)
         history_file = ohis
         history_destroy = 1B
      end
      'OBJREF': begin
         ohis = history
         history_file = history
         history_destroy = 0B
      end
      else: $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', history
   endcase

   ;; Process output file argument

   if n_elements(file_out) eq 0 then $
        message, 'Name for out file not supplied'

   ;; Process keywords

   if n_elements(depths) eq 0 then depths = mgh_moma_grid([5,5505], [10,210])

   if n_elements(variable) eq 0 then variable = 'temp'

   ;; Set miscellaneous constants

   fill_value = mgh_ncdf_fill()

   ;; Get Hslice grid data

   grid = ohis->HsliceGrid(variable)

   ;; Create output file

   message, /INFORM, string(FORMAT='(%"Creating output file %s")', file_out)

   onc = obj_new('MGHncFile', file_out, /CREATE, /CLOBBER)

   fmt = '(%"Data from ROMS dataset %s")'
   onc->AttAdd, /GLOBAL, 'long_name', string(mgh_get_property(ohis, /NAME), FORMAT=fmt)

   ;; Create dimensions

   dim = size(grid.x, /DIMENSIONS)

   onc->DimAdd, 'xi', dim[0]
   onc->DimAdd, 'eta', dim[1]

   onc->DimAdd, 'depth', n_elements(depths)

   onc->DimAdd, 'time'

   ;; Create variables

   if grid.lonlat then begin
      x_var = 'lon'
      y_var = 'lat'
   endif else begin
      x_var = 'x'
      y_var = 'y'
   endelse

   onc->VarAdd, x_var, /DOUBLE, ['xi','eta']
   onc->VarAdd, y_var, /DOUBLE, ['xi','eta']

   onc->VarAdd, 'depth', ['depth']

   onc->VarAdd, variable, ['xi','eta','depth','time']
   onc->AttAdd, variable, '_FillValue', fill_value

   onc->VarAdd, 'time', ['time']
   onc->AttAdd, 'time', 'units', 'days'

   ;; Load grid & time data

   onc->VarPut, x_var, grid.x
   onc->VarPut, y_var, grid.y

   onc->VarPut, 'depth', depths

   onc->VarPut, 'time', ohis->VarGet('ocean_time', /AUTOSCALE)

   ;; Work through records in input file, extracting and writing hslice data

   ohis->GetProperty, N_RECORDS=n_rec

   for r=0,n_rec-1 do begin

      hdata = ohis->HsliceData(variable, GRID=grid, DEPTHS=depths, RECORD=r)

      l_miss = where(~ finite(hdata), n_miss)
      if n_miss gt 0 then hdata[l_miss] = fill_value

      onc->VarPut, variable, temporary(hdata), OFFSET=[0,0,0,r]

   endfor

   obj_destroy, onc

   if history_destroy then obj_destroy, ohis

end

