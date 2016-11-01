; svn $Id$
;+
; NAME:
;   MGH_ROMS_PSI_TO_NCDF
;
; PURPOSE:
;   Get stream function data from a ROMS file and save it to a netCDF file.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   MGH_ROMS_TPORT, history, file_out
;
; POSITIONAL PARAMETERS:
;   history
;     A reference to a ROMS history sequence object or a string arrary
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
;   file_out (input, scalar string)
;     The name of a netCDF file to be created (write-only)
;
; KEYWORD PARAMETERS:
;   LOC0 (input, numeric 2-element vector)
;     Location (lon,lat) at which stream function is set to
;     zero. Default is [173.0,-42.0], in South Island.
;
;   TIME_RANGE (input, numeric 2-element vector)
;     Range of times (in days) over which values are to be averaged
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-04:
;     Written.
;-

pro MGH_ROMS_PSI_TO_NCDF, history, file_out, LOC0=loc0, TIME_RANGE=time_range

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process history argument.

   case size(history, /TNAME) of
      'STRING': begin
         ohis = obj_new('MGHromsHistory', history)
         history_file = ohis
         history_destroy = 1
      end
      'OBJREF': begin
         ohis = history
         history_file = history
         history_destroy = 0
      end
      else: $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', history
   endcase

   ;; Process output file argument

   if n_elements(file_out) eq 0 then $
        message, 'Name for out file not supplied'

   ;; Process keywords

   if n_elements(loc0) eq 0 then loc0 = [173.0,-42.0]

   if n_elements(time_range) eq 0 then time_range = [1200,2000]

   ;; Get stream function data

   data = mgh_roms_barotropic(ohis, LOC0=loc0, TIME_RANGE=time_range)

   ;; Create output file

   message, /INFORM, string(FORMAT='(%"Creating output file %s")', file_out)

   oout = obj_new('MGHncFile', file_out, /CREATE, /CLOBBER)

   fmt = '(%"Stream function data from ROMS dataset %s")'
   oout->AttAdd, /GLOBAL, 'long_name', string(mgh_get_property(ohis, /NAME), FORMAT=fmt)

   ;; Create dimensions

   dim = size(data.psi, /DIMENSIONS)

   oout->DimAdd, 'xi', dim[0]
   oout->DimAdd, 'eta', dim[1]

   ;; Can the grid be represented by 1D variables in lon-lat space?

   lonlat1D = mgh_arr_ishom(data.lon, 2) && mgh_arr_ishom(data.lat, 1)

   case lonlat1D of
      0: begin
         dim_lon = ['xi','eta']
         dim_lat = ['xi','eta']
         lon = data.lon
         lat = data.lat
      end
      1: begin
         dim_lon = 'xi'
         dim_lat = 'eta'
         lon = reform(data.lon[*,0])
         lat = reform(data.lat[0,*])
      end
   endcase

   ;; Create variables

   oout->VarAdd, 'lon', dim_lon
   oout->AttAdd, 'lon', 'long_name', 'longitude'
   oout->AttAdd, 'lon', 'units', 'degrees E'

   oout->VarAdd, 'lat', dim_lat
   oout->AttAdd, 'lat', 'long_name', 'latitude'
   oout->AttAdd, 'lat', 'units', 'degrees'

   oout->VarAdd, 'psi', ['xi','eta']
   oout->AttAdd, 'psi', 'long_name', 'stream function'
   oout->AttAdd, 'psi', 'units', 'm^3 s^-1'

   ;; Write data

   oout->VarPut, 'lon', lon
   oout->VarPut, 'lat', lat

   oout->VarPut, 'psi', data.psi

   ;; Clean up

   obj_destroy, oout

   if history_destroy then obj_destroy, ohis

end

