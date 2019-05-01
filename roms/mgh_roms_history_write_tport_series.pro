;+
; NAME:
;   MGH_ROMS_HISTORY_WRITE_TPORT_SERIES
;
; PURPOSE:
;   For a ROMS history or similar file, calculate time series of volume
;   transport and write them to a netCDF file.
;
; CALLING SEQUENCE:
;   mgh_roms_history_write_tport_series, history, section
;
; POSITIONAL PARAMETERS:
;   history
;     A reference to a ROMS history sequence object or a string array
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
;   section (input, structure vector)
;     The section(s) for which the transport is to be calculated. The
;     structure must have the tags "lon" and "lat" (both 2-element vectors)
;     and optionally "name" (string scalar). The routine constructs linear
;     Xslices (as defined in the MGHromsHistory class) and calculates
;     time series of transport for each one.
;
; KEYWORD PARAMETERS:
;   FILE_OUT (input, scalar string)
;     The name of a netCDF file to be created (write-only)
;
;   RECALC (input, switch)
;     Set this keyword to recalculate the transport time series..
;
;   TIME_RANGE (input, numeric 2-element vector)
;     Time interval (in days) over which to perform the analysis.
;     Default is [0,tmax], where tmin and tmax are minimum and
;     maximum times.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2017-08:
;     Written.
;   Mark Hadfield, 2019-04:
;     - Fixed bug: the UNPACK keyword should not have been set when
;       copying the defintions for the ocean_time variable.
;     - Added VAR_UBAR/VBAR/ZETA keywords, passed to
;       mgh_roms_series_xslice_tport.
;-
pro mgh_roms_history_write_tport_series, history, section, file_out, $
     RECALC=recalc, TIME_RANGE=time_range, $
     VAR_UBAR=var_ubar, VAR_vbar=var_vbar, VAR_ZETA=var_zeta

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process history argument.

   case size(history, /TNAME) of
      'STRING': begin
         ohis = obj_new('MGHromsHistory', history)
      end
      'OBJREF': begin
         ohis = history
      end
      else: $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', history
   endcase

   ;; Process section argument

   if n_elements(section) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', section

   n_section = n_elements(section)

   ;; Process output file argument

   if n_elements(file_out) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', file_out

   ;; Work through sections, calculating transport. Save the time vector from the first
   ;; section and check the others agree.

   data = list()

   foreach sec,section,i_sec do begin
      data->Add, mgh_roms_series_xslice_tport(ohis, TYPE='linear', /LONLAT, $
         VERTX=sec.lon, VERTY=sec.lat, TIME_RANGE=time_range, $
         VAR_UBAR=var_ubar, VAR_vbar=var_vbar, VAR_ZETA=var_zeta, RECALC=recalc)
      if i_sec eq 0 then begin
         time = data[0].time
      endif else begin
         if ~ array_equal(data[i_sec].time, time) then message, 'Time mismatch'
      endelse
   endforeach

   ;; Create output file. NetCDF4 format allows string variables.

   message, /INFORM, string(FORMAT='(%"Creating output file %s")', file_out)

   oout = obj_new('MGHncFile', file_out, /CREATE, /CLOBBER, /NETCDF4_FORMAT)

   ;; Construct & write an appropriate history attribute, appending it to any
   ;; existing such attribute

   fmt = '(%"Created by procedure mgh_roms_history_write_tport_series at %s")'
   oout->AttAdd, /GLOBAL, 'history', string(FORMAT=fmt, mgh_dt_string(mgh_dt_now()))

   ;; Create dimensions

   oout->DimAdd, 'vertex', 2
   oout->DimAdd, 'section', n_section
   oout->DimAdd, 'ocean_time'

   ;; Copy definitions and attributes for the time-coordinate variable

   oout->VarCopy, ohis, 'ocean_time', /DEFINITION, /ATTRIBUTE

   ;; Set up other variables

   oout->VarAdd, 'section', ['section'], /STRING
   oout->VarAdd, 'lon', ['vertex','section']
   oout->VarAdd, 'lat', ['vertex','section']

   oout->VarAdd, 'transport', ['section','ocean_time']

   ;; Write section data

   foreach sec,section,i_sec do begin
      oout->VarPut, 'section', sec.name, OFFSET=[i_sec], COUNT=[1]
      oout->VarPut, 'lon', sec.lon, OFFSET=[0,i_sec], COUNT=[2,1]
      oout->VarPut, 'lat', sec.lat, OFFSET=[0,i_sec], COUNT=[2,1]
   endforeach

   ;; Write time-dependent data

   oout->VarPut, 'ocean_time', time

   foreach sec,section,i_sec do begin
      tport = data[i_sec].tport[-1,*]
      oout->VarPut, 'transport', tport, OFFSET=[i_sec,0], COUNT=[1,n_elements(time)]
   endforeach

   obj_destroy, oout

end
