;+
; NAME:
;   MGH_ROMS_HISTORY_WRITE_PPORT_SERIES
;
; PURPOSE:
;   For a ROMS history or similar file, calculate time series of volume
;   transport through one or more Pslices and write them to a netCDF file.
;
;   See also MGH_ROMS_HISTORY_WRITE_TPORT_SERIES.
;
; CALLING SEQUENCE:
;   mgh_roms_history_write_pport_series, history, slice
;
; POSITIONAL PARAMETERS:
;   history
;     A reference to a ROMS history sequence object or a string array
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
;   slice (input, a list of structures)
;     The P-slices for which the transport is to be calculated.
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
;   Mark Hadfield, 2019-11:
;     Written based on MGH_ROMS_HISTORY_WRITE_TPORT_SERIES.
;-
pro mgh_roms_history_write_pport_series, history, slice, file_out, $
     RECALC=recalc, TIME_RANGE=time_range, $
     VAR_UBAR=var_ubar, VAR_vbar=var_vbar, VAR_ZETA=var_zeta, USE_ZETA=use_zeta

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

   ;; Process slice argument

   if n_elements(slice) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', slice

   n_slice = n_elements(slice)

   ;; Process output file argument

   if n_elements(file_out) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', file_out

   ;; Work through sections, calculating transport. Save the time vector from the first
   ;; section and check the others agree.

   data = list()

   foreach s,slice,i do begin
      data->Add, mgh_roms_series_pslice_tport(ohis, s, $
         TIME_RANGE=time_range, RECALC=recalc, USE_ZETA=use_zeta, $
         VAR_UBAR=var_ubar, VAR_vbar=var_vbar, VAR_ZETA=var_zeta)
      if i eq 0 then begin
         time = data[0].time
      endif else begin
         if ~ array_equal(data[i].time, time) then message, 'Time mismatch'
      endelse
   endforeach

   ;; Create output file. NetCDF4 format allows string variables.

   message, /INFORM, string(FORMAT='(%"Creating output file %s")', file_out)

   oout = obj_new('MGHncFile', file_out, /CREATE, /CLOBBER, /NETCDF4_FORMAT)

   ;; Construct & write an appropriate history attribute, appending it to any
   ;; existing such attribute

   fmt = '(%"Created by procedure mgh_roms_history_write_pport_series at %s")'
   oout->AttAdd, /GLOBAL, 'history', string(FORMAT=fmt, mgh_dt_string(mgh_dt_now()))

   ;; Create dimensions

   oout->DimAdd, 'slice', n_slice
   oout->DimAdd, 'ocean_time'

   ;; Copy definitions and attributes for the time-coordinate variable

   oout->VarCopy, ohis, 'ocean_time', /DEFINITION, /ATTRIBUTE

   ;; Set up other variables

   oout->VarAdd, 'slice', ['slice'], /STRING

   oout->VarAdd, 'transport', ['slice','ocean_time']

   ;; Write section data

   foreach s,slice,i do begin
      oout->VarPut, 'slice', s.name, OFFSET=[i], COUNT=[1]
   endforeach

   ;; Write time-dependent data

   oout->VarPut, 'ocean_time', time

   foreach s,slice,i do begin
      tport = data[i].tport[-1,*]
      oout->VarPut, 'transport', tport, OFFSET=[i,0], COUNT=[1,n_elements(time)]
   endforeach

   obj_destroy, oout

end
