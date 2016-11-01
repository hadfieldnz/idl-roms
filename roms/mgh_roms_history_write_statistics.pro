;+
; NAME:
;   MGH_ROMS_HISTORY_WRITE_STATISTICS
;
; PURPOSE:
;   For a ROMS history or similar file, calculate and save temporal
;   statistics data.
;
; CALLING SEQUENCE:
;   mgh_roms_history_write_statistics, history, file_out
;
; POSITIONAL PARAMETERS:
;   history
;     A reference to a ROMS history sequence object or a string array
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
;   file_out (input, scalar string)
;     The name of a netCDF file to be created (write-only)
;
; KEYWORD PARAMETERS:
;   TIME_RANGE (input, numeric 2-element vector)
;     Time interval (in days) over which to perform the analysis.
;     Default is [0,tmax], where tmin and tmax are minimum and
;     maximum times.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2015-11:
;     Written as MGH_ROMS_HISTORY_NCRA, adapted from MGH_ROMS_SCALAR_STATISTICS_TO_NCDF
;     which deals with Hslice data rather than full 3D data.
;   Mark Hadfield, 2015-12:
;     Renamed MGH_ROMS_HISTORY_WRITE_STATISTICS.
;   Mark Hadfield, 2015-12:
;     Added a time_range variable.
;-
function mgh_roms_history_write_statistics_rename, var, param

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   return, var+'_'+mgh_str_vanilla(param)

end

pro mgh_roms_history_write_statistics, history, file_out, $
     METHOD=method, PARAMETERS=parameters, TIME_RANGE=time_range, VARIABLES=variables

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

   ;; Process output file argument

   if n_elements(file_out) eq 0 then $
      message, 'Name for out file not supplied'

   ;; Process keywords

   if n_elements(parameters) eq 0 then parameters = 'mean'

   n_param = n_elements(parameters)

   if n_elements(method) eq 0 then begin
      ;; If all the specified parameters match an element in the following list, then data
      ;; can be retrieved record by record
      method = 0
      foreach p, parameters do begin
         found = 0B
         foreach pp, ['mean','min','max','fraction*'] do begin
            if strmatch(p, pp) then begin
               found = 1B
               break
            endif
         endforeach
         if ~ found then begin
            method = 1
            break
         endif
      endforeach
   endif

   n_param = n_elements(parameters)

   ;; Create output file

   message, /INFORM, string(FORMAT='(%"Creating output file %s")', file_out)

   oout = obj_new('MGHncFile', file_out, /CREATE, /CLOBBER)

   ;; Copy global attributes to the output file

   oout->AttCopy, ohis, /GLOBAL

   ;; Add or append a line to the history attribute

   fmt = '(%"Created by procedure mgh_roms_history_write_statistics at %s")'
   hstring = string(FORMAT=fmt, mgh_dt_string(mgh_dt_now()))
   if oout->HasAtt(/GLOBAL, 'history') then begin
      hstring = oout->AttGet(/GLOBAL, 'history') + string(10B) + hstring
   endif
   oout->AttAdd, /GLOBAL, 'history', temporary(hstring)

   ;; Copy all dimensions to the output file. The size of the unlimited
   ;; dimension will initially be zero.

   dim = ohis->DimNames()

   dim_unlim = dim[-1]

   if ~ ohis->DimInfo(dim_unlim, /IS_UNLIMITED) then $
      message, 'The final dimension in history is not unlimited'

   oout->DimCopy, ohis, dim

   ;; Create a dimension for storing the time range

   oout->DimAdd, 'two', 2

   ;; Compile a list of the variables in the file & sort them into three sets:
   ;;   - Static variables (no unlimited dimension)
   ;;   - The time coordinate variable
   ;;   - Record variables (with and unlimited dimension) other than
   ;;     the time coordinate

   var = ohis->VarNames(COUNT=n_var)

   unlim = bytarr(n_var)
   for i_var=0,n_var-1 do begin
      vdim = ohis->VarDimNames(var[i_var])
      unlim[i_var] = vdim[-1] eq dim_unlim
   endfor

   var_static = var[where(~ unlim, /NULL)]

   var_record = var[where(unlim, /NULL)]

   var_time = ''
   patt = ['ocean_time*','time*']
   for i=0,n_elements(patt) do begin
      l_match = where(strmatch(var_record, patt[i]), n_match)
      if n_match gt 0 then begin
         var_time = var_record[l_match[0]]
         break
      endif
   endfor
   if strlen(var_time) eq 0 then message, 'No time variable found'
   message, /INFORM, 'Time variable is '+var_time

   var_record = var_record[where(var_record ne var_time, /NULL)]

   ;; Determine variables for which statistics are to be calculated (by default,
   ;; all record variables)

   if n_elements(variables) eq 0 then variables = var_record

   ;; Check if all variables exist

   for i_var=0,n_elements(variables)-1 do begin
      vari = variables[i_var]
      if ~ ohis->HasVar(vari) then $
         message, 'Variable not found in history file: '+vari
   endfor

   ;; Copy definitions and attributes for static & time-coordinate variables

   oout->VarCopy, ohis, [var_static,var_time], /DEFINITION, /ATTRIBUTE, /UNPACK

   ;; Create a time range variable. This is a record variable, as I am keeping
   ;; in mind the possibility of extending this routine to store statistics for
   ;; more than one time range.

   time_units = ohis->AttGet(var_time, 'units')

   oout->VarAdd, 'time_range', ['two',dim_unlim]
   oout->AttAdd, 'time_range', 'long_name', 'time range'
   oout->AttAdd, 'time_range', 'units', mgh_str_subst(time_units, 'seconds', 'days')

   ;; For record variables, create new variables with the parameter names appended

   for i_param=0,n_param-1 do begin
      var_out = mgh_roms_history_write_statistics_rename(variables, parameters[i_param])
      oout->VarCopy, ohis, variables, RENAME=var_out, /DEFINITION, /UNPACK
   endfor

   ;; Copy data for all static variables

   oout->VarCopy, ohis, var_static, /DATA, /UNPACK

   ;; Get time and determine records over which statistics are to be calculated.

   time = ohis->VarGet(var_time, /AUTOSCALE)
   if n_elements(time_range) gt 0 then begin
      record_range = mgh_subset(time, time_range)
   endif
   if n_elements(record_range) eq 0 then begin
      n_time = ohis->DimInfo(dim_unlim, /DIMSIZE)
      record_range = [0,n_time-1]
   endif
   if n_elements(time_range) eq 0 then $
      time_range = time[record_range]

   msg = ['Processing data between records',strtrim(record_range,2),'times',mgh_format_float(time_range)]
   message, /INFORM, strjoin(temporary(msg), ' ')

   rr0 = record_range[0]
   rr1 = record_range[1]
   rrn = rr1-rr0+1

   time = time[rr0:rr1]

   ;; Write the time range and the mean time for the processed period

   oout->VarPut, 'time_range', time_range
   oout->VarPut, var_time, mgh_avg(time), /AUTOSCALE

   oout->Sync

   ;; Calculate and write statistics for all the record variables

   for i_var=0,n_elements(variables)-1 do begin

      vari = variables[i_var]

      message, /INFORM, 'Processing record variable '+vari

      mgh_tic

      ohis->VarInfo, vari, N_DIMS=n_dim, DIMENSIONS=dim

      value = ptrarr(n_param)

      case method of

         0: begin

            ;; Record-by-record.

            if n_dim eq 1 then begin
               offset = []
               count = []
            endif else begin
               offset = lonarr(n_dim-1)
               count = lonarr(n_dim-1)
            endelse

            ;; Set starting values
            for i_param=0,n_param-1 do begin
               pari = parameters[i_param]
               case 1B of
                  strmatch(pari, 'mean'): begin
                     value[i_param] = ptr_new(0)
                  end
                  strmatch(pari, 'min'):
                  strmatch(pari, 'max'):
                  strmatch(pari, 'fraction*'): begin
                     value[i_param] = ptr_new(0)
                  end
                  else: message, 'Parameter not supported for method 0: '+pari
               endcase
            endfor

            ;; Work through records accumulating satistics
            for r=0,rrn-1 do begin
               data = ohis->VarGet(vari, OFFSET=[offset,rr0+r], COUNT=[count,1], /AUTOSCALE)
               for i_param=0,n_param-1 do begin
                  pari = parameters[i_param]
                  case 1B of
                     strmatch(pari, 'mean'): begin
                        *value[i_param] += data/double(rrn)
                     end
                     strmatch(pari, 'min'): begin
                        if r eq 0 then begin
                           value[i_param] = ptr_new(data)
                        endif else begin
                           *value[i_param] = *value[i_param] < data
                        endelse
                     end
                     strmatch(pari, 'max'): begin
                        if r eq 0 then begin
                           value[i_param] = ptr_new(data)
                        endif else begin
                           *value[i_param] = *value[i_param] > data
                        endelse
                     end
                     strmatch(pari, 'fraction > *'): begin
                        pp = strsplit(pari, /EXTRACT)
                        threshold = n_elements(pp) gt 2 ? float(pp[2]) : 0
                        *value[i_param] += (data gt threshold)/double(rrn)
                     end
                  endcase
               endfor
            endfor

         end

         1: begin

            ;; Retrieve all records, stepping along the 2nd dimension (typically eta)
            ;; for 3D variables and the 2nd and 3rd (typically eta and s) for 4D variables

            case n_dim of
               1: begin
                  data = ohis->VarGet(vari, OFFSET=[rr0], COUNT=[rrn], /AUTOSCALE)
                  for i_param=0,n_param-1 do begin
                     pari = parameters[i_param]
                     case 1B of
                        strmatch(pari, 'mean'): begin
                           value[i_param] = ptr_new(mgh_avg(data))
                        end
                        strmatch(pari, 'min'): begin
                           value[i_param] = ptr_new(min(data))
                        end
                        strmatch(pari, 'max'): begin
                           value[i_param] = ptr_new(max(data))
                        end
                        strmatch(pari, 'fraction > *'): begin
                           pp = strsplit(pari, /EXTRACT)
                           threshold = n_elements(pp) gt 2 ? float(pp[2]) : 0
                           value[i_param] = total(data gt threshold)/n_elements(data)
                        end
                        strmatch(pari, 'percentile *'): begin
                           pp = strsplit(pari, /EXTRACT)
                           threshold = n_elements(pp) gt 1 ? float(pp[1]) : 50
                           value[i_param] = ptr_new(mgh_percentile(data, THRESHOLD=threshold))
                        end
                     endcase
                  endfor
                  mgh_undefine, data
               end
               2: begin
                  data = ohis->VarGet(vari, OFFSET=[0,rr0], COUNT=[0,rrn], /AUTOSCALE)
                  for i_param=0,n_param-1 do begin
                     pari = parameters[i_param]
                     case 1B of
                         strmatch(pari, 'mean'): begin
                           value[i_param] = ptr_new(mgh_avg(data, 2))
                        end
                        strmatch(pari, 'min'): begin
                           value[i_param] = ptr_new(min(data, DIMENSION=2))
                        end
                        strmatch(pari, 'max'): begin
                           value[i_param] = ptr_new(max(data, DIMENSION=2))
                        end
                        strmatch(pari, 'fraction > *'): begin
                           pp = strsplit(pari, /EXTRACT)
                           threshold = n_elements(pp) gt 2 ? float(pp[2]) : 0
                           value[i_param] = ptr_new(total(data gt threshold, 2)/rrn)
                        end
                        strmatch(pari, 'percentile *'): begin
                           pp = strsplit(pari, /EXTRACT)
                           threshold = n_elements(pp) gt 1 ? float(pp[1]) : 50
                           value[i_param] = ptr_new(cmapply('user:mgh_percentile', data, 2, FUNCTARGS={threshold: threshold}))
                        end
                     endcase
                  endfor
                  mgh_undefine, data
               end
               3: begin
                  for i_param=0,n_param-1 do $
                     value[i_param] = ptr_new(dblarr(dim[0:1]))
                  for j=0,dim[1]-1 do begin
                     data = ohis->VarGet(vari, OFFSET=[0,j,rr0], COUNT=[0,1,rrn], /AUTOSCALE)
                     for i_param=0,n_param-1 do begin
                        pari = parameters[i_param]
                        case 1B of
                           strmatch(pari, 'mean'): begin
                              (*value[i_param])[*,j] = mgh_avg(data, 3)
                           end
                           strmatch(pari, 'min'): begin
                              (*value[i_param])[*,j] = min(data, DIMENSION=3)
                           end
                           strmatch(pari, 'max'): begin
                              (*value[i_param])[*,j] = max(data, DIMENSION=3)
                           end
                           strmatch(pari, 'fraction > *'): begin
                              pp = strsplit(pari, /EXTRACT)
                              threshold = n_elements(pp) gt 2 ? float(pp[2]) : 0
                              (*value[i_param])[*,j] = total(data gt threshold, 3)/rrn
                           end
                           strmatch(pari, 'percentile *'): begin
                              pp = strsplit(pari, /EXTRACT)
                              threshold = n_elements(pp) gt 1 ? float(pp[1]) : 50
                              (*value[i_param])[*,j] = cmapply('user:mgh_percentile', data, 3, FUNCTARGS={threshold: threshold})
                           end
                        endcase
                     endfor
                     mgh_undefine, data
                  endfor
               end
               4: begin
                  for i_param=0,n_param-1 do $
                     value[i_param] = ptr_new(dblarr(dim[0:2]))
                  for k=0,dim[2]-1 do begin
                     for j=0,dim[1]-1 do begin
                        data = ohis->VarGet(vari, OFFSET=[0,j,k,rr0], COUNT=[0,1,1,rrn], /AUTOSCALE)
                        for i_param=0,n_param-1 do begin
                           pari = parameters[i_param]
                           case 1B of
                              strmatch(pari, 'mean'): begin
                                 (*value[i_param])[*,j,k] = mgh_avg(data, 4)
                              end
                              strmatch(pari, 'min'): begin
                                 (*value[i_param])[*,j,k] = min(data, DIMENSION=4)
                              end
                              strmatch(pari, 'max'): begin
                                 (*value[i_param])[*,j,k] = max(data, DIMENSION=4)
                              end
                              strmatch(pari, 'fraction > *'): begin
                                 pp = strsplit(pari, /EXTRACT)
                                 threshold = n_elements(pp) gt 2 ? float(pp[2]) : 0
                                 (*value[i_param])[*,j,k] = total(data gt threshold, 4)/rrn
                              end
                              strmatch(pari, 'percentile *'): begin
                                 pp = strsplit(pari, /EXTRACT)
                                 threshold = n_elements(pp) gt 1 ? float(pp[1]) : 50
                                 (*value[i_param])[*,j,k] = cmapply('user:mgh_percentile', data, 4, FUNCTARGS={threshold: threshold})
                              end
                           endcase
                        endfor
                     mgh_undefine, data
                     endfor
                  endfor
               end
            endcase

         end

      endcase

      for i_param=0,n_param-1 do begin
         var_out = mgh_roms_history_write_statistics_rename(vari, parameters[i_param])
         oout->VarPut, var_out, *value[i_param], /AUTOSCALE
      endfor

      ptr_free, value

      mgh_toc

      oout->Sync

   endfor

   ;; Close files

   obj_destroy, oout

end
