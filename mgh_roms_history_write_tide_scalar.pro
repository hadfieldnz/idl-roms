;+
; NAME:
;   MGH_ROMS_HISTORY_WRITE_TIDE_sCALAR
;
; PURPOSE:
;   For a ROMS history or similar file, calculate and save scalar tidal analysis data.
;
; CALLING SEQUENCE:
;   mgh_roms_history_write_tide_scalar, history, file_out
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
;   Mark Hadfield, 2015-12:
;     Written, merging functionality from MGH_ROMS_HISTORY_WRITE_ANNUAL,
;     MGH_ROMS_NAIVE_SCALAR_TO_NCDF and MGH_ROMS_TIDE_SCALAR_TO_NCDF.
;-
function mgh_roms_history_write_tide_scalar_cname, var, i_coeff, CONSTITUENTS=constituents

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Generate variable names appropriate for the coefficients
   ;; as implied by MGH_TIDE_SCALAR_ANALYSIS. Note that var can be a vector.

   case 1B of
      i_coeff eq 0: begin
         return, var+'_mean'
      end
      i_coeff gt 0: begin
         s0 = constituents[(i_coeff-1)/2]
         s1 = i_coeff eq 2*(i_coeff/2) ? 'sin' : 'cos'
         return, var+'_'+s0+'_'+s1
      end
   endcase

end

pro mgh_roms_history_write_tide_scalar, history, file_out, $
     CONSTITUENTS=constituents, TIME_RANGE=time_range, TIME_STAMPED=time_stamped, VARIABLES=variables

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

   if n_elements(constituents) eq 0 then constituents = ['M2','M4','M6']

   n_con = n_elements(constituents)

   ;; Create output file

   message, /INFORM, string(FORMAT='(%"Creating output file %s")', file_out)

   oout = obj_new('MGHncFile', file_out, /CREATE, /CLOBBER)

   ;; Copy global attributes to the output file

   oout->AttCopy, ohis, /GLOBAL

   ;; Add or append a line to the history attribute

   fmt = '(%"Created by procedure mgh_roms_history_write_tide_scalar at %s")'
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
      message, 'Last dimension in history is not unlimited'

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

   ;; For record variables, create new variables for the tidal coefficients

   for i_coeff=0,2*n_con do begin
      var_out = mgh_roms_history_write_tide_scalar_cname(variables, i_coeff, CONSTITUENTS=constituents)
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

   ;; Parse the time units. Since time data has been autoscaled above, we're
   ;; only interested in the offset.

   time_units = {scale: 1D/(24D*3600D), offset: 0D}
   if ohis->MGHncSequence::HasAtt(var_time, 'units') then $
      time_units = mgh_dt_units(ohis->AttGet(var_time, 'units'))
   time_ref = time_units.offset
   mgh_undefine, time_units

   ;; Determine whether a time-stamped analysis is to be carried out.

   if n_elements(time_stamped) eq 0 then $
      time_stamped = time_ref gt 0

   ;; If time-stamped, create an MGHtideParameters object to use in the
   ;; analysis; otherwise determine the frequencies associated with each
   ;; constituent

   if time_stamped then begin
      oprm = obj_new('MGHtideParameters', /DEFAULT_PARAMETERS)
   endif else begin
      n_con = n_elements(constituents)
      const = mgh_t_getconsts()
      l_con = lonarr(n_con)
      for c=0,n_con-1 do begin
         l_con[c] = where(strmatch(const.name, constituents[c], /FOLD_CASE))
         if l_con[c] lt 0 then $
            message, 'Could not find constituent: '+constituents[c]
      endfor
      ;; Convert to d^-1
      frequency = 24.0D*const.freq[l_con]
   endelse

   ;; Work through variables, calculating and writing coefficients

   for i_var=0,n_elements(variables)-1 do begin

      vari = variables[i_var]

      message, /INFORM, 'Processing record variable '+vari

      mgh_tic

      ohis->VarInfo, vari, N_DIMS=n_dim, DIMENSIONS=dim

      ;; Retrieve all records, stepping along penultimate (typically
      ;; vertical) dimension for 4D variables

      case n_dim of
         1: begin
            coeff = [!values.f_nan]
            data = ohis->VarGet(vari, OFFSET=[rr0], COUNT=[rrn], /AUTOSCALE)
            if min(finite(data)) eq 1 then begin
               if keyword_set(time_stamped) then begin
                  mgh_tide_scalar_analysis, $
                     oprm, time+time_ref, data, cc, CONSTITUENTS=constituents
               endif else begin
                  mgh_naive_scalar_analysis, $
                     time, data, cc, FREQUENCY=frequency
               endelse
               coeff = cc
            endif
            mgh_undefine, data
         end
         2: begin
            data = ohis->VarGet(vari, OFFSET=[0,rr0], COUNT=[0,rrn], /AUTOSCALE)
            coeff = replicate(!values.f_nan, [2*n_con+1,dim[0]])
            for i=0,dim[0]-1 do begin
               dd = reform(data[i,*])
               if min(finite(dd)) eq 1 then begin
                  if keyword_set(time_stamped) then begin
                     mgh_tide_scalar_analysis, $
                        oprm, time+time_ref, dd, cc, CONSTITUENTS=constituents
                  endif else begin
                     mgh_naive_scalar_analysis, $
                        time, dd, cc, FREQUENCY=frequency
                  endelse
                  coeff[*,i] = cc
               endif
            endfor
            mgh_undefine, data
         end
         3: begin
            coeff = replicate(!values.f_nan, [2*n_con+1,dim[0:1]])
            for j=0,dim[1]-1 do begin
               data = ohis->VarGet(vari, OFFSET=[0,j,rr0], COUNT=[0,1,rrn], /AUTOSCALE)
               for i=0,dim[0]-1 do begin
                  dd = reform(data[i,0,*])
                  if min(finite(dd)) eq 1 then begin
                     if keyword_set(time_stamped) then begin
                        mgh_tide_scalar_analysis, $
                           oprm, time+time_ref, dd, cc, CONSTITUENTS=constituents
                     endif else begin
                        mgh_naive_scalar_analysis, $
                           time, dd, cc, FREQUENCY=frequency
                     endelse
                     coeff[*,i,j] = cc
                  endif
               endfor
            endfor
            mgh_undefine, data
         end
         4: begin
            coeff = replicate(!values.f_nan, [2*n_con+1,dim[0:2]])
            for k=0,dim[2]-1 do begin
               for j=0,dim[1]-1 do begin
                  data = ohis->VarGet(vari, OFFSET=[0,j,k,rr0], COUNT=[0,1,1,rrn], /AUTOSCALE)
                  for i=0,dim[0]-1 do begin
                     dd = reform(data[i,0,0,*])
                     if min(finite(dd)) eq 1 then begin
                        if keyword_set(time_stamped) then begin
                           mgh_tide_scalar_analysis, $
                              oprm, time+time_ref, dd, cc, CONSTITUENTS=constituents
                        endif else begin
                           mgh_naive_scalar_analysis, $
                              time, dd, cc, FREQUENCY=frequency
                        endelse
                        coeff[*,i,j,k] = cc
                     endif
                  endfor
               endfor
               mgh_undefine, data
            endfor
         end
      endcase

      for i_coeff=0,2*n_con do begin
         var_out = mgh_roms_history_write_tide_scalar_cname(vari, i_coeff, CONSTITUENTS=constituents)
         case n_dim of
            1: oout->VarPut, var_out, coeff[i_coeff], /AUTOSCALE
            2: oout->VarPut, var_out, reform(coeff[i_coeff,*]), /AUTOSCALE
            3: oout->VarPut, var_out, reform(coeff[i_coeff,*,*]), /AUTOSCALE
            4: oout->VarPut, var_out, reform(coeff[i_coeff,*,*,*]), /AUTOSCALE
         end
      endfor

      mgh_toc

      oout->Sync

   endfor

   ;; Close files

   obj_destroy, oout

end
