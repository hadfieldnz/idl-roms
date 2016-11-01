;+
; FUNCTION NAME:
;   MGH_ROMS_HSLICE_STATISTICS
;
; PURPOSE:
;   This function calculates and returns statistics (mean, percentiles,
;   min/max and others) from a series of Hslices extracted from a ROMS
;   history file (MGHromsHistory) object.
;
; CALLING SEQUENCE
;   result = mgh_roms_hslice_statistics(ofile, variable)
;
; RETURN VALUE:
;   The function returns the results
;   a structure.
;
; POSITIONAL ARGUMENTS:
;   ofile (input, object)
;     A reference to an MGHromsHistory object.
;
;   variable (input)
;     The name of a variable recognised by the MGHromsHistory object.
;
; KEYWORD ARGUMENTS:
;   DATA_MULTIPLIER (input, numeric scalar)
;     Number by which data values are multiplied before statistical processing.
;     The default depends on the variable and is calculated by MGH_ROMS_RESOLVE_DATA
;
;   DATA_TRANSFORMATION (input, numeric string)
;     Function applied to data values before statistical processing.
;
;   DEPTH
;     Set this keyword to a scalar numeric value to specify the
;     depth of a z surface on which data are to be plotted. This
;     keyword should be specified only for variables having a depth
;     coordinate and it cannot be used together with LEVEL or SIGMA.
;
;   LAYER
;     Set this keyword to a scalar integer to specify the bed layer to
;     be plotted.  This keyword should be specified only for variables
;     having a bed-layer dimension.
;
;   LEVEL
;     Set this keyword to a scalar integer to specify the
;     s-coordinate level to be plotted.  This keyword should be
;     specified only for variables having a depth coordinate and it
;     cannot be used together with DEPTH or SIGMA.
;
;   RECORD_RANGE (input, numeric 2-element vector)
;     Range of record numbers (0-based indices in the time dimension) over
;     which statistics are to be calculated. The record range may also be
;     specified indirectly via the TIME_RANGE keyword.
;
;   SIGMA
;     Set this keyword to a scalar numeric value to specify the
;     sigma values of a constant-sigma surface on which data are to be
;     plotted. This keyword should be specified only for variables having
;     a depth coordinate and it cannot be used together with DEPTH or LEVEL.
;
;   TIME_RANGE (input, 2-element numeric)
;     Time range (days) over which statistics are to be calculated. The
;     record range may also be specified directly via the RECORD_RANgE keyword.
;
;   XI_RANGE
;   ETA_RANGE
;     Use these keywords to display a subset of the domain. They are
;     interpreted relative to the rho grid.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2016-03:
;     Written, based on code extracted from mgh_roms_plot_hstats__define.
;-
function mgh_roms_hslice_statistics, ohis, variable, RECALC=recalc, $
     DATA_MULTIPLIER=data_multiplier, $
     DATA_TRANSFORMATION=data_transformation, $
     DEPTH=depth, LEVEL=level, SIGMA=sigma, $
     PARAMETER=parameter, $
     RECORD_RANGE=record_range, TIME_RANGE=time_range, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process variable name argument

   if n_elements(variable) eq 0 then variable = 'zeta'

   if ~ (isa(variable, 'STRING') || isa(variable, 'STRUCT'))  then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'variable'

   if n_elements(variable) ne 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'variable'

   ;; Specify statistical parameter to be calculated

   if n_elements(parameter) eq 0 then parameter = 'mean'

   ;; Process history file argument

   if n_elements(ohis) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'ohis'

   if n_elements(ohis) gt 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'ohis'

   if ~ obj_valid(ohis) gt 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_objref_bad', 'ohis'

   ;; Construct the temporary file name. Quite an operation!

   ohis->GetProperty, FILE_NAME=fname

   str_file = mgh_format_integer(mgh_hashcode(fname))
   str_var = mgh_format_integer(mgh_hashcode(variable))

   tmpfile = list('mgh_roms_hslice_statistics', str_file, str_var, parameter)
   if n_elements(depth) gt 0 then $
      tmpfile->Add, ['depth',mgh_format_float(depth)], /EXTRACT
   if n_elements(level) gt 0 then $
      tmpfile->Add, ['level',mgh_format_integer(level)], /EXTRACT
   if n_elements(sigma) gt 0 then $
      tmpfile->Add, ['sigma',mgh_format_float(sigma)], /EXTRACT
   if n_elements(xi_range) gt 0 then $
      tmpfile->Add, ['xr',mgh_format_integer(xi_range)], /EXTRACT
   if n_elements(eta_range) gt 0 then $
      tmpfile->Add, ['er',mgh_format_integer(eta_range)], /EXTRACT
   if n_elements(record_range) gt 0 then $
      tmpfile->Add, ['rr',mgh_format_integer(record_range)], /EXTRACT
   if n_elements(time_range) gt 0 then $
      tmpfile->Add, ['tr',mgh_format_float(time_range)], /EXTRACT
   if n_elements(data_multiplier) gt 0 then $
      tmpfile->Add, ['dm',mgh_format_float(data_multiplier)], /EXTRACT
   if n_elements(data_transformation) gt 0 then $
      tmpfile->Add, ['dt',data_transformation], /EXTRACT
   tmpfile = filepath(strjoin(tmpfile->ToArray(), '_')+'.idl_data', /TMP)

   if ~ file_test(tmpfile) then recalc = !true

   ;; Set ETA_RANGE and XI_RANGE relative to the rho grid.
   ;; The default is to show all  points.

   dim_rho = [ohis->DimInfo('xi_rho', /DIMSIZE),ohis->DimInfo('eta_rho', /DIMSIZE)]

   if n_elements(xi_range) eq 0 then xi_range = [0,dim_rho[0]-1]
   if n_elements(eta_range) eq 0 then eta_range = [0,dim_rho[1]-1]

   if xi_range[0] lt 0 then xi_range[0] += dim_rho[0]
   if xi_range[1] lt 0 then xi_range[1] += dim_rho[0]

   if eta_range[0] lt 0 then eta_range[0] += dim_rho[1]
   if eta_range[1] lt 0 then eta_range[1] += dim_rho[1]

   xra0 = xi_range[0]
   xra1 = xi_range[1]
   xran = xra1-xra0+1
   era0 = eta_range[0]
   era1 = eta_range[1]
   eran = era1-era0+1

   ;; Get grid data required for horizontal slice retrievals

   var_xi_range = xi_range
   var_eta_range = eta_range

   var_dims = ohis->VarDims(variable)

   vdh0 = var_dims.horizontal[0]
   if vdh0 eq 'xi_u' || vdh0 eq 'xi_psi' then $
        var_xi_range += [0,-1]

   vdh1 = var_dims.horizontal[1]
   if vdh1 eq 'eta_v' || vdh1 eq 'eta_psi' then $
        var_eta_range += [0,-1]

   grid = ohis->HsliceGrid(variable, ETA_RANGE=var_eta_range, XI_RANGE=var_xi_range)

   ;; Establish records to be processed

   has_time = strlen(grid.dims.time) gt 0

   if has_time then begin
      time_var = ohis->TimeVarName(grid.dims.time)
      if isa(time_var, /NULL) then message, 'Time variable not found'
      if ohis->HasAtt(time_var, 'units') then begin
         time_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
      endif else begin
         time_units = {scale: 1, offset: 0}
      endelse
      time = ohis->VarGet(time_var, AUTOSCALE=0)*time_units.scale
      if n_elements(time_range) gt 0 then begin
         record_range = mgh_subset(time, time_range)
      endif
      if n_elements(record_range) eq 0 then begin
         n_time = ohis->DimInfo(grid.dims.time, /DIMSIZE)
         record_range = [0,n_time-1]
      endif
      if n_elements(time_range) eq 0 then $
         time_range = time[record_range]
      fmt = '(%"Getting %s data between records %d %d, times %s %s")'
      message, /INFORM, string(FORMAT=fmt, string(variable, /PRINT), record_range, mgh_format_float(time_range))
      rra0 = record_range[0]
      rra1 = record_range[1]
      rran = rra1-rra0+1
   endif else begin
      rran = 1
      fmt = '(%"Variable %s does not vary with time")'
      message, /INFORM, string(FORMAT=fmt, string(variable, /PRINT))
   endelse

   ;; Calculate values to be plotted if necessary

   if keyword_set(recalc) then begin

      message, /INFORM, 'Calculating '+parameter+' data'

      ;; Get data. The DATA_MULTIPLIER and DATA_TRANSFORMATION arguments are applied *before*
      ;; statistical processing

      dim = [grid.xi_range[1]-grid.xi_range[0]+1,grid.eta_range[1]-grid.eta_range[0]+1]

      ;; For some parameters, we need to get all the data before
      ;; calculating the result

      get_all = ~ (strmatch(parameter, 'mean', /FOLD_CASE) || $
         strmatch(parameter, 'max', /FOLD_CASE) || strmatch(parameter, 'min', /FOLD_CASE) || $
         strmatch(parameter, 'fraction*', /FOLD_CASE))

      if get_all then begin
         if has_time then begin
            data = fltarr([dim,rran])
            for r=0,rran-1 do begin
               data[0,0,r] = ohis->HsliceData(variable, GRID=grid, RECORD=rra0+r, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
            endfor
         endif else begin
            data = ohis->HsliceData(variable, GRID=grid, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
         endelse
         if n_elements(data_multiplier) gt 0 then $
            data *= data_multiplier
         if n_elements(data_transformation) gt 0 then $
            data = call_function(data_transformation, data)
      endif

      case !true of

         strmatch(parameter, 'mean', /FOLD_CASE): begin
            if has_time then begin
               values = dblarr(dim)
               for r=0,rran-1 do begin
                  data = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                     RECORD=rra0+r, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
                  if n_elements(data_multiplier) gt 0 then $
                     data *= data_multiplier
                  if n_elements(data_transformation) gt 0 then $
                     data = call_function(data_transformation, data)
                  values += temporary(data)
               endfor
               values /= rran
            endif else begin
               data = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                  DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
               if n_elements(data_multiplier) gt 0 then $
                  data *= data_multiplier
               if n_elements(data_transformation) gt 0 then $
                  data = call_function(data_transformation, data)
               values = temporary(data)
            endelse
         end

         strmatch(parameter, 'min', /FOLD_CASE): begin
            if has_time then begin
               for r=0,rran-1 do begin
                  data = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                     RECORD=rra0+r, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
                  if n_elements(data_multiplier) gt 0 then $
                     data *= data_multiplier
                  if n_elements(data_transformation) gt 0 then $
                     data = call_function(data_transformation, data)
                  values = r gt 0 ? values < temporary(data) : temporary(data)
               endfor
            endif else begin
               data = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                  DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
               if n_elements(data_multiplier) gt 0 then $
                  data *= data_multiplier
               if n_elements(data_transformation) gt 0 then $
                  data = call_function(data_transformation, data)
               values = temporary(data)
            endelse
         end

         strmatch(parameter, 'max', /FOLD_CASE): begin
            if has_time then begin
               for r=0,rran-1 do begin
                  data = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                     RECORD=rra0+r, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
                  if n_elements(data_multiplier) gt 0 then $
                     data *= data_multiplier
                  if n_elements(data_transformation) gt 0 then $
                     data = call_function(data_transformation, data)
                  values = r gt 0 ? values > temporary(data) : temporary(data)
               endfor
            endif else begin
               data = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                  DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
               if n_elements(data_multiplier) gt 0 then $
                  data *= data_multiplier
               if n_elements(data_transformation) gt 0 then $
                  data = call_function(data_transformation, data)
               values = temporary(data)
            endelse
         end

         strmatch(parameter, 'fraction > ?*', /FOLD_CASE): begin
            pp = strsplit(parameter, /EXTRACT)
            threshold = float(pp[2])
            if has_time then begin
               values = dblarr(dim)
               for r=0,rran-1 do begin
                  data = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                     RECORD=rra0+r, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
                  if n_elements(data_multiplier) gt 0 then $
                     data *= data_multiplier
                  if n_elements(data_transformation) gt 0 then $
                     data = call_function(data_transformation, data)
                  values += temporary(data) gt threshold
               endfor
               values /= rran
            endif else begin
               data = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                  DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
               if n_elements(data_multiplier) gt 0 then $
                  data *= data_multiplier
               if n_elements(data_transformation) gt 0 then $
                  data = call_function(data_transformation, data)
               values = temporary(data) gt threshold
            endelse
         end

         strmatch(parameter, 'median', /FOLD_CASE): begin
            dim = size(data, /DIMENSIONS)
            values = make_array(dim[0:1], VALUE=!values.f_nan)
            for j=0,dim[1]-1 do begin
               for i=0,dim[0]-1 do begin
                  d = data[i,j,*]
                  l_good = where(finite(d), n_good)
                  if n_good gt 0 then $
                     values[i,j] = median(d[l_good])
               endfor
            endfor
         end

         strmatch(parameter, 'median (DJF)', /FOLD_CASE): begin
            dt_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
            ts = mgh_dt_caldat(time[rra0:rra1]+dt_units.offset)
            dim = size(data, /DIMENSIONS)
            values = make_array(dim[0:1], VALUE=!values.f_nan)
            for j=0,dim[1]-1 do begin
               for i=0,dim[0]-1 do begin
                  d = data[i,j,*]
                  l_good = where(finite(d) and (ts.month eq 12 or ts.month eq 1 or ts.month eq 2), n_good)
                  if n_good gt 0 then $
                     values[i,j] = median(d[l_good])
               endfor
            endfor
         end

         strmatch(parameter, 'median (MAM)', /FOLD_CASE): begin
            dt_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
            ts = mgh_dt_caldat(time[rra0:rra1]+dt_units.offset)
            dim = size(data, /DIMENSIONS)
            values = make_array(dim[0:1], VALUE=!values.f_nan)
            for j=0,dim[1]-1 do begin
               for i=0,dim[0]-1 do begin
                  d = data[i,j,*]
                  l_good = where(finite(d) and (ts.month eq 3 or ts.month eq 4 or ts.month eq 5), n_good)
                  if n_good gt 0 then $
                     values[i,j] = median(d[l_good])
               endfor
            endfor
         end

         strmatch(parameter, 'median (JJA)', /FOLD_CASE): begin
            dt_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
            ts = mgh_dt_caldat(time[rra0:rra1]+dt_units.offset)
            dim = size(data, /DIMENSIONS)
            values = make_array(dim[0:1], VALUE=!values.f_nan)
            for j=0,dim[1]-1 do begin
               for i=0,dim[0]-1 do begin
                  d = data[i,j,*]
                  l_good = where(finite(d) and (ts.month eq 6 or ts.month eq 7 or ts.month eq 8), n_good)
                  if n_good gt 0 then $
                     values[i,j] = median(d[l_good])
               endfor
            endfor
         end

         strmatch(parameter, 'median (SON)', /FOLD_CASE): begin
            dt_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
            ts = mgh_dt_caldat(time[rra0:rra1]+dt_units.offset)
            dim = size(data, /DIMENSIONS)
            values = make_array(dim[0:1], VALUE=!values.f_nan)
            for j=0,dim[1]-1 do begin
               for i=0,dim[0]-1 do begin
                  d = data[i,j,*]
                  l_good = where(finite(d) and (ts.month eq 9 or ts.month eq 10 or ts.month eq 11), n_good)
                  if n_good gt 0 then $
                     values[i,j] = median(d[l_good])
               endfor
            endfor
         end

         strmatch(parameter, 'standard deviation', /FOLD_CASE): begin
            values = mgh_stdev(data, 3, /NAN)
         end

         strmatch(parameter, 'percentile*', /FOLD_CASE): begin
            p = strsplit(parameter, /EXTRACT)
            threshold = n_elements(p) gt 1 ? float(p[1]) : 50
            message, /INFORM, string(FORMaT='(%"Calculating %sth percentile")', mgh_format_float(threshold))
            values = cmapply('user:mgh_percentile', data, 3, FUNCTARGS={threshold: threshold, nan: 1B})
         end

         strmatch(parameter, 'annual amplitude', /FOLD_CASE): begin
            dim = size(data, /DIMENSIONS)
            tim = time[rra0:rra1]
            annual_period = ohis->HasAtt(time_var, 'cycle_length') ? ohis->AttGet(time_var, 'cycle_length') : 365.
            n_harm = 1
            coeff = replicate(!values.f_nan, [dim[0:1],2*n_harm+1])
            for j=0,dim[1]-1 do begin
               for i=0,dim[0]-1 do begin
                  d = reform(data[i,j,*])
                  if min(finite(d)) eq 0 then continue
                  a = fltarr(2*n_harm+1)
                  !null = mgh_curvefit(tim, temporary(d), replicate(1., dim[2]), a, FUNCTION_NAME='mgh_fit_annual', CYCLE_LENGTH=annual_period, N_HARMONICS=n_harm)
                  coeff[i,j,*] = temporary(a)
               endfor
            endfor
            values = sqrt(coeff[*,*,1]^2+coeff[*,*,2]^2)
         end

         strmatch(parameter, 'annual phase', /FOLD_CASE): begin
            dim = size(data, /DIMENSIONS)
            tim = time[rra0:rra1]
            annual_period = ohis->HasAtt(time_var, 'cycle_length') ? ohis->AttGet(time_var, 'cycle_length') : 365.
            n_harm = 1
            coeff = replicate(!values.f_nan, [dim[0:1],2*n_harm+1])
            for j=0,dim[1]-1 do begin
               for i=0,dim[0]-1 do begin
                  d = reform(data[i,j,*])
                  if min(finite(d)) eq 0 then continue
                  a = fltarr(2*n_harm+1)
                  !null = mgh_curvefit(tim, temporary(d), replicate(1., dim[2]), a, FUNCTION_NAME='mgh_fit_annual', CYCLE_LENGTH=annual_period, N_HARMONICS=n_harm)
                  coeff[i,j,*] = temporary(a)
               endfor
            endfor
            values = atan(coeff[*,*,2], coeff[*,*,1])/(2*!pi)
            values = (values + 1) mod 1
         end

         strmatch(parameter, 'annual phase (months)', /FOLD_CASE): begin
            dim = size(data, /DIMENSIONS)
            tim = time[rra0:rra1]
            annual_period = ohis->HasAtt(time_var, 'cycle_length') ? ohis->AttGet(time_var, 'cycle_length') : 365.
            n_harm = 1
            coeff = replicate(!values.f_nan, [dim[0:1],2*n_harm+1])
            for j=0,dim[1]-1 do begin
               for i=0,dim[0]-1 do begin
                  d = reform(data[i,j,*])
                  if min(finite(d)) eq 0 then continue
                  a = fltarr(2*n_harm+1)
                  !null = mgh_curvefit(tim, temporary(d), replicate(1., dim[2]), a, FUNCTION_NAME='mgh_fit_annual', CYCLE_LENGTH=annual_period, N_HARMONICS=n_harm)
                  coeff[i,j,*] = temporary(a)
               endfor
            endfor
            values = atan(coeff[*,*,2], coeff[*,*,1])/(2*!pi)
            values = (values + 1) mod 1
            values = 12*values
         end

         strmatch(parameter, 'semi-annual amplitude', /FOLD_CASE): begin
            dim = size(data, /DIMENSIONS)
            tim = time[rra0:rra1]
            annual_period = ohis->HasAtt(time_var, 'cycle_length') ? ohis->AttGet(time_var, 'cycle_length') : 365.
            n_harm = 2
            coeff = replicate(!values.f_nan, [dim[0:1],2*n_harm+1])
            for j=0,dim[1]-1 do begin
               for i=0,dim[0]-1 do begin
                  d = reform(data[i,j,*])
                  if min(finite(d)) eq 0 then continue
                  a = fltarr(2*n_harm+1)
                  !null = mgh_curvefit(tim, temporary(d), replicate(1., dim[2]), a, FUNCTION_NAME='mgh_fit_annual', CYCLE_LENGTH=annual_period, N_HARMONICS=n_harm)
                  coeff[i,j,*] = temporary(a)
               endfor
            endfor
            values = sqrt(coeff[*,*,3]^2+coeff[*,*,4]^2)
         end

         strmatch(parameter, 'semi-annual phase', /FOLD_CASE): begin
            dim = size(data, /DIMENSIONS)
            tim = time[rra0:rra1]
            annual_period = ohis->HasAtt(time_var, 'cycle_length') ? ohis->AttGet(time_var, 'cycle_length') : 365.
            n_harm = 2
            coeff = replicate(!values.f_nan, [dim[0:1],2*n_harm+1])
            for j=0,dim[1]-1 do begin
               for i=0,dim[0]-1 do begin
                  d = reform(data[i,j,*])
                  if min(finite(d)) eq 0 then continue
                  a = fltarr(2*n_harm+1)
                  !null = mgh_curvefit(tim, temporary(d), replicate(1., dim[2]), a, FUNCTION_NAME='mgh_fit_annual', CYCLE_LENGTH=annual_period, N_HARMONICS=n_harm)
                  coeff[i,j,*] = temporary(a)
               endfor
            endfor
            values = atan(coeff[*,*,4], coeff[*,*,3])/(2*!pi)
            l_neg = where(values lt 0, n_neg)
            if n_neg gt 0 then values[l_neg] += 1.0
         end

      endcase

      if has_time then begin
         result = {grid: grid, values: values, time_range: time_range, time_offset: time_units.offset}
      endif else begin
         result = {grid: grid, values: values}

      endelse

      message, /INFORM, 'Saving '+parameter+' data to '+tmpfile
      save, FILE=tmpfile, result

   endif else begin

      message, /INFORM, 'Restoring '+parameter+' data from '+tmpfile
      restore, FILE=tmpfile

   endelse

   return, result

end
