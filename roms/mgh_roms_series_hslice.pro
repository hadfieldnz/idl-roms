;+
; FUNCTION NAME:
;   MGH_ROMS_SERIES_HSLICE
;
; PURPOSE:
;   This function extracts and returns a time series of hslice data from
;   a ROMS history file
;
; CALLING SEQUENCE
;   result = mgh_roms_series_hslice(file)
;
; RETURN VALUE:
;   The function returns the time series and ancillary information in
;   a structure.
;
; POSITIONAL PARAMETERS:
;   ofile (input, object)
;     A reference to an MGHromsHistory or MGHromsStation object.
;
; KEYWORD ARGUMENTS:
;   DEPTH (input, numeric scalar or vector)
;     This keyword specifies the depth(s) at which data are to be
;     extracted. This keyword may be specified only for variables
;     having a depth coordinate and it cannot be used together with
;     LEVEL or SIGMA.
;
;   LAYER (input, integer scalar or vector)
;     This keyword specifies the bed layer(s) at which data are to be
;     extracted. This keyword should be specified only for variables
;     having a bed-layer dimension. It is currently supported only
;     for history files.
;
;   LEVEL (input, integer scalar or vector)
;     This keyword specifies the s-coordinate level(s) at which data are
;     to be extracted.  This keyword may be specified only for
;     variables having a depth coordinate and it cannot be used
;     together with DEPTH, PROFILE or SIGMA.

;   RECORD_RANGE (integer 2-element vector)
;   TIME_RANGE (numeric 2-element vector)
;     Specify the range of records to be extracted either as either indices
;     (RECORD_RANGE) or times (in days) from the simulation's reference time
;     (TIME_RANGE).
;
;   SIGMA (integer scalar or vector)
;     This keyword specifies the sigma-coordinate level(s) at which data
;     are to be extracted.  This keyword may be specified only for
;     variables having a depth coordinate and it cannot be used
;     together with DEPTH, LEVEL or PROFILE.
;
;   VARIABLE (input, string or structure scalar)
;     A variable descriptor, must be a string or a structure that can
;     be interpreted by the MGHromsHistory or MGHromsStation object's
;     methods. Default is 'zeta'.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2014-03:
;     Written.
;   Mark Hadfield, 2016-02:
;     Added temporary file support.
;-
function mgh_roms_series_hslice, ohis, $
     DEPTH=depth, LEVEL=level, LAYER=layer, SIGMA=sigma, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
     RECORD_RANGE=record_range, TIME_RANGE=time_range, $
     VARIABLE=variable, RECALC=recalc

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process file argument

   if n_elements(ohis) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'ohis'

   if n_elements(ohis) gt 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'ohis'

   if ~ obj_valid(ohis) gt 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_objref_bad', 'ohis'

   ;; Process VARIABLE argument

   if n_elements(variable) eq 0 then variable = 'zeta'

   ;; Construct a useful temporary file name. Quite an operation!

   ohis->GetProperty, FILE_NAME=fname

   str_file = mgh_format_integer(mgh_hashcode(fname))
   str_var = mgh_format_integer(mgh_hashcode(variable))

   tmpfile = ['mgh_roms_series_hslice',str_file,str_var]
   if n_elements(depth) gt 0 then $
      tmpfile = [tmpfile,'depth',mgh_format_float(depth)]
   if n_elements(level) gt 0 then $
      tmpfile = [tmpfile,'level',mgh_format_integer(level)]
   if n_elements(layer) gt 0 then $
      tmpfile = [tmpfile,'layer',mgh_format_integer(layer)]
   if n_elements(sigma) gt 0 then $
      tmpfile = [tmpfile,'sigma',mgh_format_integer(sigma)]
   if n_elements(xi_range) gt 0 then $
      tmpfile = [tmpfile,'xr',mgh_format_integer(xi_range)]
   if n_elements(eta_range) gt 0 then $
      tmpfile = [tmpfile,'xr',mgh_format_integer(eta_range)]
   if n_elements(record_range) gt 0 then $
      tmpfile = [tmpfile,'rr',mgh_format_integer(record_range)]
   if n_elements(time_range) gt 0 then $
      tmpfile = [tmpfile,'tr',mgh_format_float(time_range)]
   tmpfile = filepath(strjoin(tmpfile, '_')+'.idl_data', /TMP)

   if ~ file_test(tmpfile) then recalc = !true

   if ~ keyword_set(recalc) then begin
      message, /INFORM, 'Restoring ROMS time series data from '+tmpfile
      restore, FILE=tmpfile
      return, result
   endif

   ;; Set ETA_RANGE and XI_RANGE relative to the rho grid.
   ;; The default is to show all interior points.

   dim_rho = [ohis->DimInfo('xi_rho', /DIMSIZE),ohis->DimInfo('eta_rho', /DIMSIZE)]

   if n_elements(xi_range) eq 0 then xi_range = [0,dim_rho[0]-1]
   if n_elements(eta_range) eq 0 then eta_range = [0,dim_rho[1]-1]

   if xi_range[0] lt 0 then xi_range[0] += dim_rho[0]
   if xi_range[1] lt 0 then xi_range[1] += dim_rho[0]

   if eta_range[0] lt 0 then eta_range[0] += dim_rho[1]
   if eta_range[1] lt 0 then eta_range[1] += dim_rho[1]

   ;; Get grid data required for horizontal slice retrievals

   var_xi_range = xi_range
   var_eta_range = eta_range

   var_dims = ohis->VarDims(variable)

   vdh0 = var_dims.horizontal[0]
   if vdh0 eq 'xi_u' || vdh0 eq 'xi_psi' then $
      var_xi_range += [-1,0]

   vdh1 = var_dims.horizontal[1]
   if vdh1 eq 'eta_v' || vdh1 eq 'eta_psi' then $
      var_eta_range += [-1,0]

   grid = ohis->HsliceGrid(variable, XI_RANGE=var_xi_range, ETA_RANGE=var_eta_range)

   ;; Establish records to be processed

   has_time = strlen(grid.dims.time) gt 0

   if has_time then begin
      case 1B of
         ohis->HasVar(grid.dims.time): $
            time_var = grid.dims.time
         ohis->HasVar('ocean_time'): $
            time_var = 'ocean_time'
         ohis->HasVar('scrum_time'): $
            time_var = 'scrum_time'
         else: $
            message, 'Time variable not found'
      endcase
      time = ohis->VarGet(time_var, AUTOSCALE=0)
      time_units = {scale: 1D/(24D*3600D), offset: 0D}
      if ohis->HasAtt(time_var, 'units') then $
         time_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
      time *= time_units.scale
      time_offset = time_units.offset
      mgh_undefine, time_units
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
      message, /INFORM, string(FORMAT=fmt, string(variable, /PRINT), $
         mgh_format_integer(record_range), mgh_format_float(time_range))
      rra0 = record_range[0]
      rra1 = record_range[1]
      rran = rra1-rra0+1
      time = time[rra0:rra1]
   endif else begin
      rran = 1
      msg = ['Variable', string(variable, /PRINT), 'does not vary with time']
      message, /INFORM, strjoin(temporary(msg), ' ')
   endelse

   ;; Set up result structure

   has_vertical = strlen(grid.dims.vertical) gt 0

   dim = [grid.xi_range[1]-grid.xi_range[0]+1,grid.eta_range[1]-grid.eta_range[0]+1]

   result = create_struct('grid', grid)
   if has_time then result = create_struct(result, 'time', time, 'time_offset', time_offset)

   case !true of
      (~ has_vertical): begin
         result = create_struct(result, 'value', fltarr([dim,1,rran]))
      end
      n_elements(level) gt 0: begin
         result = create_struct(result, 'level', level, 'value', fltarr([dim,n_elements(level),rran]))
      end
      n_elements(depth) gt 0: begin
         result = create_struct(result, 'depth', depth, 'value', fltarr([dim,n_elements(depth),rran]))
      end
      n_elements(sigma) gt 0: begin
         result = create_struct(result, 'sigma', sigma, 'value', fltarr([dim,n_elements(sigma),rran]))
      end
      n_elements(layer) gt 0: begin
         result = create_struct(result, 'layer', layer, 'value', fltarr([dim,n_elements(layer),rran]))
      end
      else: begin
         result = create_struct(result, 'value', fltarr([dim,1,rran]))
      end
   endcase

   if has_time then begin
      for r=rra0,rra1 do $
         result.value[*,*,*,r-rra0] = ohis->HSliceData(variable, GRID=grid, DEPTH=depth, LEVEL=level, SIGMA=sigma, LAYER=layer, RECORD=r)
   endif else begin
      result.value = ohis->HSliceData(variable, GRID=grid, DEPTH=depth, LEVEL=level, SIGMA=sigma, LAYER=layer)
   endelse

   message, /INFORM, 'Saving ROMS time series data to '+tmpfile
   save, FILE=tmpfile, result

   return, result

end
