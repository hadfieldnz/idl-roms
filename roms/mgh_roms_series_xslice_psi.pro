;+
; FUNCTION NAME:
;   MGH_ROMS_SERIES_XSLICE_PSI
;
; PURPOSE:
;   This function extracts and returns a time series of barotropic stream function
;   along a specified section (Xslice) using data from a ROMS history file.
;
;   See also: MGH_ROMS_SERIES_XSLICE_TPORT, MGH_ROMS_SERIES_PSLICE_TPORT
;
; CALLING SEQUENCE
;   result = mgh_roms_series_xslice_psi(ohis, grid)
;
; RETURN VALUE:
;   The function returns the time series and ancillary information in
;   a structure.
;
; POSITIONAL PARAMETERS:
;   ohis (input, object)
;     A reference to an MGHromsHistory object.
;
;   grid (input, structure)
;     The Xslice grid structure.
;
; KEYWORD ARGUMENTS:
;
;   VAR_UBAR (input, string scalar)
;   VAR_VBAR (input, string scalar)
;   VAR_ZETA (input, string scalar)
;     Names for the variables representing depth-average velocity and sea surface
;     height. Defaults are "ubar", "vbar" and "zeta".
;
;   USE_ZETA (input, switch)
;     This keyword controls whether zeta information from the history file is used
;     to calculate the transport. The default is !true (use zeta).
;
;   RECORD_RANGE (integer 2-element vector)
;   TIME_RANGE (numeric 2-element vector)
;     Specify the range of records to be extracted either as either indices
;     (RECORD_RANGE) or times (in days) from the simulation's reference time
;     (TIME_RANGE).
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2016-09:
;     Written.
;   Mark Hadfield, 2019-11:
;     Fixed bug: the function destroyed its ohis object.
;   Mark Hadfield, 2019-11:
;     The function is now supplied with a grid structure via its "grid"
;     argument rather than constructing one internally.
;-
function mgh_roms_series_xslice_psi, ohis, grid, RECALC=recalc, $
     RECORD_RANGE=record_range, TIME_RANGE=time_range, $
     VAR_UBAR=var_ubar, VAR_vbar=var_vbar, VAR_ZETA=var_zeta, USE_ZETA=use_zeta, $
     WHOLE_DOMAIN=whole_domain

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process file argument

   if n_elements(ohis) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'ohis'
   if ~ obj_valid(ohis) gt 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_objref_bad', 'ohis'
   if n_elements(grid) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'grid'
   if ~ isa(grid, 'STRUCT') then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'grid'

   ;; Process arguments controlling variable selection

   if n_elements(var_ubar) eq 0 then var_ubar = 'ubar'
   if n_elements(var_vbar) eq 0 then var_vbar = 'vbar'
   if n_elements(var_zeta) eq 0 then var_zeta = 'zeta'

   if n_elements(use_zeta) eq 0 then use_zeta = !true

   ;; Get & check variable dimensions.

   dim_ubar = ohis->VarDims(var_ubar)
   dim_vbar = ohis->VarDims(var_vbar)

   if min([strlen(dim_ubar.horizontal),strlen(dim_ubar.time)]) eq 0 then $
      message, 'The Ubar variable must have horizontal & time dimensions'
   if min([strlen(dim_vbar.horizontal),strlen(dim_vbar.time)]) eq 0 then $
      message, 'The Vbar variable must have horizontal & time dimensions'
   if dim_ubar.time ne dim_vbar.time then $
      message, 'The Ubar & Vbar variables have different time dimensions'
   if ~ array_equal(dim_ubar.horizontal, ['xi_u','eta_u']) then $
      message, 'The U variable has incorrect horizontal dimensions'
   if ~ array_equal(dim_vbar.horizontal, ['xi_v','eta_v']) then $
      message, 'The V variable has incorrect horizontal dimensions'

   ;; Read grid dimensions

   dim_rho = ohis->DimRho()

   ;; Construct temporary file name & return data from the temporary file
   ;; if recalculation is not required

   ohis->GetProperty, FILE_NAME=fname

   str_file = mgh_format_integer(mgh_hashcode(fname))
   str_var = mgh_format_integer(mgh_hashcode([var_ubar,var_vbar,var_zeta]))
   str_use = mgh_format_integer(keyword_set(use_zeta))

   tmpfile = list('mgh_roms_series_pslice_psi', 'file', str_file, 'var', str_var, $
      'use', str_use, 'grid', mgh_format_integer(mgh_hashcode(grid)))
   if n_elements(time_range) gt 0 then $
      tmpfile->Add, ['tr',mgh_format_float(time_range)], /EXTRACT
   tmpfile = filepath(strjoin(tmpfile->ToArray(), '_')+'.idl_data', /TMP)

   if ~ file_test(tmpfile) then recalc = !true

   if ~ keyword_set(recalc) then begin
      message, /INFORM, 'Restoring ROMS Xslice stream function time series data from '+tmpfile
      restore, FILE=tmpfile
      return, result
   endif

   ;; Establish time variable name and records to be extracted. Get time data

   n_time = ohis->DimInfo(dim_ubar.time, /DIMSIZE)
   mgh_resolve_indices, n_time, record_range

   case !true of
      ohis->HasVar(dim_ubar.time): $
         time_var = dim_ubar.time
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
      n_time = ohis->DimInfo(dim_ubar.time, /DIMSIZE)
      record_range = [0,n_time-1]
   endif

   rr0 = record_range[0]
   rr1 = record_range[1]
   rrn = rr1 - rr0 + 1

   time = time[rr0:rr1]

   ;; Process xi & eta ranges

   if keyword_set(whole_domain) then begin
      ;; With x_range and eta_range unset, the mgh_roms_barotropic function
      ;; below will process the whole domain.
   endif else begin
      xi_range = [floor(min(grid.xi)-1) > 0, ceil(max(grid.xi)+1) < (dim_rho[0]-1)]
      eta_range = [floor(min(grid.eta)-1) > 0, ceil(max(grid.eta)+1) > dim_rho[1]-1]
   endelse

   ;; Calculate transport

   tport = reform(dblarr(grid.n_points*rrn), [grid.n_points,rrn])

   for r=rr0,rr1 do begin
      b = mgh_roms_barotropic(ohis, RECALC=recalc, RECORD_RANGE=[r,r], XI_RANGE=xi_range, ETA_RANGE=eta_range, $
         VAR_UBAR=var_ubar, VAR_VBAR=var_vbar, VAR_ZETA=var_zeta, UsE_ZETA=use_zeta)
      if r eq rr0 then begin
         loc = mgh_locate2(b.lon, b.lat, XOUT=grid.x, YOUT=grid.y)
         loc_xi = reform(loc[0,*])
         loc_eta = reform(loc[1,*])
      endif
      p = mgh_interpolate(b.psi, loc_xi, loc_eta)
      tport[*,r-rr0] = p - p[0]
      mgh_undefine, b, p
   endfor

   result = {time: time, tport: tport}

   message, /INFORM, 'Saving ROMS Xslice stream function time series data to '+tmpfile
   save, FILE=tmpfile, result

   return, result

end
