;+
; NAME:
;   MGH_ROMS_BAROTROPIC
;
; PURPOSE:
;   Calculated barotropic stream function from time-averaged data from
;   a ROMS simulation.
;
;###########################################################################
; Copyright (c) 2003-2012 NIWA:
;   http://www.niwa.co.nz/
; Licensed under the MIT open source license:
;   http://www.opensource.org/licenses/mit-license.php
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-09:
;     Written.
;   Mark Hadfield, 2004-02:
;     Stripped out code to handle multiple simulations. This is better
;     handled at a higher level.
;   Mark Hadfield, 2015-07:
;     Added temporary file capability.
;   Mark Hadfield, 2016-01:
;     - Removed default time_range.
;     - Now accepts RECORD_RANGE or TIME_RANGE keywords.
;     - The names of the ubar, vbar and zeta variables can now be specified,
;       to allow statistics like 'ubar_mean', 'vbar_mean' and 'zeta_mean'
;     - Added the ability to process subsets of the domain with
;       XI_RANGE and ETA_RANGE keywords.
;-
function mgh_roms_barotropic, ohis, $
     LOC0=loc0, RECALC=recalc, $
     RECORD_RANGE=record_range, TIME_RANGE=time_range, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
     VAR_UBAR=var_ubar, VAR_VBAR=var_vbar, VAR_ZETA=var_zeta, $
     USE_ZETA=use_zeta

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(ohis) eq 0 then message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'ohis'

   if n_elements(var_ubar) eq 0 then var_ubar = 'ubar'
   if n_elements(var_vbar) eq 0 then var_vbar = 'vbar'
   if n_elements(var_zeta) eq 0 then var_zeta = 'zeta'

   if n_elements(use_zeta) eq 0 then use_zeta = !true

   ohis->GetProperty, FILE_NAME=fname

   str_file = mgh_format_integer(mgh_hashcode(fname))
   str_var = mgh_format_integer(mgh_hashcode([var_ubar,var_vbar,var_zeta]))
   str_use = mgh_format_integer(keyword_set(use_zeta))

   name = ['mgh_roms_barotropic',str_file,'var',str_var,'use',str_use]
   if n_elements(time_range) gt 0 then $
      name = [name,string(FORMaT='(%"tr_%0.3f_%0.3f")', time_range)]
   if n_elements(loc0) gt 0 then $
      name = [name,string(FORMaT='(%"loc0_%0.6f_%0.6f")', loc0)]
   if n_elements(xi_range) gt 0 then $
      name = [name,string(FORMaT='(%"xr_%d_%d")', xi_range)]
   if n_elements(eta_range) gt 0 then $
      name = [name,string(FORMaT='(%"er_%d_%d")', eta_range)]
   tmpfile = filepath(strjoin(name, '_')+'.idl_data', /TMP)

   if ~ file_test(tmpfile) then recalc = !true

  if ~ keyword_set(recalc) then begin
     message, /INFORM, 'Retrieving data from temporary file: '+tmpfile
     restore, FILE=tmpfile
     return, result
  endif

   ;; Get grid info

   dim = ohis->DimRho()

   if n_elements(xi_range) eq 0 then xi_range = [0,dim[0]-1]
   if n_elements(eta_range) eq 0 then eta_range = [0,dim[1]-1]

   xra0 = xi_range[0]
   xra1 = xi_range[1]
   xran = xra1 - xra0 + 1

   era0 = eta_range[0]
   era1 = eta_range[1]
   eran = era1 - era0 + 1

   offset = [xra0,era0]
   count = [xran,eran]

   pm = ohis->VarGet('pm', OFFSET=offset, COUNT=count)
   pn = ohis->VarGet('pn', OFFSET=offset, COUNT=count)

   if ohis->HasVar('mask_rho') then begin
      mask_rho = ohis->VarGet('mask_rho', OFFSET=offset, COUNT=count)
      mask_psi = mask_rho[0:-2,0:-2] * $
                 mask_rho[0:-2,1:-1] * $
                 mask_rho[1:-1,0:-2] * $
                 mask_rho[1:-1,1:-1]
   end

   lon = ohis->VarGet('lon_rho', OFFSET=offset, COUNT=count)
   lat = ohis->VarGet('lat_rho', OFFSET=offset, COUNT=count)

   lon_psi = mgh_stagger(lon, DELTA=[-1,-1])
   lat_psi = mgh_stagger(lat, DELTA=[-1,-1])

   h = ohis->VarGet('h', OFFSET=offset, COUNT=count)

   h_u = mgh_stagger(h, DELTA=[-1,0])
   h_v = mgh_stagger(h, DELTA=[0,-1])
   mgh_undefine, h

   dims_ubar = ohis->VarDims(var_ubar)
   dims_vbar = ohis->VarDims(var_vbar)

   ;; Establish records to be processed

   has_time = strlen(dims_ubar.time) gt 0

   if has_time then begin
      time_var = ohis->TimeVarName(dims_ubar.time)
      if isa(time_var, /NULL) then message, 'Time variable not found'
      time = ohis->VarGet(time_var, /AUTOSCALE)
      if n_elements(time_range) gt 0 then begin
         record_range = mgh_subset(time, time_range)
      endif
      if n_elements(record_range) eq 0 then begin
         n_time = ohis->DimInfo(dims_ubar.time, /DIMSIZE)
         record_range = [0,n_time-1]
      endif
      if n_elements(time_range) eq 0 then $
         time_range = time[record_range]
      msg = ['Getting barotropic velocity data between records', $
         strtrim(record_range,2),'times', $
         mgh_format_float(time_range)]
      message, /INFORM, strjoin(temporary(msg), ' ')
      rra0 = record_range[0]
      rra1 = record_range[1]
      rran = rra1-rra0+1
   endif else begin
      rran = 1
      msg = ['Variable', string(*self.variable, /PRINT), 'does not vary with time']
      message, /INFORM, strjoin(temporary(msg), ' ')
   endelse

   ubar = 0
   vbar = 0
   zeta = 0

   for r=0,rran-1 do begin

      offset_ubar = offset
      count_ubar = count-[1,0]
      offset_vbar = offset
      count_vbar = count-[0,1]
      offset_zeta = offset
      count_zeta = count
      if has_time then begin
         offset_ubar = [offset_ubar,rra0+r]
         count_ubar = [count_ubar,1]
         offset_vbar = [offset_vbar,rra0+r]
         count_vbar = [count_vbar,1]
         offset_zeta = [offset_zeta,rra0+r]
         count_zeta = [count_zeta,1]
      endif

      ubar += ohis->VarGet(var_ubar, OFFSET=offset_ubar, COUNT=count_ubar, /AUTOSCALE)/float(rran)
      vbar += ohis->VarGet(var_vbar, OFFSET=offset_vbar, COUNT=count_vbar, /AUTOSCALE)/float(rran)

      if keyword_set(use_zeta) then begin
         zeta += ohis->VarGet(var_zeta, OFFSET=offset_zeta, COUNT=count_zeta, /AUTOSCALE)/float(rran)
      endif

   endfor

   ubar[where(~ finite(ubar), /NULL)] = 0
   vbar[where(~ finite(vbar), /NULL)] = 0

   if keyword_set(use_zeta) then begin
      zeta = mgh_fill2d(zeta)
      zeta_u = mgh_stagger(zeta, DELTA=[-1,0])
      zeta_v = mgh_stagger(zeta, DELTA=[0,-1])
   endif else begin
      zeta_u = 0
      zeta_v = 0
   endelse

   psi = mgh_roms_psi(ubar*(h_u+zeta_u), vbar*(h_v+zeta_v), pm, pn, MASK=mask_rho)

   vort = mgh_roms_curl(ubar, vbar, pm, pn)

   if n_elements(loc0) gt 0 then begin
      loc = round(mgh_locate2a(lon, lat, XOUT=loc0[0], YOUT=loc0[1]))
      psi -= psi[loc[0],loc[1]]
   endif

   if n_elements(mask_psi) gt 0 then begin
      result = {psi: psi, vort: vort, lon: lon_psi, lat: lat_psi, mask: mask_psi}
   endif else begin
      result = {psi: psi, vort: vort, lon: lon_psi, lat: lat_psi}
   endelse

  message, /INFORM, 'Saving data to temporary file: '+tmpfile
  save, FILE=tmpfile, result

  return, result

end
