; svn $Id$
;+
; NAME:
;   MGH_ROMS_GPAN
;
; PURPOSE:
;   Calculate geopotential anomaly from fields in a ROMS history file
;
; POSITIONAL PARAMETERS:
;   ohis
; KEYWORD PARAMETERS:
;   DEPTH (input, numeric scalar)
;     Depth to which geopotential anomaly is referenced.
;
;   FILL (input, numeric scalar)
;     Depth to which geopotential anomaly is referenced.
;
; RETURN VALUE:
;   The function returns a 2D array containing surface geopotential
;   anomaly referenced to DEPTH.
;
; PROCEDURE:
;   Interpolate to horizontal planes, optionally fill, then calculate
;
;###########################################################################
;
; This software is provided subject to the following conditions:
;
; 1.  NIWA makes no representations or warranties regarding the
;     accuracy of the software, the use to which the software may
;     be put or the results to be obtained from the use of the
;     software.  Accordingly NIWA accepts no liability for any loss
;     or damage (whether direct of indirect) incurred by any person
;     through the use of or reliance on the software.
;
; 2.  NIWA is to be acknowledged as the original author of the
;     software where the software is used or presented in any form.
;
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2004-02:
;     Written.
;   Mark Hadfield, 2009-04:
;     Updated for new ROMS vertical transform and vertical stretching
;     options, but introduced a bug: VTRANSFORM not passed to
;     MGH_ROMS_S_TO_Z.
;     from file!
;   Mark Hadfield, 2009-10:
;     - Removed calls to widget_event(/NOWAIT).
;     - Fixed bug introduced in 2009-04.
;-
function mgh_roms_gpan, ohis, $
     DEPTH=depth, FILL=fill, RECALC=recalc, TIME_RANGE=time_range

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(ohis) eq 0 then message, 'No ROMS file'

   if n_elements(depth) eq 0 then depth = 2000

   if n_elements(fill) eq 0 then fill = 1B

   if n_elements(time_range) eq 0 then time_range = [20,50]

   ;; Construct name for temporary file

   ohis->GetProperty, NAME=name

   tmpfile = ['mgh_roms_gpan', $
              !version.release, $
              name, $
              string(depth, FORMAT='(F0.0)'), $
              string(keyword_set(fill), FORMAT='(I0)'), $
              string(time_range, FORMAT='(F0.0)')]

   tmpfile = filepath(strjoin(tmpfile, '_')+'.idl_data', /TMP)

   ;; If applicable, get data from temporary file and return it

   if ~ file_test(tmpfile, /READ) then recalc = 1B

   if ~ keyword_set(recalc) then begin
      fmt = '(%"Restoring geopotential anomaly data from temporary file %s")'
      message, /INFORM, string(FORMAT=temporary(fmt), tmpfile)
      restore, FILENAME=tmpfile, /RELAX
      return, result
   endif

   ;; Recalculate

   fmt = '(%"ROMS geopotential anomaly data for %s will be calculated & ' + $
         'saved to temporary file %s")'
   message, /INFORM, string(FORMAT=temporary(fmt), name, tmpfile)

   ;; Get grid info from first history file

   dim = ohis->DimRho()

   lon = ohis->VarGet('lon_rho')
   lat = ohis->VarGet('lat_rho')

   h = ohis->VarGet('h')

   if ohis->HasVar('mask_rho') then mask = ohis->VarGet('mask_rho')

   s = ohis->HasVar('s_rho') ? ohis->VarGet('s_rho') : ohis->VarGet('sc_r')

   theta_s = ohis->VarGet('theta_s')
   theta_b = ohis->VarGet('theta_b')

   hc = ohis->VarGet('hc')

   vstretch = ohis->HasVar('Vstretching') ? ohis->VarGet('Vstretching') : 1
   vtransform = ohis->HasVar('Vtransform') ? ohis->VarGet('Vtransform') : 1

   cs = mgh_roms_s_to_cs(s, $
                         THETA_S=theta_s, THETA_B=theta_b, $
                         VSTRETCH=vstretch)

   ;; Set up z-grid depths

   if n_elements(depth) eq 0 then depth = 2000

   n_z = round(depth/25.)

   z = mgh_range(0, depth, N_ELEMENTS=n_z)

   ;; Retrieve time-averaged temperature and salinity fields. Don't
   ;; retrieve all data, then average over time, because this would
   ;; use way too much memory.

   temp = fltarr(dim)
   salt = fltarr(dim)

   tdims = ohis->VarDims('temp')
   tday = ohis->VarGet(tdims.time)

   time_pad = 1.E-3

   rrange = mgh_minmax(where(tday gt time_range[0]-time_pad and $
                             tday lt time_range[1]+time_pad))

   msg = ['Getting data between records', strtrim(rrange,2), $
          'times', mgh_format_float(tday[rrange])]

   message, /INFORM, strjoin(temporary(msg), ' ')

   rra0 = rrange[0]  &  rra1 = rrange[1]  &  rran = rra1-rra0+1

   for r=rra0,rra1 do begin
      temp += ohis->VarGet('temp', OFFSET=[0,0,0,r], COUNT=[0,0,0,1])/float(rran)
      salt += ohis->VarGet('salt', OFFSET=[0,0,0,r], COUNT=[0,0,0,1])/float(rran)
   endfor

   ;; Interpolate temperature & salinity to z-grid

   temp_z = replicate(!values.f_nan, [dim[0:1],n_z])
   salt_z = replicate(!values.f_nan, [dim[0:1],n_z])

   for i=0,dim[0]-1 do begin
      for j=0,dim[1]-1 do begin
         if n_elements(mask) eq 0 || mask[i,j] gt 0 then begin
            zz = mgh_roms_s_to_z(s, h[i,j], HC=hc, CS=cs, $
                                 VTRANSFORM=vtransform)
            tzz = interpol(temp[i,j,*], zz, -z)
            szz = interpol(salt[i,j,*], zz, -z)
            bot = where(z gt h[i,j], n_bot)
            if n_bot gt 0 then begin
               tzz[bot] = !values.f_nan
               szz[bot] = !values.f_nan
            endif
            temp_z[i,j,*] = tzz
            salt_z[i,j,*] = szz
         endif
      endfor
   endfor

   ;; Optionally fill data

   if keyword_set(fill) then begin
      for k=0,n_z-1 do begin
         temp_z[*,*,k] = mgh_fill2d(temp_z[*,*,k], METHOD=1)
         salt_z[*,*,k] = mgh_fill2d(salt_z[*,*,k], METHOD=1)
      endfor
   endif

   ;; Calculate GPAN. Ignore latitudinal dependence of gravity, as
   ;; this is ignored in the model. Convert model temperature (actually
   ;; potential temperature) to true temperature

   case n_elements(mask) gt 0 of
      0B: result = {gpan: replicate(!values.f_nan, dim[0:1]), lon: lon, lat: lat}
      1B: result = {gpan: replicate(!values.f_nan, dim[0:1]), lon: lon, lat: lat, $
                    mask: mask}
   endcase

   p = mgh_sw_pres(z, /PASCAL)

   for i=0,dim[0]-1 do begin
      for j=0,dim[1]-1 do begin
         t = reform(temp_z[i,j,*])
         s = reform(salt_z[i,j,*])
         if min(finite(t)) gt 0 && min(finite(s)) gt 0 then begin
            t = mgh_sw_temp(t, s, p, /PASCAL)
            result.gpan[i,j] = - mgh_sw_gpan(t, s, p, /PASCAL)
         endif
      endfor
   endfor

   save, result, FILENAME=tmpfile

   return, result

end
