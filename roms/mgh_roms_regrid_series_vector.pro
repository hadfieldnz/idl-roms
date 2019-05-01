;+
; FUNCTION NAME:
;   MGH_ROMS_REGRID_SERIES_VECTOR
;
; PURPOSE:
;   This function extracts and returns a time series of vector data from
;   a ROMS "regrid" forcing file (i.e. one which is regridded by ROMS)
;
; CALLING SEQUENCE
;   result = mgh_roms_regrid_series_vector(file, position)
;
; RETURN VALUE:
;   The function returns the time series and ancillary information in
;   a structure.
;
; POSITIONAL PARAMETERS:
;   ofile (input, object)
;     A reference to an MGHncSequence object.
;
;   location (input, numeric vector)
;     A 2-element numeric vector specifying the longitude and latitude
;     of the point where the data are to be extracted.
;
; KEYWORD ARGUMENTS:
;   NEAREST (input, switch)
;     This keyword controls horizontal interpolation for history file
;     series. The default is 0 (interpolate linearly).
;
;   VARIABLE (input, string or structure vector)
;     A 2-element vector, must be a string or a structure, that can
;     be interpreted by the MGHromsHistory or MGHromsStation object's
;     methods. Default is ['u','v'] if available, otherwise
;     ['ubar','vbar'].
;
;   RECORD_RANGE (integer 2-element vector)
;   TIME_RANGE (numeric 2-element vector)
;     Specify the range of records to be extracted either as either indices
;     (RECORD_RANGE) or times (in days) from the simulation's reference time
;     (TIME_RANGE).
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2009-02:
;     Written.
;   Mark Hadfield, 2009-09:
;     Added SIGMA keyword.
;   Mark Hadfield, 2010-04:
;     Fixed bug: getting whole profile when LEVEL keyword set.
;   Mark Hadfield, 2013-10:
;     - Added HEIGHT (above bottom) keyword.
;     - Modernised indentation.
;   Mark Hadfield, 2014-04:
;     - Ensured that the code allows for vector values of DEPTH, HEIGHT,
;       LEVEL and SIGMA. Cleaned up a few minor bugs in this and updated
;       documentation.
;   Mark Hadfield, 2016-03:
;     - Added temporary-file support.
;     - Added NEAREST keyword.
;-
function mgh_roms_regrid_series_vector, ofile, location, $
     NEAREST=nearest, RECALC=recalc, $
     RECORD_RANGE=record_range, TIME_RANGE=time_range, $
     VARIABLE=variable

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process file argument

   if n_elements(ofile) eq 0 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'ofile'

   if n_elements(ofile) gt 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'ofile'

   if ~ obj_valid(ofile) gt 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_objref_bad', 'ofile'

   ;; Process VARIABLE argument

   if n_elements(variable) eq 0 then begin
      has_uv = ofile->HasVar('u') && ofile->HasVar('v')
      variable = has_uv ? ['u','v'] : ['ubar','vbar']
   endif

   if n_elements(variable) ne 2 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'variable'

   ;; Determine file type

   case !true of
      obj_isa(ofile, 'MGHromsHistory'): ftype = 'history'
      obj_isa(ofile, 'MGHromsStation'): ftype = 'station'
      else: message, 'Unknown file type'
   endcase

   ;; Construct a useful temporary file name. Quite an operation!

   ofile->GetProperty, FILE_NAME=fname

   str_file = mgh_format_integer(mgh_hashcode(fname))
   str_var = mgh_format_integer(mgh_hashcode(variable))

   tmpfile = ['mgh_roms_series_vector',ftype,str_file,str_var]
   if n_elements(position) gt 0 then $
      tmpfile = [tmpfile,'position',mgh_format_float(position)]
   if n_elements(depth) gt 0 then $
      tmpfile = [tmpfile,'depth',mgh_format_float(depth)]
   if n_elements(height) gt 0 then $
      tmpfile = [tmpfile,'height',mgh_format_float(height)]
   if n_elements(level) gt 0 then $
      tmpfile = [tmpfile,'level',mgh_format_integer(level)]
   if n_elements(sigma) gt 0 then $
      tmpfile = [tmpfile,'sigma',mgh_format_integer(sigma)]
   if n_elements(profile) gt 0 then $
      tmpfile = [tmpfile,'profile',mgh_format_integer(profile)]
   if n_elements(nearest) gt 0 then $
      tmpfile = [tmpfile,'nearest',mgh_format_integer(nearest)]
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

   ;; Does the file have (lon,lat) or (x,y) data?

   lonlat = ofile->HasVar('lon_rho') && ofile->HasVar('lat_rho')

   ;; Does the file have dynamic bathymetry? At the moment we read only
   ;; the initial value.

   use_bath = ofile->HasVar('bath')

   ;; Handling of dimensions & variables depends on whether
   ;; the file is a history or stations file

   case ftype of

      'history': begin

         ;; Get & check variable dimensions.

         dim_u = ofile->VarDims(variable[0])
         dim_v = ofile->VarDims(variable[1])

         if min([strlen(dim_u.horizontal),strlen(dim_u.time)]) eq 0 then $
            message, 'The U variable must have horizontal & time dimensions'
         if min([strlen(dim_v.horizontal),strlen(dim_v.time)]) eq 0 then $
            message, 'The V variable must have horizontal & time dimensions'
         if dim_u.time ne dim_v.time then $
            message, 'The U & V variables have different time dimensions'
         if dim_u.vertical ne dim_v.vertical then $
            message, 'The U & V variables have different vertical dimensions'
         if ~ array_equal(dim_u.horizontal, ['xi_u','eta_u']) then $
            message, 'The U variable has incorrect horizontal dimensions'
         if ~ array_equal(dim_v.horizontal, ['xi_v','eta_v']) then $
            message, 'The V variable has incorrect horizontal dimensions'

         ;; Read grid dimensions

         dim_rho = ofile->DimRho()

         ;; Establish position of point to be extracted and calculate
         ;; water depth, location, etc.

         if n_elements(position) eq 0 then $
            position = 0.5*(dim_rho[0:1]-1)

         if n_elements(position) ne 2 then $
            message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'position'

         if keyword_set(nearest) then begin

            ipos = round(position)

            if use_bath then begin
               h = ofile->VarGet('bath', COUNT=[1,1,1], OFFSET=[ipos,0])
            endif else begin
               h = ofile->VarGet('h', COUNT=[1,1], OFFSET=ipos)
            endelse

            if lonlat then begin
               x = ofile->VarGet('lon_rho', COUNT=[1,1], OFFSET=ipos)
               y = ofile->VarGet('lat_rho', COUNT=[1,1], OFFSET=ipos)
            endif else begin
               x = ofile->VarGet('x_rho', COUNT=[1,1], OFFSET=ipos)
               y = ofile->VarGet('y_rho', COUNT=[1,1], OFFSET=ipos)
            endelse

            angle = 0
            if ofile->HasVar('angle') then $
               angle = ofile->VarGet('angle', COUNT=[1,1], OFFSET=ipos)

            fmt = '(%"Nearest point %d, %d at %s, %s, angle %s deg in a water depth of %s m")'
            message, /INFORM, string(FORMAT=fmt, ipos, mgh_format_float([x,y,angle*!radeg,h]))

         endif else begin

            ipos = floor(position)
            fpos = position - ipos

            if use_bath then begin
               h = interpolate(ofile->VarGet('bath', COUNT=[2,2,1], OFFSET=[ipos,0]), fpos[0], fpos[1])
            endif else begin
               h = interpolate(ofile->VarGet('h', COUNT=[2,2], OFFSET=ipos), fpos[0], fpos[1])
            endelse

            if lonlat then begin
               x = interpolate(ofile->VarGet('lon_rho', COUNT=[2,2], OFFSET=ipos), fpos[0], fpos[1])
               y = interpolate(ofile->VarGet('lat_rho', COUNT=[2,2], OFFSET=ipos), fpos[0], fpos[1])
            endif else begin
               x = interpolate(ofile->VarGet('x_rho', COUNT=[2,2], OFFSET=ipos), fpos[0], fpos[1])
               y = interpolate(ofile->VarGet('y_rho', COUNT=[2,2], OFFSET=ipos), fpos[0], fpos[1])
            endelse

            angle = 0
            if ofile->HasVar('angle') then $
               angle = interpolate(ofile->VarGet('angle', COUNT=[2,2], OFFSET=ipos), fpos[0],fpos[1])

            fmt = '(%"Point %s, %s at %s, %s, angle %s deg in a water depth of %s m")'
            message, /INFORM, string(FORMAT=fmt, mgh_format_float([position,x,y,angle*!radeg,h]))

         endelse

      end

      'station': begin

         if n_elements(nearest) gt 0 then $
            message, 'NEAREST keyword may not be specified for station data'

         ;; Get & check variable dimensions.

         dim_u = ofile->VarDims(variable[0])
         dim_v = ofile->VarDims(variable[1])

         if min([strlen(dim_u.station),strlen(dim_u.time)]) eq 0 then $
            message, 'The U variable must have station & time dimensions'
         if min([strlen(dim_v.station),strlen(dim_v.time)]) eq 0 then $
            message, 'The V variable must have station & time dimensions'
         if dim_u.time ne dim_v.time then $
            message, 'The U & V variables have different time dimensions'
         if dim_u.vertical ne dim_v.vertical then $
            message, 'The U & V variables have different vertical dimensions'

         ;; Establish position of point to be extracted and calculate
         ;; water depth, location, etc.

         if n_elements(position) eq 0 then position = 0

         if n_elements(position) ne 1 then $
            message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'position'

         if use_bath then begin
            h = ofile->VarGet('bath', COUNT=[1,1], OFFSET=[position,0])
         endif else begin
            h = ofile->VarGet('h', COUNT=[1], OFFSET=position)
         endelse

         angle = 0
         if ofile->HasVar('angle') then $
            angle = ofile->VarGet('angle', COUNT=[1], OFFSET=position)

         ipos = ofile->VarGet('Ipos', COUNT=[1], OFFSET=position)
         jpos = ofile->VarGet('Jpos', COUNT=[1], OFFSET=position)

         if lonlat then begin
            x = ofile->VarGet('lon_rho', COUNT=[1], OFFSET=position)
            y = ofile->VarGet('lat_rho', COUNT=[1], OFFSET=position)
         endif else begin
            x = ofile->VarGet('x_rho', COUNT=[1], OFFSET=position)
            y = ofile->VarGet('y_rho', COUNT=[1], OFFSET=position)
         endelse

         fmt = '(%"Station %d position %s, %s, angle %s deg at %s, %s in a water depth of %s m")'
         message, /INFORM, $
            string(FORMAT=fmt, position, mgh_format_float([ipos,jpos,angle*!radeg,x,y,h]))

      end

   endcase

   ;; Establish time variable name and records to be extracted. Get
   ;; time data

   n_time = ofile->DimInfo(dim_u.time, /DIMSIZE)
   mgh_resolve_indices, n_time, record_range

   case !true of
      ofile->HasVar(dim_u.time): $
         time_var = dim_u.time
      ofile->HasVar('ocean_time'): $
         time_var = 'ocean_time'
      ofile->HasVar('scrum_time'): $
         time_var = 'scrum_time'
      else: $
         message, 'Time variable not found'
   endcase

   time = ofile->VarGet(time_var, AUTOSCALE=0)

   time_units = {scale: 1D/(24D*3600D), offset: 0D}
   if ofile->HasAtt(time_var, 'units') then $
      time_units = mgh_dt_units(ofile->AttGet(time_var, 'units'))
   time *= time_units.scale
   time_offset = time_units.offset
   mgh_undefine, time_units

   if n_elements(time_range) gt 0 then begin
      record_range = mgh_subset(time, time_range)
   endif
   if n_elements(record_range) eq 0 then begin
      n_time = ofile->DimInfo(dim_u.time, /DIMSIZE)
      record_range = [0,n_time-1]
   endif

   rr0 = record_range[0]
   rr1 = record_range[1]
   rrn = rr1 - rr0 + 1

   time = time[rr0:rr1]

   ;; Establish level to be extracted

   dim_vertical = dim_u.vertical
   has_vertical = strlen(dim_vertical) gt 0

   if has_vertical then begin
      nk = (n_elements(depth) gt 0) + (n_elements(height) gt 0) + (n_elements(level) gt 0) + $
         (n_elements(sigma) gt 0) + keyword_set(profile)
      if nk gt 1 then $
         message, 'The DEPTH, HEIGHT, LEVEL, PROFILE and SIGMA keywords must not be specified together'
      if ofile->HasVar(dim_vertical) then begin
         svar = ofile->VarGet(dim_vertical)
      endif else begin
         case dim_vertical of
            's_rho': svar = ofile->VarGet('sc_r')
            's_w': svar = ofile->VarGet('sc_w')
         endcase
      endelse
      get_vertical = n_elements(depth) gt 0 || n_elements(height) gt 0 || n_elements(sigma) gt 0 || $
         size(level, /N_DIMENSIONS) gt 0 || keyword_set(profile)
      if get_vertical then begin
         theta_s = ofile->VarGet('theta_s')
         theta_b = ofile->VarGet('theta_b')
         hc = ofile->VarGet('hc')
         vstretch = ofile->HasVar('Vstretching') ? ofile->VarGet('Vstretching') : 1
         vtransform = ofile->HasVar('Vtransform') ? ofile->VarGet('Vtransform') : 1
         cs = mgh_roms_s_to_cs(svar, THETA_S=theta_s, THETA_B=theta_b, VSTRETCH=vstretch)
         z = mgh_roms_s_to_z(svar, h, HC=hc, CS=cs, VTRANSFORM=vtransform)
      endif else begin
         if n_elements(level) eq 0 then level = n_elements(svar) - 1
      endelse
   endif else begin
      if n_elements(depth) gt 0 then $
         message, 'The DEPTH keyword is not required or allowed ' + $
         'because the variables have no vertical dimension'
      if n_elements(height) gt 0 then $
         message, 'The HEIGHT keyword is not required or allowed ' + $
         'because the variables have no vertical dimension'
      if n_elements(level) gt 0 then $
         message, 'The LEVEL keyword is not required or allowed ' + $
         'because the variables have no vertical dimension'
      if keyword_set(profile) then $
         message, 'The PROFILE keyword may not be set because the ' + $
         'variables have no vertical dimension'
      if n_elements(sigma) gt 0 then $
         message, 'The SIGMA keyword is not required or allowed ' + $
         'because the variables have no vertical dimension'
      get_vertical = !false
   endelse

   ;; Get vector data

   case ftype of
      'history': begin
         if keyword_set(nearest) then begin
            offset = round(position)
            count = [1,1]
         endif else begin
            offset = floor(position)
            count = [2,2]
         endelse
         if has_vertical then begin
            if get_vertical then begin
               offset = [offset,0]
               count = [count,0]
            endif else begin
               offset = [offset,my_level]
               count = [count,1]
            endelse
         endif
         offset = [offset,rr0]
         count = [count,rrn]
         uv = reform(ofile->VectorGet(variable, OFFSET=offset, COUNT=count))
         if ~ keyword_set(nearest) then begin
            ;; Interpolate horizontally
            fpos = position - floor(position)
            c00 = (1-fpos[0])*(1-fpos[1])
            c10 = fpos[0]*(1-fpos[1])
            c01 = (1-fpos[0])*fpos[1]
            c11 = fpos[0]*fpos[1]
            case size(uv, /N_DIMENSIONS) of
               3: begin
                  uv = reform(c00*uv[0,0,*]+c10*uv[1,0,*]+c01*uv[0,1,*]+c11*uv[1,1,*])
               end
               4: begin
                  uv = reform(c00*uv[0,0,*,*]+c10*uv[1,0,*,*]+c01*uv[0,1,*,*]+c11*uv[1,1,*,*])
               end
            endcase
         endif
      end
      'station': begin
         offset = position
         count = 1
         if has_vertical then begin
            if get_vertical then begin
               offset = [0,offset]
               count = [0,count]
            endif else begin
               offset = [level,offset]
               count = [1,count]
            endelse
         endif
         offset = [offset,rr0]
         count = [count,rrn]
         uv = complex(reform(ofile->VarGet(variable[0], COUNT=count, OFFSET=offset)), $
            reform(ofile->VarGet(variable[1], COUNT=count, OFFSET=offset)))
      end
   endcase

   ;; Handle the vertical dimension

   case !true of
      (~ has_vertical): begin
         result ={time: time, time_offset: time_offset, x: x, y: y, h: h, angle: angle, uv: uv}
      end
      keyword_set(profile): begin
         result ={time: time, time_offset: time_offset, x: x, y: y, h: h, angle: angle, z: z, uv: uv}
      end
      n_elements(level) gt 0: begin
         if get_vertical then begin
            result ={time: time, time_offset: time_offset, x: x, y: y, h: h, angle: angle, level: level, uv: uv[level,*]}
         endif else begin
            result ={time: time, time_offset: time_offset, x: x, y: y, h: h, angle: angle, level: level, uv: uv}
         endelse
      end
      n_elements(depth) gt 0: begin
         kk = fltarr(n_elements(depth))
         for k=0,n_elements(depth)-1 do begin
            ;; Return surface data for non-finite depth values
            kk[k] = finite(depth[k]) ? mgh_locate(z, XOUT=[-depth[k]]) : n_elements(svar)-1
         endfor
         uv = interpolate(uv, kk, findgen(rrn), /GRID)
         if size(depth, /N_DIMENSIONS) eq 0 then uv = reform(uv)
         result ={time: time, time_offset: time_offset, x: x, y: y, h: h, angle: angle, depth: depth, uv: uv}
      end
      n_elements(height) gt 0: begin
         kk = fltarr(n_elements(height))
         for k=0,n_elements(height)-1 do begin
            ;; Return surface data for non-finite height values
            kk[k] = finite(height[k]) ? mgh_locate(z, XOUT=[height[k]-h]) : n_elements(svar)-1
         endfor
         uv = interpolate(uv, kk, findgen(rrn), /GRID)
         if size(height, /N_DIMENSIONS) eq 0 then uv = reform(uv)
         result ={time: time, time_offset: time_offset, x: x, y: y, h: h, angle: angle, height: height, uv: uv}
      end
      n_elements(sigma) gt 0: begin
         kk = fltarr(n_elements(sigma))
         for k=0,n_elements(sigma)-1 do begin
            ;; Return surface data if sigma = 1, bottom data if sigma = 0
            if sigma[k] eq -1 then begin
               kk[k] = 0
            endif else if sigma[k] eq 0 then begin
               kk[k] = n_elements(svar)-1
            endif else begin
               kk[k] = mgh_locate(z, XOUT=[(sigma[k]-1)*h])
            endelse
         endfor
         uv = interpolate(uv, kk, findgen(rrn), /GRID)
         if size(sigma, /N_DIMENSIONS) eq 0 then uv = reform(uv)
         result ={time: time, x: x, y: y, h: h, angle: angle, sigma: sigma, uv: uv}
      end
   endcase

   message, /INFORM, 'Saving ROMS time series data to '+tmpfile
   save, FILE=tmpfile, result

   return, result

end
