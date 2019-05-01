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
;     A reference to an MGHromsRegrid object.
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
;   Mark Hadfield, 2019-05:
;     Written, based on mgh_roms_series_vector
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

   if ~ obj_isa(ofile, 'MGHromsRegrid') then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_objwrgclass', 'ofile'

   ;; Process VARIABLE argument

   if n_elements(variable) eq 0 then begin
      has_uv = ofile->HasVar('u') && ofile->HasVar('v')
      variable = has_uv ? ['u','v'] : ['ubar','vbar']
   endif

   if n_elements(variable) ne 2 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'variable'

   ;; Construct a useful temporary file name. Quite an operation!

   ofile->GetProperty, FILE_NAME=fname

   str_file = mgh_format_integer(mgh_hashcode(fname))
   str_var = mgh_format_integer(mgh_hashcode(variable))

   tmpfile = ['mgh_roms_regrid_series_vector',str_file,str_var]
   if n_elements(location) gt 0 then $
      tmpfile = [tmpfile,'location',mgh_format_float(location)]
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

   ;; Get & check variable dimensions.

   dim = ofile->VarDimNames(variable[0])

   if n_elements(dim) ne 3 then $
      message, 'Unexpected number of dimensions

   if ~ array_equal(dim, ofile->VarDimNames(variable[1])) then $
      message, 'U and V dimensions are mismatched'

   ;; Retrieve longitude and latitude and determine the grid position
   ;; of the specified location

   lon = ofile->VarGet('lon')
   lat = ofile->VarGet('lat')

   case !true of

      size(lon, /N_DIMENSIONS) eq 1 && size(lat, /N_DIMENSIONS) eq 1: begin
         lonlat2d = !false
         pos = [mgh_locate(lon, XOUT=location[0]), $
                mgh_locate(lat, XOUT=location[1])]
      end

      size(lon, /N_DIMS) eq 2 && size(lat, /N_DIMS) eq 2: begin
         if ~ array_equal(size(lon, /DIMENSIONS), size(lat, /DIMENSIONS)) then $
            message, 'Mismatched dimensions'
         lonlat2d = !true
         pos = mgh_locate2(lon, lat, XOUT=location[0], YOUT=location[2])
      end

   endcase

   if keyword_set(nearest) then begin

      ipos = round(pos)

      if lonlat2d then begin
         x = lon[ipos[0],ipos[1]]
         y = lat[ipos[0],ipos[1]]
      endif else begin
         x = lon[ipos[0]]
         y = lat[ipos[1]]
      endelse

      fmt = '(%"Nearest point %d, %d at %s, %s")'
      message, /INFORM, string(FORMAT=fmt, ipos, mgh_format_float([x,y]))

   endif else begin

      if lonlat2d then begin
         x = interpolate(lon, pos[0], pos[1])
         y = interpolate(lat, position[0], position[1])
      endif else begin
         x = interpolate(lon, pos[0])
         y = interpolate(lat, pos[1])
      endelse

      ;;; This should print the location originally specified
      fmt = '(%"Point %s, %s at %s, %s")'
      message, /INFORM, string(FORMAT=fmt, mgh_format_float([pos,x,y]))

   endelse

   ;; Gettime data

   n_time = ofile->DimInfo(dim[2], /DIMSIZE)
   mgh_resolve_indices, n_time, record_range

   case !true of
      ofile->HasVar(dim[2]): $
         time_var = dim[2]
      ofile->HasVar('ocean_time'): $
         time_var = 'ocean_time'
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

   ;; Get vector data

   if keyword_set(nearest) then begin
      offset = round(pos)
      count = [1,1]
   endif else begin
      offset = floor(pos)
      count = [2,2]
   endelse
   offset = [offset,rr0]
   count = [count,rrn]
   uv = complex(reform(ofile->VarGet(variable[0], OFFSET=offset, COUNT=count)), $
                reform(ofile->VarGet(variable[1], OFFSET=offset, COUNT=count)))
   if ~ keyword_set(nearest) then begin
      ;; Interpolate horizontally
      fpos = pos - floor(pos)
      c00 = (1-fpos[0])*(1-fpos[1])
      c10 = fpos[0]*(1-fpos[1])
      c01 = (1-fpos[0])*fpos[1]
      c11 = fpos[0]*fpos[1]
      uv = reform(c00*uv[0,0,*]+c10*uv[1,0,*]+c01*uv[0,1,*]+c11*uv[1,1,*])
   endif

   result ={time: time, x: x, y: y, uv: uv}

   message, /INFORM, 'Saving ROMS time series data to '+tmpfile
   save, FILE=tmpfile, result

   return, result

end
