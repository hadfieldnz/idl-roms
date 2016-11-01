;+
; CLASS NAME:
;   MGHromsStation
;
; PURPOSE:
;   This class wraps a ROMS station files. Each MGHromsStation object
;   has associated with it an optional ROMS grid file, which can be
;   searched for grid information if the station files lack it.
;
; OBJECT CREATION SEQUENCE
;   obj = obj_new('MGHromsStation', Files)
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2001-03:
;     Written.
;   Mark Hadfield, 2010-07:
;     The VarGet method's AUTOSCALE keyword now defaults to 1, as
;     for the MGHromsHistory class.
;   Mark Hadfield, 2012-04:
;     - Removed references to grid_destroy: unnecessary with automatic
;       garbage collection.
;     - Fixed an apparently long-standing bug: Cleanup method did not call
;       MGHncSequence::Cleanup.
;     - Added support for summed variables (eg. "dye_01+dye_02").
;   Mark Hadfield, 2012-04:
;     - Begin adding supporting for sediment (the bed dimension).
;-
function MGHromsStation::Init, files, GRID_FILE=grid_file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

    case size(grid_file, /TNAME) of
        'OBJREF': begin
            self.grid_file = grid_file
        end
        'STRING': begin
            self.grid_file = obj_new('MGHncReadFile', grid_file)
        end
        'UNDEFINED':
        else: message, 'The argument is of the wrong data type'
    endcase

    return, self->MGHncSequence::Init(files, _STRICT_EXTRA=extra)

end


pro MGHromsStation::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGHncSequence::Cleanup

end


; MGHromsStation::GetCppOptions
;
;   Return a string array containing the list of CPP options
;
function MGHromsStation::GetCppOptions, COUNT=count

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  s = self->AttGetData('CPP-options', EXISTS=exists)
  if ~ exists then $
    s = self->AttGetData('CPP_options', EXISTS=exists)
  if ~ exists then $
    message, 'Cannot find attribute'

  if exists then begin
    result = strtrim(strsplit(s, ',', /EXTRACT, COUNT=count), 2)
  endif else begin
    count = 0
    result = ''
  endelse

  return, result

end

; MGHromsHistory::GetScoord
;
function MGHromsStation::GetScoord, dim_vertical

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(dim_vertical) eq 0 then dim_vertical = 's_rho'

   ;; Retrieve s-coordinate parameters and s values for the specified
   ;; vertical dimension. The following code is complicated by a the various
   ;; flavours of ROMS output files:
   ;;  - In older files, the bottom value of s_w is omitted.
   ;;  - In ROMS 2.1 and earlier, the names of the s-coordinate variables did
   ;;    not match the corresponding dimensions--heavens knows why. In ROMS 2.2
   ;;    this was fixed.
   ;;  - In ROMS-AGRIF files, theta_s, theta_b and hc are stored as global
   ;;    attributes, not as variables
   ;;  - In Nicolas Gruber's version of ROMS, the s coordinate values are not stored.

   if strlen(dim_vertical) eq 0 then return, !null

   if ~ self->HasDim('s_rho') then return, !null

   result = dictionary()

   n_s_rho = self->DimInfo('s_rho', /DIMSIZE)
   result.n = n_s_rho + (dim_vertical eq 's_w')

   result.theta_s = self->HasVar('theta_s') ? self->VarGet('theta_s') : self->AttGet('theta_s', /GLOBAL)
   result.theta_b = self->HasVar('theta_b') ? self->VarGet('theta_b') : self->AttGet('theta_b', /GLOBAL)

   result.hc = self->HasVar('hc') ? self->VarGet('hc') : self->AttGet('hc', /GLOBAL)

   result.vstretch = self->HasVar('Vstretching') ? self->VarGet('Vstretching'): 1

   result.vtransform = self->HasVar('Vtransform') ? self->VarGet('Vtransform'): 1

   if self->HasVar(dim_vertical) then begin
      s = self->VarGet(dim_vertical)
   endif else begin
      case dim_vertical of
         's_rho': begin
            if self->HasVar('sc_r') then begin
               s =  self->VarGet('sc_r')
            endif else begin
               s = mgh_stagger(mgh_range(-1, 0, STRIDE=1.0/n_s_rho), DELTA=-1)
            endelse
         end
         's_w': begin
            if self->HasVar('sc_w') then begin
               s = self->VarGet('sc_w')
            endif else begin
               s = mgh_range(-1, 0, STRIDE=1.0D/n_s_rho)
            endelse
            if n_elements(s) eq n_s_rho then s = [-1,s]
         end
      endcase
   endelse
   ;; Clip s-coordinate values; this may be necessary
   ;; when netCDF data have been packed.
   result.s = (s > (-1)) < 0

   return, result->ToStruct(/RECURSIVE)

end

; MGHromsStation::TimeVarName
;
function MGHromsStation::TimeVarName, var

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   for i=0,n_elements(var)-1 do begin
      if self->HasVar(var[i]) then return, var[i]
   endfor

   if self->HasVar('ocean_time') then return, 'ocean_time'

   if self->HasVar('scrum_time') then return, 'scrum_time'

   return, !null

end
; MGHromsStation::GetZGrid
;
;   Return an array of Z values at the station locations
;
function MGHromsStation::GetZGrid, VERTICAL=vertical, ZETA=zeta

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  if n_elements(vertical) eq 0 then vertical = 'rho'

  n_sta = self->DimInfo('station', /DIMSIZE)

  if n_elements(zeta) eq 0 then zeta = fltarr(n_sta)

  h = self->VarGet('h')

  ;; Read s-coordinate vectors. Although s-coordinate data at w levels (sc_w)
  ;; is usually stored in ROMS files, the names & dimensions vary between
  ;; file types, so use only rho-level info.

  if strlowcase(vertical) eq 'w' then begin
    n_s_rho = self->DimInfo('s_rho', /DIMSIZE)
    s = self->HasVar('s_w') ? self->VarGet('s_w') : self->VarGet('sc_w')
    if n_elements(s) eq n_s_rho then s = [-1,s]
    s = (s < 0) > (-1)
  endif else begin
    s= self->VarGet('s_rho')
  endelse

  n_s = n_elements(s)

  ;; Read s-coordinate parameters

  theta_s = self->VarGet('theta_s')
  theta_b = self->VarGet('theta_b')
  hc = self->VarGet('hc')
  vstretch = self->VarGet('Vstretching') ? self->VarGet('Vstretching'): 1
  vtransform = self->VarGet('Vtransform') ? self->VarGet('Vtransform'): 1

  ;; Create output array and calculate heights.

  result = fltarr([n_s,n_sta])

  cs = mgh_roms_s_to_cs(s, THETA_S=theta_s, THETA_B=theta_b, VSTRETCH=vstretch)

  for i=0,n_sta-1 do begin
    result[*,i] = mgh_roms_s_to_z(s, h[i], ZETA=zeta[i], CS=cs, HC=hc, VTRANSFORM=vtransform)
  endfor

  return, result

end

; MGHromsStation::VarDimNames
;
function MGHromsStation::VarDimNames, var, COUNT=count

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(var) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'var'

   vartype = size(var, /TYPE)

   case 1B of

      isa(var, 'STRING') && strmatch(var, '*+*'): begin
         vars = strsplit(var, '+', /EXTRACT)
         result = self->MGHncSequence::VarDimNames(vars[0], COUNT=count)
      end

      isa(var, 'STRING') && (var eq 'sbar'): begin
         result = self->MGHncSequence::VarDimNames('ubar', COUNT=count)
      end

      isa(var, 'STRING') && (var eq 'spd'): begin
         result = self->MGHncSequence::VarDimNames('u', COUNT=count)
      end

      isa(var, 'STRING') && strmatch(var, 'uv@*'): begin
        result = self->MGHncSequence::VarDimNames('u', COUNT=count)
      end

      isa(var, 'STRING') && (var eq 'Sbot'): begin
        result = self->MGHncSequence::VarDimNames('Ubot', COUNT=count)
      end

      isa(var, 'STRING') && (var eq 'Dsbl'): begin
         result = self->MGHncSequence::VarDimNames('Hsbl', COUNT=count)
      end

      isa(var, 'STRING') && (var eq 'Dbbl'): begin
         result = self->MGHncSequence::VarDimNames('Hbbl', COUNT=count)
      end

      isa(var, 'STRING'): begin
         result = self->MGHncSequence::VarDimNames(var, COUNT=count)
      end

      else: message, 'var argument invalid'

   endcase

   return, result

end

; MGHromsStation::VarDims

function MGHromsStation::VarDims, var

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Get dimension names as a string array.

   dim = self->VarDimNames(var)

   ;; Set up the output structure

   result = {vertical: '', station: '', time: '', bed: ''}

   ;; Load variable dimensions into output

   for d=0,n_elements(dim)-1 do begin
      if strmatch(dim[d], 'station', /FOLD_CASE) then $
           result.station[0] = dim[d]
      if strmatch(dim[d], '*time*', /FOLD_CASE) then $
           result.time = dim[d]
      if strmatch(dim[d], 's_*', /FOLD_CASE) then $
           result.vertical = dim[d]
      if strmatch(dim[d], '*bed', /FOLD_CASE) then $
           result.bed = dim[d]
   endfor

   return, result

end

; MGHromsStation::VarGet
;
function MGHromsStation::VarGet, var, $
     AUTOSCALE=autoscale, COUNT=count, OFFSET=offset, _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(var) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'var'

   if n_elements(autoscale) eq 0 then autoscale = 1B

   vartype = size(var, /TYPE)

   case 1B of

     isa(var, 'STRING') && strmatch(var, '*time*'): begin
       result = self->MGHncSequence::VarGet(var, AUTOSCALE=0, $
         COUNT=count, OFFSET=offset, $
         _STRICT_EXTRA=_extra)
       ;; If the AUTOSCALE keyword is set, apply the scale factor
       ;; specified in the time variable's "units" attribute. Ignore the
       ;; offset as this is problematic.
       if keyword_set(autoscale) then begin
         time_units = {scale: 1D/(24D*3600D), offset: 0D}
         if self->MGHncSequence::HasAtt(var, 'units') then $
           time_units = mgh_dt_units(self->AttGet(var, 'units'))
         result = result*time_units.scale
       endif
     end

      isa(var, 'STRING') && strmatch(var, '*+*'): begin
         vars = strsplit(var, '+', /EXTRACT)
         result = 0
         foreach v,vars do begin
            result += self->VarGet(v, AUTOSCALE=autoscale, $
                                   COUNT=count, OFFSET=offset)
         endforeach
      end

      isa(var, 'STRING') && (var eq 'sbar'): begin
         ubar = self->MGHncSequence::VarGet('ubar', AUTOSCALE=autoscale, $
                                            COUNT=count, OFFSET=offset, $
                                            _STRICT_EXTRA=_extra)
         vbar = self->MGHncSequence::VarGet('vbar', AUTOSCALE=autoscale, $
                                            COUNT=count, OFFSET=offset, $
                                            _STRICT_EXTRA=_extra)
         result = sqrt(temporary(ubar)^2 + temporary(vbar)^2)
      end

      isa(var, 'STRING') && (var eq 'spd'): begin
         u = self->MGHncSequence::VarGet('u', AUTOSCALE=autoscale, $
                                         COUNT=count, OFFSET=offset, $
                                         _STRICT_EXTRA=_extra)
         v = self->MGHncSequence::VarGet('v', AUTOSCALE=autoscale, $
                                         COUNT=count, OFFSET=offset, $
                                         _STRICT_EXTRA=_extra)
         result = sqrt(temporary(u)^2 + temporary(v)^2)
      end

      isa(var, 'STRING') && strmatch(var, 'uv@*'): begin
        cj = complex(0, 1)
        ;; Velocity component direction, geographic convention, degrees.
        direction = float(strmid(var, 3))
        u = self->MGHncSequence::VarGet('u', AUTOSCALE=autoscale, COUNT=count, OFFSET=offset, _STRICT_EXTRA=_extra)
        v = self->MGHncSequence::VarGet('v', AUTOSCALE=autoscale, COUNT=count, OFFSET=offset, _STRICT_EXTRA=_extra)
        uv = complex(temporary(u), temporary(v))
        if self->HasVar('angle') then begin
          angle = self->VarGet('angle', COUNT=count[1], OFFSET=offset[1])
        endif else begin
          angle = 0
        endelse
        n_angle = n_elements(angle)
        if n_angle gt 1 then begin
          for i=0,n_angle-1 do uv[*,i,*] *= exp(cj*(angle[i]-!dtor*(90-direction)))
        endif else begin
          uv  *= exp(cj*(angle-!dtor*(90-direction)))
        endelse
        return, real_part(uv)
      end

      isa(var, 'STRING') && (var eq 'Sbot'): begin
        ubot = self->MGHncSequence::VarGet('Ubot', AUTOSCALE=autoscale, $
          COUNT=count, OFFSET=offset, $
          _STRICT_EXTRA=_extra)
        vbot = self->MGHncSequence::VarGet('Vbot', AUTOSCALE=autoscale, $
          COUNT=count, OFFSET=offset, $
          _STRICT_EXTRA=_extra)
        result = sqrt(temporary(ubot)^2 + temporary(vbot)^2)
      end

      isa(var, 'STRING') && (var eq 'Dsbl'): begin
         result = - self->MGHncSequence::VarGet('Hsbl', AUTOSCALE=autoscale, $
                                                COUNT=count, OFFSET=offset, $
                                                _STRICT_EXTRA=_extra)
      end

      isa(var, 'STRING') && (var eq 'Dbbl'): begin
         my_offset = n_elements(offset) gt 0 ? offset : [0,0]
         my_count = n_elements(count) gt 0 ? count : [0,0]
         result = self->MGHncSequence::VarGet('Hbbl', $
                                              AUTOSCALE=autoscale, COUNT=my_count, $
                                              OFFSET=my_offset, _STRICT_EXTRA=_extra)
         h = self->MGHncSequence::VarGet('h', $
                                         AUTOSCALE=autoscale, COUNT=my_count[0], $
                                         OFFSET=my_offset[0], _STRICT_EXTRA=_extra)
         if n_elements(h) gt 1 then begin
            h = rebin(reform([h], [size(h, /DIMENSIONS),1]), $
                      [size([h], /DIMENSIONS),my_count[1]])
         endif
         result += h
      end

      vartype eq 7: begin
         result = self->MGHncSequence::VarGet(var, AUTOSCALE=autoscale, $
                                              COUNT=count, OFFSET=offset, $
                                              _STRICT_EXTRA=_extra)
      end

      else: message, 'var argument invalid'

   endcase

   return, result

end

pro MGHromsStation__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        { MGHromsStation, inherits MGHncSequence, $
          grid_file: obj_new()}

end
