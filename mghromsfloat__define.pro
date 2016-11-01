;+
; CLASS NAME:
;   MGHromsFloat
;
; PURPOSE:
;   This class wraps a series of SCRUM/ROMS floats files.
;
; OBJECT CREATION SEQUENCE
;   obj = obj_new('MGHromsFloat', Files)
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 1999-08:
;     Written.
;   Mark Hadfield, 2012-01:
;     Removed GRID_DESTROY functionality
;   Mark Hadfield, 2012-05:
;     For the VarGetFloat method, AUTOSCALE now defaults to 1.
;     This is required to handle packed float data.
;   Mark Hadfield, 2016-10:
;     - Removed obsolete references to keyword arguments in the HasVar method.
;     - Removed support for STRUCT variables in HasVar and VarGet.
;-
function MGHromsFloat::Init, files, GRID_FILE=grid_file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR

   case size(grid_file, /TNAME) of
      'OBJREF': begin
         self.grid_file = grid_file
      end
      'STRING': begin
         self.grid_file = obj_new('MGHromsHistory', grid_file)
      end
      'UNDEFINED':
      else: message, 'The argument is of the wrong data type'
   endcase

   return, self->MGHncSequence::Init(files, _STRICT_EXTRA=extra)

end


pro MGHromsFloat::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR

end


pro MGHromsFloat::GetProperty, GRID_FILE=grid_file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR

   grid_file = self.grid_file

   self->MGHncSequence::GetProperty, _STRICT_EXTRA=extra

end

function MGHromsFloat::TimeVarName

   compile_opt DEFINT32
   compile_opt STRICTARR

   vars = self->VarNames()

   case 1 of
      max(vars eq 'ocean_time') gt 0: return, 'ocean_time'
      else: return, ''
   endcase

end

; MGHromsFloat::GetCppOptions
;
;   Return a string array containing the list of CPP options
;
function MGHromsFloat::GetCppOptions, COUNT=count

   compile_opt DEFINT32
   compile_opt STRICTARR

   att = 'CPP-options'

   if self->HasAtt(att, /GLOBAL) then begin
      s = self->AttGet(att, /GLOBAL)
      s = strtrim(strsplit(s, ',', /EXTRACT), 2)
      count = n_elements(s)
      return, s
   endif else begin
      count = 0
      return, ''
   endelse

end

; MGHromsFloat::GetZGrid
;
;   Return a 3-D array of Z values on the model grid
;
function MGHromsFloat::GetZGrid, HORIZONTAL=horizontal, VERTICAL=vertical

   compile_opt DEFINT32
   compile_opt STRICTARR

   if n_elements(horizontal) eq 0 then horizontal = 'rho'

   if n_elements(vertical) eq 0 then vertical = 'rho'

   ;; Read/calculate depths on horizontal grid

   h = self.grid_file->VarGet('h')

   h = mgh_roms_stagger(h, TO=horizontal)

   numh = size(h, /DIMENSIONS)

   ;; Read s-coordinate vectors

   case strlowcase(vertical) of
      'rho': s = self->VarGet('sc_r')
      'w'  : s = (self->VarGet('sc_w') < 0) > (-1)
   endcase

   nums = n_elements(s)

   ;; Read s-coordinate parameters

   theta_s = self->VarGet('theta_s')
   theta_b = self->VarGet('theta_b')
   hc = self->VarGet('hc')
   vstretch = self->HasVar('Vstretching') ? self->VarGet('Vstretching') : 1
   vtransform = self->HasVar('Vtransform') ? self->VarGet('Vtransform') : 1

   ;; Create output array and calculate heights.

   result = fltarr(numh[0], numh[1], nums)

   cs = mgh_roms_s_to_cs(s, THETA_S=theta_s, THETA_B=theta_b, VSTRETCH=vstretch)
   for i=0,numh[0]-1 do for j=0,numh[1]-1 do $
        result[i,j,*] = mgh_roms_s_to_z(s, h[i,j,*], CS=cs, HC=hc, VTRANSFORM=vtransform)

   return, result

end

; MGHromsFloat::HasVar
;
function MGHromsFloat::HasVar, varname

   compile_opt DEFINT32
   compile_opt STRICTARR

   case size(varname, /TNAME) of

      'STRING': begin

         if self->MGHncSequence::HasVar(VarName) then return, !true

         if obj_valid(self.grid_file) then begin

            if self.grid_file->HasVar(VarName) then return, !true

            if max(strmatch(['x','y','lon','lat'],varname)) gt 0 then begin

               if ~ self->HasVar('Xgrid') then return, !false
               if ~ self->HasVar('Ygrid') then return, !false

               return, !true

            endif

         endif

         return, !false

      end

      else: begin

         message, BLOCK='mgh_mblk_motley', NAME='wrongtype', 'varname'

      end

   endcase

end

; MGHromsFloat::VarGet
;
function MGHromsFloat::VarGet, var, AUTOSCALE=autoscale, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR

   case 1B of

      isa(var, 'STRUCT'): begin
         message, 'No structure-type varnames have been implemented yet'
      end

      isa(var, 'STRING') && strmatch(var, '*time*'): begin
         result = self->MGHncSequence::VarGet(var, $
                                              AUTOSCALE=0, _STRICT_EXTRA=extra)
         ;; If the AUTOSCALE keyword is set, apply the scale factor
         ;; specified in the time variable's "units" attribute. Ignore the
         ;; offset as this is problematic.
         if keyword_set(autoscale) then begin
            time_units = {scale: 1D/(24D*3600D), offset: 0D}
            if self->MGHncSequence::HasAtt(var, 'units') then $
                 time_units = mgh_dt_units(self->AttGet(var, 'units'))
            result = result*time_units.scale
         endif
         return, result
      end

      isa(var, 'STRING'): begin
         if self->MGHncSequence::HasVar(var) then $
              return, self->MGHncSequence::VarGet(var, AUTOSCALE=autoscale, _STRICT_EXTRA=extra)
         if obj_valid(self.grid_file) then begin
            if self.grid_file->HasVar(var) then $
                 return, self.grid_file->VarGet(var, AUTOSCALE=autoscale, _STRICT_EXTRA=extra)
            if max(strmatch(['x','y','lon','lat'], var)) gt 0 then begin
               pos = self->VarGet(var+'_rho')
               xgrid = self->VarGet('Xgrid', AUTOSCALE=autoscale, _STRICT_EXTRA=extra)
               ygrid = self->VarGet('Ygrid', AUTOSCALE=autoscale, _STRICT_EXTRA=extra)
               return, mgh_interpolate(pos, xgrid, ygrid, MISSING=!values.f_nan)
            endif
         endif
         message, 'Variable not found: programmer error!'
      end

   endcase

end

; MGHromsFloat::VarGetFloat
;
;   Specialised "varget" function for float data
;
function MGHromsFloat::VarGetFloat, varname, $
     AUTOSCALE=autoscale, $
     FLOAT_RANGE=float_range, $
     FLOAT_STRIDE=float_stride, $
     FLOATS=floats, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     RECORDS=records

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Check variable

   if ~ self->HasVar(varname) then $
        message, 'Variable not found'

   var_dims = self->VarDimNames(varname, COUNT=n_var_dims)

   if n_var_dims ne 2 then $
        message, 'This method works only for 2D variables'

   ;; Process keywords

   if n_elements(autoscale) eq 0 then autoscale = 1

   ;; Select floats & records. The first dimension should be "drifter"
   ;; and is interpreted as a float index; the second dimension is
   ;; interpreted as time

   n = self->DimInfo(var_dims[0], /DIMSIZE)

   mgh_resolve_indices, n, float_range, float_stride, floats

   n_floats = n_elements(floats)

   n = self->DimInfo(var_dims[1], /DIMSIZE)

   mgh_resolve_indices, n, record_range, record_stride, records

   n_records = n_elements(records)

   ;; Construct output array

   template = self->VarGet(varname, AUTOSCALE=autoscale, COUNT=[1,1])

   result = replicate(template[0], n_floats, n_records)

   for r=0,n_elements(records)-1 do begin

      rec = records[r]

      ;; Get float position data

      f0 = min(floats, MAX=f1)

      data = self->VarGet(varname, AUTOSCALE=autoscale, $
                          COUNT=[f1-f0+1,1], OFFSET=[f0,rec])

      ;; Assigning multiple values to a single array location.
      ;; Values are written along the matching dimension. See article
      ;; on comp.lang.idl-pvwave by Wayne Landsman on 2002-02-05 entitled
      ;; "RE: Fast Shear"
      result[0,r] = data[floats-f0]

   endfor

   return, result


end

pro MGHromsFloat__Define

   compile_opt DEFINT32
   compile_opt STRICTARR

   struct_hide, {MGHromsFloat, inherits MGHncSequence, grid_file: obj_new()}

end
