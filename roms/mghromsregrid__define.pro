;+
; CLASS NAME:
;   MGHromsRegrid
;
; PURPOSE:
;   This class wraps a series of ROMS forcing files with 1D or 2D
;   longitude & latitude data, which are regridded by ROMS
;
; OBJECT CREATION SEQUENCE
;   obj = obj_new('MGHromsRegrid', Files)
;
; SLICES:
;   Several of the methods below are used to support extraction of
;   data on spatial slices through the ROMS grid. The slice types are:
;
;     L-slice
;       A 2D horizontal slice on a rectilinear lon-lat grid as used for forcing
;       data.
;
; PATCHES:
;   A Patch is a polygonal area over which integrals (and maybe other statistics)
;   can be calculated. Its vertices are defined in physical space.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2019-04:
;     Written as a severely stripped-down form of MGHromsHistory.
;-
;+
; METHOD NAME:
;   MGHromsRegrid::VarGet
;
; CALLING SEQUENCE:
;   slice = history->VarGet(arg)
;
; POSITIONAL PARAMETERS:
;   arg (input, scalar string)
;     The name of a variable. In addition to the variables actually in
;     the history file, this class's VarGet method supports several
;     synthetic variables, eg "spd" for speed, calculated from
;     variables u and v.
;
; RETURN VALUE:
;   A numeric or character
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2000-12:
;       Written.
;   Mark Hadfield, 2008-01:
;       The AUTOSCALE keyword is now intercepted by this method and
;       passed to MGHncSequence::VarGet appropriately. This was done
;       so that auto-scaling of time variables could be enabled or
;       disabled using this switch.
;-
;+
; METHOD NAME:
;   MGHromsRegrid::VarDims
;
; CALLING SEQUENCE:
;   slice = history->VarDims(arg)
;
; POSITIONAL PARAMETERS:
;   arg (input, scalar string)
;     The name of a variable in the history file. Currently this has
;     to be a space-time variable, ie. one with 2 horizontal
;     dimensions, an optional vertical dimension and an optional time
;     dimension, in that order. Examples of this are h, zeta, temp,
;     u. In future, this requirement may be relaxed and the method
;     made more general.
;
; RETURN VALUE:
;   The method returns a structure with tags "horizontal", "vertical"
;   and "time" containing the respective dimensions.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, Dec 2000:
;       Written.
;-
function MGHromsRegrid::Init, files, _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   return, self->MGHncSequence::Init(files, _STRICT_EXTRA=_extra)

end

pro MGHromsRegrid::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGHncSequence::Cleanup

end

; MGHromsRegrid::HasVar
;
function MGHromsRegrid::HasVar, var

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case !true of

      isa(var, 'STRING') && strmatch(var, 'abs(*,*)'): begin
         pp = stregex(var, '(^abs\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         result = self->HasVar(vv[0]) && self->HasVar(vv[1])
      end

      isa(var, 'STRING') && strmatch(var, 'harmonic_amp(*)'): begin
         pp = stregex(var, '(^harmonic_amp\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         result = self->HasVar(vv[0]) && self->HasVar(vv[1])
      end

      isa(var, 'STRING') && strmatch(var, 'harmonic_pha(*)'): begin
         pp = stregex(var, '(^harmonic_pha\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         result = self->HasVar(vv[0]) && self->HasVar(vv[1])
      end

      isa(var, 'STRING') && strmatch(var, 'tide_scalar_amp(*)'): begin
         pp = stregex(var, '(^tide_scalar_amp\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         result = self->HasVar(vv[0]) && self->HasVar(vv[1])
      end

      isa(var, 'STRING') && strmatch(var, 'tide_scalar_pha(*)'): begin
         pp = stregex(var, '(^tide_scalar_pha\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         result = self->HasVar(vv[0]) && self->HasVar(vv[1])
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_sma(*,*)'): begin
         pp = stregex(var, '(^tide_vector_sma\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         result = self->HasVar(vv[0]+'_cos') && self->HasVar(vv[0]+'_sin') && $
                  self->HasVar(vv[1]+'_cos') && self->HasVar(vv[1]+'_sin')
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_ecc(*,*)'): begin
         pp = stregex(var, '(^tide_vector_ecc\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         result = self->HasVar(vv[0]+'_cos') && self->HasVar(vv[0]+'_sin') && $
            self->HasVar(vv[1]+'_cos') && self->HasVar(vv[1]+'_sin')
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_inc(*,*)'): begin
         pp = stregex(var, '(^tide_vector_inc\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         result = self->HasVar(vv[0]+'_cos') && self->HasVar(vv[0]+'_sin') && $
            self->HasVar(vv[1]+'_cos') && self->HasVar(vv[1]+'_sin')
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_pha(*,*)'): begin
         pp = stregex(var, '(^tide_vector_pha\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         result = self->HasVar(vv[0]+'_cos') && self->HasVar(vv[0]+'_sin') && $
            self->HasVar(vv[1]+'_cos') && self->HasVar(vv[1]+'_sin')
      end

      else: begin
         result = self->MGHncSequence::HasVar(var)
      end

   endcase

   return, result

end

; ** L-slice methods *******************************************
;
; MGHromsRegrid::LsliceData
;
function MGHromsRegrid::LsliceData, var, $
   GRID=grid, I_RANGE=i_range, J_RANGE=j_range, RECORD=record

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Check for problems with inputs

   if n_elements(var) ne 1 then $
      message, 'The name of a variable must be supplied'

   if ~ (isa(var, 'STRING') || isa(var, 'STRUCT'))  then $
      message, 'A variable identifier must be supplied'

   ;; If no grid information is supplied, then get it. Otherwise check
   ;; grid is consistent with the current variable name & function
   ;; arguments.

   if n_elements(grid) eq 0 then $
      grid = self->LsliceGrid(I_RANGE=i_range, J_RANGE=j_range)

   ;; Get the names of the time dimension & variable, if any, associated
   ;; with the variable. This code requires that any variable with a time
   ;; dimension must have a time attribute.

   if self->HasAtt(var, 'time') then begin
      var_hastime = !true
      var_time = self->AttGet(var, 'time')
      self->VarInfo, var, DIM_NAMES=dim_time
   endif else begin
      var_hastime = !false
   endelse

   ;; Check consistency of function arguments with variable dimensions

   if var_hastime then begin
      if n_elements(record) eq 0 then record = 0
      if record lt 0 then $
         record = self->DimInfo(dim_time, /DIMSIZE) + record
   endif else begin
      if n_elements(record) gt 0 then begin
         message, 'The RECORD keyword is not required or allowed when ' + $
            'the variable '+var+' has no time dimension'
      endif
   endelse

   ;; Abbreviations for i and j ranges

   ira0 = grid.i_range[0]
   ira1 = grid.i_range[1]
   iran = ira1-ira0+1
   jra0 = grid.j_range[0]
   jra1 = grid.j_range[1]
   jran = jra1-jra0+1

   ;; Read & unpack variable

   result = self->VarGet(var, OFFSET=[ira0,jra0,record], COUNT=[iran,jran,1], /AUTOSCALE)

   return, result

end

; MGHromsRegrid::LsliceGrid
;
function MGHromsRegrid::LsliceGrid, $
   I_RANGE=i_range, J_RANGE=j_range

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Get dimensions of the lon-lat grid

   vinfo_lon = self->VarInfo('lon')
   vinfo_lat = self->VarInfo('lat')

   case !true of

      (vinfo_lon.n_dims eq 1) && (vinfo_lat.n_dims eq 1): begin
         lonlat2d = !false
         dim = [self->DimInfo(vinfo_lon.dim_names[0], /DIMSIZE), $
                self->DimInfo(vinfo_lat.dim_names[0], /DIMSIZE)]
      endcase

      (vinfo_lon.n_dims eq 2) && (vinfo_lat.n_dims eq 2): begin
         dim_names = vinfo_lon.dim_names
         if ~ array_equal(vinfo_lat.dim_names, dim_names) then $
            message, 'Longitude and latitude dimensions are mismatched'
         lonlat2d = !true
         dim = [self->DimInfo(dim_names[0], /DIMSIZE), $
                self->DimInfo(dim_names[1], /DIMSIZE)]
      endcase

      else: $
         message, 'Longitude and latitude dimensions are mismatched'

   endcase

   ;; Establish horizontal range. Default is to retrieve all data.

   if n_elements(i_range) eq 0 then i_range = [0,dim[0]-1]
   if n_elements(j_range) eq 0 then j_range = [0,dim[1]-1]

   ;; Interpret negative values in I_RANGE and J_RANGE as offsets
   ;; from the end of the grid.

   if i_range[0] lt 0 then i_range[0] += dim[0]
   if i_range[1] lt 0 then i_range[1] += dim[0]
   if j_range[0] lt 0 then j_range[0] += dim[1]
   if j_range[1] lt 0 then j_range[1] += dim[1]

   ;; Calculate parameters for retrieving & processing horizontal grid
   ;; data, which are defined on rho points.

   offset = [i_range[0],j_range[0]]
   count = [i_range[1]-i_range[0]+1,j_range[1]-j_range[0]+1]

   ;; Retrieve & interpolate horizontal position

   if lonlat2d then begin
      lon = self->VarGet('lon', OFFSET=offset, COUNT=count)
      lat = self->VarGet('lat', OFFSET=offset, COUNT=count)
   endif else begin
      lon = self->VarGet('lon', OFFSET=offset[0], COUNT=count[0])
      lat = self->VarGet('lat', OFFSET=offset[1], COUNT=count[1])
   endelse


   result = dictionary()

   result.lonlat2d = lonlat2d
   result.dim = dim
   result.i_range = i_range
   result.j_range = j_range
   result.lon = lon
   result.lat = lat

   return, result->ToStruct(/RECURSIVE, /NO_COPY)

end

; MGHromsRegrid::LocateXY
;
function MGHromsRegrid::LocateXY, x, y

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(lonlat) eq 0 then $
        lonlat = self->HasVar('lon_rho') && self->HasVar('lat_rho')

   if lonlat then begin
      xr = self->VarGet('lon_rho')
      yr = self->VarGet('lat_rho')
   endif else begin
      xr = self->VarGet('x_rho')
      yr = self->VarGet('y_rho')
   endelse

   return, mgh_locate2(xr, yr, XOUT=[x], YOUT=[y])

end

; MGHromsRegrid::TimeVarName
;
function MGHromsRegrid::TimeVarName, var

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

; MGHromsRegrid::VarDims

function MGHromsRegrid::VarDims, var

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Get dimension names as a string array.

   dim = self->VarDimNames(var)

   ;; Set up the output structure

   result = {horizontal: strarr(2), vertical: '', bed: '', time:''}

   ;; Load variable dimensions into output

   for d=0,n_elements(dim)-1 do begin
      if strmatch(dim[d], 'xi_*') || strmatch(dim[d], 'lon') then $
         result.horizontal[0] = dim[d]
      if strmatch(dim[d], 'eta_*') || strmatch(dim[d], 'lat') then $
         result.horizontal[1] = dim[d]
      if strmatch(dim[d], '*time*') then $
         result.time = dim[d]
      if strmatch(dim[d], 's_*') then $
         result.vertical = dim[d]
      if strmatch(dim[d], '*bed') then $
         result.bed = dim[d]
   endfor

   return, result

end

; MGHromsRegrid::VarDimNames
;
function MGHromsRegrid::VarDimNames, var, COUNT=count

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(var) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', var

   case !true of

      isa(var, 'STRING') && (var eq 'psi'): begin
         result = self->VarDimNames('psi(ubar,vbar,zeta)', COUNT=count)
      end

      isa(var, 'STRING') && strmatch(var, 'psi(*,*,*)'): begin
         pp = stregex(var, '(^psi\()(.+)(,)(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4,6]], ll[[2,4,6]])
         result = self->MGHncSequence::VarDimNames(vv[0], COUNT=count)
         result[0] = 'xi_psi'
         result[1] = 'eta_psi'
      end

      isa(var, 'STRING') && strmatch(var, 'abs(*,*)'): begin
         pp = stregex(var, '(^abs\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         result = self->MGHncSequence::VarDimNames(vv[0], COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && strmatch(var, 'harmonic_amp(*)'): begin
         pp = stregex(var, '(^harmonic_amp\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         result = self->MGHncSequence::VarDimNames(vv[0], COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && strmatch(var, 'harmonic_pha(*)'): begin
         pp = stregex(var, '(^harmonic_pha\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         result = self->MGHncSequence::VarDimNames(vv[0], COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && strmatch(var, 'tide_scalar_amp(*)'): begin
         pp = stregex(var, '(^tide_scalar_amp\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         result = self->MGHncSequence::VarDimNames(vv[0], COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && strmatch(var, 'tide_scalar_pha(*)'): begin
         pp = stregex(var, '(^tide_scalar_pha\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         result = self->MGHncSequence::VarDimNames(vv[0], COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_sma(*,*)'): begin
         pp = stregex(var, '(^tide_vector_sma\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         result = self->MGHncSequence::VarDimNames(vv[0]+'_cos', COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_ecc(*,*)'): begin
         pp = stregex(var, '(^tide_vector_ecc\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         result = self->MGHncSequence::VarDimNames(vv[0]+'_cos', COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_inc(*,*)'): begin
         pp = stregex(var, '(^tide_vector_inc\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         result = self->MGHncSequence::VarDimNames(vv[0]+'_cos', COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_pha(*,*)'): begin
         pp = stregex(var, '(^tide_vector_pha\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         result = self->MGHncSequence::VarDimNames(vv[0]+'_cos', COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && strmatch(var, 'component(*,*)'): begin
         pp = stregex(var, '(^component\()(.+)(,)(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4,6]], ll[[2,4,6]])
         result = self->MGHncSequence::VarDimNames(vv[0], COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && strmatch(var, 'uvbar@*'): begin
         result = self->MGHncSequence::VarDimNames('ubar', COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && strmatch(var, 'uv@*'): begin
         result = self->MGHncSequence::VarDimNames('u', COUNT=count)
         result[0] = 'xi_rho'
         result[1] = 'eta_rho'
      end

      isa(var, 'STRING') && (var eq 'ugs'): begin
         result = self->MGHncSequence::VarDimNames('zeta', COUNT=count)
         result[0] = 'xi_v'
         result[1] = 'eta_v'
      end

      isa(var, 'STRING') && (var eq 'vgs'): begin
         result = self->MGHncSequence::VarDimNames('zeta', COUNT=count)
         result[0] = 'xi_u'
         result[1] = 'eta_u'
      end

      isa(var, 'STRING') && strmatch(var, '*+*'): begin
         vars = strsplit(var, '+', /EXTRACT)
         result = self->VarDimNames(vars[0], COUNT=count)
      end

      isa(var, 'STRING'): begin
         result =self->MGHncSequence::VarDimNames(var, COUNT=count)
      end

      else: message, 'Invalid argument: var'

   endcase

   return, result

end

; MGHromsRegrid::VarGet
;
function MGHromsRegrid::VarGet, var, $
     AUTOSCALE=autoscale, COUNT=count, OFFSET=offset

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(var) eq 0 then $
      message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', var

   if n_elements(autoscale) eq 0 then autoscale = 1B

   case 1B of

      isa(var, 'STRING') && strmatch(var, 'abs(*,*)'): begin
         pp = stregex(var, '(^abs\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         data = self->VectorGet(vv, COUNT=count, OFFSET=offset)
         return, abs(temporary(data))
      end

      isa(var, 'STRING') && strmatch(var, 'harmonic_amp(*)'): begin
         pp = stregex(var, '(^harmonic_amp\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         data = self->VectorGet(vv, COUNT=count, OFFSET=offset)
         return, abs(temporary(data))
      end

      isa(var, 'STRING') && strmatch(var, 'harmonic_pha(*)'): begin
         pp = stregex(var, '(^harmonic_pha\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         data = self->VectorGet(vv, COUNT=count, OFFSET=offset)
         result = atan(temporary(data), /PHASE)/(2*!pi)
         result = (result + 1) mod 1
         return, result
      end

      isa(var, 'STRING') && strmatch(var, 'tide_scalar_amp(*)'): begin
         pp = stregex(var, '(^tide_scalar_amp\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         data = self->VectorGet(vv, COUNT=count, OFFSET=offset)
         return, abs(temporary(data))
      end

      isa(var, 'STRING') && strmatch(var, 'tide_scalar_pha(*)'): begin
         pp = stregex(var, '(^tide_scalar_pha\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         data = self->VectorGet(vv, COUNT=count, OFFSET=offset)
         return, (360+!radeg*atan(temporary(data), /PHASE)) mod 360
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_sma(*,*)'): begin
         pp = stregex(var, '(^tide_vector_sma\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         ;; Get u and v tidal coefficients. The complex plane represents (xi,eta) space.
         ;; This method returns all data on the same (RHO) grid automatically
         uv_cos = self->VectorGet(vv+'_cos', COUNT=count, OFFSET=offset, MASK_VALUE=0)
         uv_sin = self->VectorGet(vv+'_sin', COUNT=count, OFFSET=offset, MASK_VALUE=0)
         ;; No rotation to geographic coordiantes required for semi-major axis
         ;; Convert to u and v tidal coefficients as returned by MGH_TIDE_VECTOR_ANALYSIS. The
         ;; complex plane now represents temporal variation, I think.
         u_coeff = complex(real_part(uv_cos), real_part(uv_sin))
         v_coeff = complex(imaginary(uv_cos), imaginary(uv_sin))
         mgh_undefine, uv_cos, uv_sin
         ;; Convert to u and v amplitude and phase, as represented, for example, in the
         ;; NIWA EEZ tidal model output files.
         u_amp = abs(u_coeff)
         u_pha = (atan(u_coeff, /PHASE)*!radeg+360) mod 360
         v_amp = abs(v_coeff)
         v_pha = (atan(v_coeff, /PHASE)*!radeg+360) mod 360
         ;; Convert to tidal ellipse form
         ellipse = mgh_tide_ap2ep(u_amp, u_pha, v_amp, v_pha)
         mgh_undefine, u_amp, u_pha, v_amp, v_pha
         ;; Return the semi-major axis
         return, ellipse.sma
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_ecc(*,*)'): begin
         pp = stregex(var, '(^tide_vector_ecc\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         ;; Get u and v tidal coefficients. The complex plane represents (xi,eta) space.
         ;; This method returns all data on the same (RHO) grid automatically
         uv_cos = self->VectorGet(vv+'_cos', COUNT=count, OFFSET=offset, MASK_VALUE=0)
         uv_sin = self->VectorGet(vv+'_sin', COUNT=count, OFFSET=offset, MASK_VALUE=0)
         ;; No rotation to geographic coordinates required for eccentricity
         ;; Convert to u and v tidal coefficients as returned by MGH_TIDE_VECTOR_ANALYSIS. The
         ;; complex plane now represents temporal variation, I think.
         u_coeff = complex(real_part(uv_cos), real_part(uv_sin))
         v_coeff = complex(imaginary(uv_cos), imaginary(uv_sin))
         mgh_undefine, uv_cos, uv_sin
         ;; Convert to u and v amplitude and phase, as represented, for example, in the
         ;; NIWA EEZ tidal model output files.
         u_amp = abs(u_coeff)
         u_pha = (atan(u_coeff, /PHASE)*!radeg+360) mod 360
         v_amp = abs(v_coeff)
         v_pha = (atan(v_coeff, /PHASE)*!radeg+360) mod 360
         ;; Convert to tidal ellipse form
         ellipse = mgh_tide_ap2ep(u_amp, u_pha, v_amp, v_pha)
         mgh_undefine, u_amp, u_pha, v_amp, v_pha
         ;; Return the eccentricity
         return, ellipse.ecc
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_inc(*,*)'): begin
         pp = stregex(var, '(^tide_vector_inc\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         ;; Get u and v tidal coefficients. The complex plane represents (xi,eta) space.
         ;; This method returns all data on the same (RHO) grid automatically.
         uv_cos = self->VectorGet(vv+'_cos', COUNT=count, OFFSET=offset, MASK_VALUE=0)
         uv_sin = self->VectorGet(vv+'_sin', COUNT=count, OFFSET=offset, MASK_VALUE=0)
         ;; Rotate to geographic coordinates.
         if self->HasVar('angle') then begin
            angle = self->VarGet('angle', COUNT=count[0:1], OFFSET=offset[0:1])
            mgh_uv_rotate, uv_cos, angle
            mgh_uv_rotate, uv_sin, angle
         endif
         ;; Convert to u and v tidal coefficients as returned by MGH_TIDE_VECTOR_ANALYSIS. The
         ;; complex plane now represents temporal variation, I think.
         u_coeff = complex(real_part(uv_cos), real_part(uv_sin))
         v_coeff = complex(imaginary(uv_cos), imaginary(uv_sin))
         mgh_undefine, uv_cos, uv_sin
         ;; Convert to u and v amplitude and phase, as represented, for example, in the
         ;; NIWA EEZ tidal model output files.
         u_amp = abs(u_coeff)
         u_pha = (atan(u_coeff, /PHASE)*!radeg+360) mod 360
         v_amp = abs(v_coeff)
         v_pha = (atan(v_coeff, /PHASE)*!radeg+360) mod 360
         ;; Convert to tidal ellipse form
         ellipse = mgh_tide_ap2ep(u_amp, u_pha, v_amp, v_pha)
         mgh_undefine, u_amp, u_pha, v_amp, v_pha
         ;; Return the inclination
         return, ellipse.inc
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_pha(*,*)'): begin
         pp = stregex(var, '(^tide_vector_pha\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         ;; Get u and v tidal coefficients. The complex plane represents (xi,eta) space.
         ;; This method returns all data on the same (RHO) grid automatically.
         uv_cos = self->VectorGet(vv+'_cos', COUNT=count, OFFSET=offset, MASK_VALUE=0)
         uv_sin = self->VectorGet(vv+'_sin', COUNT=count, OFFSET=offset, MASK_VALUE=0)
         ;; Rotate to geographic coordinates.
         if self->HasVar('angle') then begin
            angle = self->VarGet('angle', COUNT=count[0:1], OFFSET=offset[0:1])
            mgh_uv_rotate, uv_cos, angle
            mgh_uv_rotate, uv_sin, angle
         endif
         ;; Convert to u and v tidal coefficients as returned by MGH_TIDE_VECTOR_ANALYSIS. The
         ;; complex plane now represents temporal variation, I think.
         u_coeff = complex(real_part(uv_cos), real_part(uv_sin))
         v_coeff = complex(imaginary(uv_cos), imaginary(uv_sin))
         mgh_undefine, uv_cos, uv_sin
         ;; Convert to u and v amplitude and phase, as represented, for example, in the
         ;; NIWA EEZ tidal model output files.
         u_amp = abs(u_coeff)
         u_pha = (atan(u_coeff, /PHASE)*!radeg+360) mod 360
         v_amp = abs(v_coeff)
         v_pha = (atan(v_coeff, /PHASE)*!radeg+360) mod 360
         ;; Convert to tidal ellipse form
         ellipse = mgh_tide_ap2ep(u_amp, u_pha, v_amp, v_pha)
         mgh_undefine, u_amp, u_pha, v_amp, v_pha
         ;; Return the inclination
         return, ellipse.pha
      end

      isa(var, 'STRING') && strmatch(var, 'component(*,*,*)'): begin
         pp = stregex(var, '(^component\()(.+)(,)(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4,6]], ll[[2,4,6]])
         ;; Velocity component direction, geographic convention, degrees.
         direction = float(vv[2])
         uv = self->VectorGet(vv[0:1], COUNT=count, OFFSET=offset, MASK_VALUE=0)
         if self->HasVar('angle') then begin
            angle = self->VarGet('angle', COUNT=count[0:1], OFFSET=offset[0:1])
            mgh_uv_rotate, uv, angle-!dtor*(90-direction)
         endif
         return, real_part(temporary(uv))
      end

      isa(var, 'STRING') && strmatch(var, 'uvbar@*'): begin
         ;; Velocity component direction, geographic convention, degrees.
         direction = float(strmid(var, 6))
         uv = self->VectorGet(['ubar','vbar'], COUNT=count, OFFSET=offset, MASK_VALUE=0)
         if self->HasVar('angle') then begin
            angle = self->VarGet('angle', COUNT=count[0:1], OFFSET=offset[0:1])
            mgh_uv_rotate, uv, angle-!dtor*(90-direction)
         endif
         return, real_part(temporary(uv))
      end

      isa(var, 'STRING') && strmatch(var, 'uv@*'): begin
         ;; Velocity component direction, geographic convention, degrees.
         direction = float(strmid(var, 3))
         uv = self->VectorGet(['u','v'], COUNT=count, OFFSET=offset, MASK_VALUE=0)
         if self->HasVar('angle') then begin
            angle = self->VarGet('angle', COUNT=count[0:1], OFFSET=offset[0:1])
            mgh_uv_rotate, uv, angle-!dtor*(90-direction)
         endif
         return, real_part(temporary(uv))
      end

      isa(var, 'STRING') && (var eq 'ugs'): begin
         ;; The following works only if zeta is 2-dimensional--nightmare!
         ;; See the code for "psi" for the right way to do it!!
         dim_rho = self->DimRho()
         my_offset = n_elements(offset) gt 0 ? offset : [0,0,0]
         my_count = n_elements(count) gt 0 ? count : [dim_rho[0:1]+[0,-1],0]
         pn = self->MGHncSequence::VarGet('pn', COUNT=mY_count[0:1]+[0,1], OFFSET=my_offset[0:1])
         f = self->MGHncSequence::VarGet('f', COUNT=my_count[0:1]+[0,1], OFFSET=my_offset[0:1])
         if self->MGHncSequence::HasVar('mask_rho') then $
            mask = self->MGHncSequence::VarGet('mask_rho', COUNT=my_count[0:1]+[0,1], OFFSET=my_offset[0:1])
         zeta = self->MGHncSequence::VarGet('zeta', AUTOSCALE=autoscale, COUNT=my_count+[0,1,0], OFFSET=my_offset)
         result = mgh_roms_gsvel(temporary(zeta), temporary(pn), temporary(f), DIRECTION=0, MASK=mask)
      end

      isa(var, 'STRING') && (var eq 'vgs'): begin
         ;; The following works only if zeta is 2-dimensional--nightmare!
         ;; See the code for "psi" for the right way to do it!!
         dim_rho = self->DimRho()
         my_offset = n_elements(offset) gt 0 ? offset : [0,0,0]
         my_count = n_elements(count) gt 0 ? count : [dim_rho[0:1]+[-1,0],0]
         pm = self->MGHncSequence::VarGet('pm', COUNT=my_count[0:1]+[1,0], OFFSET=my_offset[0:1])
         f = self->MGHncSequence::VarGet('f', COUNT=my_count[0:1]+[1,0], OFFSET=my_offset[0:1])
         if self->MGHncSequence::HasVar('mask_rho') then begin
            mask = self->MGHncSequence::VarGet('mask_rho', COUNT=my_count[0:1]+[1,0], OFFSET=my_offset[0:1])
         endif
         zeta = self->MGHncSequence::VarGet('zeta', AUTOSCALE=autoscale, COUNT=my_count+[1,0,0], OFFSET=my_offset)
         result = mgh_roms_gsvel(temporary(zeta), temporary(pm), temporary(f), DIRECTION=1, MASK=mask)
      end

      isa(var, 'STRING') && strmatch(var, '*time*'): begin
         result = self->MGHncSequence::VarGet(var, AUTOSCALE=0, COUNT=count, OFFSET=offset)
         ;; If the AUTOSCALE keyword is set, apply the scale factor specified in the time
         ;; variable's "units" attribute. Ignore the offset as this is problematic.
         if keyword_set(autoscale) then begin
            time_units = {scale: 1D/(24D*3600D), offset: 0D}
            if self->MGHncSequence::HasAtt(var, 'units') then begin
               tu = self->MGHncSequence::AttGet(var, 'units')
               time_units = mgh_dt_units(tu)
            endif
            result *= time_units.scale
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

      isa(var, 'STRING'): begin
         result = self->MGHncSequence::VarGet(var, AUTOSCALE=autoscale, COUNT=count, OFFSET=offset)
      end

      else: message, 'Invalid argument: var'

   endcase

   return, result

end

pro MGHromsRegrid::VarInfo, var, $
   ALL=all, ATT_NAMES=att_names, DATATYPE=datatype, DIM_NAMES=dim_names, $
   DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_ATTS=n_atts, N_DIMS=n_dims

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case !true of

      isa(var, 'STRING') && strmatch(var, 'abs(*,*)'): begin
         pp = stregex(var, '(^abs\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, vv[0], $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
         dim_names[0] = 'xi_rho'
         dim_names[1] = 'eta_rho'
         dimensions[0] = self->MGHncSequence::DimInfo('xi_rho', /DIMSIZE)
         dimensions[1] = self->MGHncSequence::DimInfo('eta_rho', /DIMSIZE)
      end

      isa(var, 'STRING') && strmatch(var, 'harmonic_amp(*)'): begin
         pp = stregex(var, '(^harmonic_amp\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, vv[0], $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
         dim_names[0] = 'xi_rho'
         dim_names[1] = 'eta_rho'
         dimensions[0] = self->MGHncSequence::DimInfo('xi_rho', /DIMSIZE)
         dimensions[1] = self->MGHncSequence::DimInfo('eta_rho', /DIMSIZE)
      end

      isa(var, 'STRING') && strmatch(var, 'harmonic_pha(*)'): begin
         pp = stregex(var, '(^harmonic_pha\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, vv[0], $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
         dim_names[0] = 'xi_rho'
         dim_names[1] = 'eta_rho'
         dimensions[0] = self->MGHncSequence::DimInfo('xi_rho', /DIMSIZE)
         dimensions[1] = self->MGHncSequence::DimInfo('eta_rho', /DIMSIZE)
      end

      isa(var, 'STRING') && strmatch(var, 'tide_scalar_amp(*)'): begin
         pp = stregex(var, '(^tide_scalar_amp\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, vv[0], $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
         dim_names[0] = 'xi_rho'
         dim_names[1] = 'eta_rho'
         dimensions[0] = self->MGHncSequence::DimInfo('xi_rho', /DIMSIZE)
         dimensions[1] = self->MGHncSequence::DimInfo('eta_rho', /DIMSIZE)
      end

      isa(var, 'STRING') && strmatch(var, 'tide_scalar_pha(*)'): begin
         pp = stregex(var, '(^tide_scalar_pha\()(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[2], ll[2])+['_cos','_sin']
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, vv[0], $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
         dim_names[0] = 'xi_rho'
         dim_names[1] = 'eta_rho'
         dimensions[0] = self->MGHncSequence::DimInfo('xi_rho', /DIMSIZE)
         dimensions[1] = self->MGHncSequence::DimInfo('eta_rho', /DIMSIZE)
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_sma(*,*)'): begin
         pp = stregex(var, '(^tide_vector_sma\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, vv[0], $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
         dim_names[0] = 'xi_rho'
         dim_names[1] = 'eta_rho'
         dimensions[0] = self->MGHncSequence::DimInfo('xi_rho', /DIMSIZE)
         dimensions[1] = self->MGHncSequence::DimInfo('eta_rho', /DIMSIZE)
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_ecc(*,*)'): begin
         pp = stregex(var, '(^tide_vector_ecc\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, vv[0], $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
         dim_names[0] = 'xi_rho'
         dim_names[1] = 'eta_rho'
         dimensions[0] = self->MGHncSequence::DimInfo('xi_rho', /DIMSIZE)
         dimensions[1] = self->MGHncSequence::DimInfo('eta_rho', /DIMSIZE)
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_inc(*,*)'): begin
         pp = stregex(var, '(^tide_vector_inc\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, vv[0], $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
         dim_names[0] = 'xi_rho'
         dim_names[1] = 'eta_rho'
         dimensions[0] = self->MGHncSequence::DimInfo('xi_rho', /DIMSIZE)
         dimensions[1] = self->MGHncSequence::DimInfo('eta_rho', /DIMSIZE)
      end

      isa(var, 'STRING') && strmatch(var, 'tide_vector_pha(*,*)'): begin
         pp = stregex(var, '(^tide_vector_pha\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, vv[0], $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
         dim_names[0] = 'xi_rho'
         dim_names[1] = 'eta_rho'
         dimensions[0] = self->MGHncSequence::DimInfo('xi_rho', /DIMSIZE)
         dimensions[1] = self->MGHncSequence::DimInfo('eta_rho', /DIMSIZE)
      end

      else: begin
         self->MGHncSequence::VarInfo, var, $
            ATT_NAMES=att_names, DATATYPE=datatype, DIM_NAMES=dim_names, $
            DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_ATTS=n_atts, N_DIMS=n_dims
      end

   endcase

   if arg_present(all) then begin
      all = {att_names: att_names, datatype: datatype, dim_names: dim_names, $
         dimensions: dimensions, fill_value: fill_value, n_atts: n_atts, n_dims: n_dims}
   endif

end
; MGHromsRegrid::VectorGet
;
function MGHromsRegrid::VectorGet, var, $
     COUNT=count, OFFSET=offset, MASK_VALUE=mask_value

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(var) eq 0 then var = ['u','v']

   if n_elements(var) ne 2 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_WRGNUMELEM', var

   if size(var, /TYPE) ne 7 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_WRONGTYPE', var

   if n_elements(mask_value) eq 0 then mask_value = 0

   dim_u = self->VarDims(var[0])
   dim_v = self->VarDims(var[1])

   case !true of
      array_equal(dim_u.horizontal, ['xi_rho','eta_rho']) && $
         array_equal(dim_v.horizontal, ['xi_rho','eta_rho']): begin
         grid_type = 'rho'
      end
      array_equal(dim_u.horizontal, ['xi_u','eta_u']) && $
         array_equal(dim_v.horizontal, ['xi_v','eta_v']): begin
         grid_type = 'uv'
      end
      else: begin
         message, 'Unrecognised horizontal grid type'
      end
   endcase

   fmt = '(%"Variables %s & %s have mismatched netCDF dimensions")'
   if dim_u.vertical ne dim_v.vertical then $
        message, string(FORMAT=fmt, var)
   if dim_u.time ne dim_v.time then $
        message, string(FORMAT=fmt, var)

   case grid_type of

      'rho': begin

         u = self->VarGet(var[0], COUNT=count, OFFSET=offset, /AUTOSCALE)
         v = self->VarGet(var[1], COUNT=count, OFFSET=offset, /AUTOSCALE)

         l_miss = where(~ finite(u), n_miss)
         if n_miss gt 0 then u[l_miss] = mask_value

         l_miss = where(~ finite(v), n_miss)
         if n_miss gt 0 then v[l_miss] = mask_value

      end

      'uv': begin

         dim = [self->DimInfo('xi_rho', /DIMSIZE),self->DimInfo('eta_rho', /DIMSIZE)]

         ;; Set up default values for offset and count. The offset is relative to
         ;; the full RHO grid

         n_var_dim = 2
         default_offset = [0,0]
         default_count = dim
         if strlen(dim_u.vertical) gt 0 then begin
            n_var_dim += 1
            default_offset = [default_offset,0]
            default_count = [default_count,self->DimInfo(dim_u.vertical, /DIMSIZE)]
         endif
         if strlen(dim_u.time) gt 0 then begin
            n_var_dim += 1
            default_offset = [default_offset,0]
            default_count = [default_count,self->DimInfo(dim_u.time, /DIMSIZE)]
         endif

         mgh_undefine, dim_u, dim_v

         my_offset = n_elements(offset) gt 0 ? offset : default_offset
         my_count = n_elements(count) gt 0 ? count : default_count

         if my_count[0] eq 0 then my_count[0] = dim[0]-my_offset[0]
         if my_count[1] eq 0 then my_count[1] = dim[1]-my_offset[1]

         ;; Vector data are calculated only on the rho grid interior points but
         ;; are padded as necessary with NaNs to occupy the full rho grid.

         pad_west = my_offset[0] lt 1
         pad_east = (my_offset[0]+my_count[0]-1) gt (dim[0]-2)
         pad_south = my_offset[1] lt 1
         pad_north = (my_offset[1]+my_count[1]-1) gt (dim[1]-2)

         ;; Get u data, fill mask values and interpolate to RHO grid

         offset_u = my_offset
         offset_u[0] = offset_u[0]-1+keyword_set(pad_west)
         count_u = my_count
         count_u[0] = count_u[0]+1-keyword_set(pad_west)-keyword_set(pad_east)

         u = self->MGHncSequence::VarGet(var[0], /AUTOSCALE, OFFSET=offset_u, COUNT=count_u)

         l_miss = where(~ finite(u), n_miss)
         if n_miss gt 0 then u[l_miss] = mask_value

         if keyword_set(pad_west) || keyword_set(pad_east) then begin
            dim_u = size(u, /DIMENSIONS)
            dim_u[0] = dim_u[0] + keyword_set(pad_west) + keyword_set(pad_east)
            uu = replicate(!values.f_nan, dim_u)
            case n_var_dim of
               2: uu[keyword_set(pad_west),0] = temporary(u)
               3: uu[keyword_set(pad_west),0,0] = temporary(u)
               4: uu[keyword_set(pad_west),0,0,0] = temporary(u)
            endcase
            u = temporary(uu)
         endif

         delta = replicate(0, n_var_dim)
         delta[0] = -1
         u = mgh_stagger(temporary(u), DELTA=delta)

         ;; Get v data, fill mask values and interpolate to RHO grid

         offset_v = my_offset
         offset_v[1] = offset_v[1]-1+keyword_set(pad_south)
         count_v = my_count
         count_v[1] = count_v[1]+1-keyword_set(pad_south)-keyword_set(pad_north)

         v = self->MGHncSequence::VarGet(var[1], /AUTOSCALE, OFFSET=offset_v, COUNT=count_v)

         l_miss = where(~ finite(v), n_miss)
         if n_miss gt 0 then v[l_miss] = mask_value

         if keyword_set(pad_south) || keyword_set(pad_north) then begin
            dim_v = size(v, /DIMENSIONS)
            dim_v[1] = dim_v[1] + keyword_set(pad_south) + keyword_set(pad_north)
            vv = replicate(!values.f_nan, dim_v)
            case n_var_dim of
               2: vv[0,keyword_set(pad_south),0] = temporary(v)
               3: vv[0,keyword_set(pad_south),0,0] = temporary(v)
               4: vv[0,keyword_set(pad_south),0,0,0] = temporary(v)
            endcase
            v = temporary(vv)
         endif

         delta = replicate(0, n_var_dim)
         delta[1] = -1
         v = mgh_stagger(temporary(v), DELTA=delta)

      end

   endcase

   ;; Return data as complex number

   return, complex(temporary(u), temporary(v))

end

pro MGHromsRegrid__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, {MGHromsRegrid, inherits MGHncSequence}

end
