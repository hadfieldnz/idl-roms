;+
; CLASS NAME:
;   MGHromsHistory
;
; PURPOSE:
;   This class wraps a series of ROMS history (or similar) files.
;
; OBJECT CREATION SEQUENCE
;   obj = obj_new('MGHromsHistory', Files)
;
; SLICES:
;   Several of the methods below are used to support extraction of
;   data on spatial slices through the ROMS grid. The slice types are:
;
;     C-slice
;       Vertical along the grid lines in the direction of increasing
;       xi (DIRECTION=0) or eta (DIRECTION=1).
;
;     H-slice
;       Horizontal, i.e. constant depth, constant s, constant sigma,
;       constant height above the bottom, or on a sediment layer.
;
;     P-slice
;       A vertical slice along the grid's normal velocity surfaces,
;       defined by a series of vertices on the PSI grid. P-slices
;       are intended to allow the exact calculation of volume and
;       tracer fluxes.
;
;     R-slice
;       A 2D horizontal slice on a rectilinear lon-lat grid. Not yet
;       implemented (2017-05-03).
;
;     X-slice
;       A vertical slice along a straight or curved line (in x,y or lon,lat
;       space). The line is defined by at least two points and can be
;       (piecewise) linear or curved, with intermediate points determined
;       by a spline function.
;
;   The vertical slices (Cslice, Pslice and Xslice) naturally reduce in the
;   vertical to horizontal (constant z, constant s, or depth-invariant)
;   transects along the same path. These are called C-transect, P-transect
;   and X-transect, abbreviated Ctran, Ptran and Xtran.
;
;   Hslices and Cslices are specified relative to the grid associated
;   with a particular variable. So in constructing an Hslice or Cslice
;   one must specify the associated variable (or the variable dimensions.
;
;   Pslices are defined relative to the PSI grid and Xslices are defined in
;   physical space, so no associated variable is required.
;
;   P-slices can be thought of as generalisations of the slice geometry implied
;   in the GetTransportSlice method, though in the latter case we (perversely, I think)
;   define the slice on the RHO grid and then draw it through the normal-velocity
;   faces on the south or west sides of the cells.
;
; PATCHES:
;   A Patch is a polygonal area over which integrals (and maybe other statistics)
;   can be calculated. Its vertices are defined in physical space.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 1999-08:
;     Written.
;   Mark Hadfield, 2003-10:
;     Deactivate code that keeps track of a grid file, from which
;     data can be retrieved if not found in the history file.  -
;   Mark Hadfield, 2004-07:
;     - Fixed bug: depths returned by CsliceZ method for w grids
;       were incorrect for old-style ROMS history files, in which the
;       lowest w level is not written to the file.
;   Mark Hadfield, 2005-10:
;     - Added HsliceMean method, returning time-averaged data for an
;       Hslice.
;   Mark Hadfield, 2009-03:
;     - Fixed a long-standing bug in the GetTransportBox method: slices were
;       mislocated by 1, leading to small errors and imbalances in transport
;       integrated around the box.
;   Mark Hadfield, 2010-10:
;     - Variable type checking now done with new ISA function.
;   Mark Hadfield, 2010-11:
;     - Added support for synthetic variable bstr (scalar bottom stress).
;   Mark Hadfield, 2011-03:
;     - Calls to ISA function removed for compatibility with IDL 7.1.
;   Mark Hadfield, 2011-03:
;     - Begin adding supporting for sediment (the bed dimension) and
;       time-varying bathymetry.
;   Mark Hadfield, 2011-07:
;     - Removed references to grid_destroy: unnecessary with automatic
;       garbage collection.
;     - Added support for synthetic variables sstr (scalar surface stress)
;       and Sbot (bottom orbital speed).
;   Mark Hadfield, 2012-02:
;     - Adapted for ROMS-AGRIF files.
;   Mark Hadfield, 2012-08:
;     - Fixed up some more issues related to ROMS-AGRIF files.
;   Mark Hadfield, 2012-09:
;     - Default bed layer for Hslices is 0.
;     - Made the test for masked values in the HsliceData method more
;       tolerant of inexact 0 & 1 values introduced by packing.
;     - S-coordinate values read from history files are now clipped
;       to [-1,0] (another packing fix).
;   Mark Hadfield, 2012-12:
;     - Fix bug in code to calculate Dbbl: depth(h) not autoscaled.
;   Mark Hadfield, 2013-05:
;     - Cslice methods overhauled: a Cslice is now on the rho grid by
;       definition; the SLICE keyowrd has been renamed INDEX.
;     - GetTransportCslice method renamed GetTransportSlice.
;   Mark Hadfield, 2014-06:
;     - The end-points of an Xslice are now specified in (x,y) or
;       (lon,lat) space. The keywords END0 and END1 have been
;       replaced by VERTX and VERTY, with a view to allowing
;       the number of vertices to exceed 2 in the future.
;   Mark Hadfield, 2015-06:
;     - The GetZgrid method is deprecated in favour of the stand-alone function
;       MGH_ROMS_ZGRID.
;   Mark Hadfield, 2015-11:
;     - Adapted the HasVar and VarInfo methods to handle synthetic variables.
;     - Overhauled synthetic variables. The synthetic variable names now typically
;       encode the native variables from which they are calculated. A large
;       number of synthetic variables (eg. spd, sbar, sstr, Sbot) have been replaced
;       by a single pattern, "abs(*,*)".
;     - Stripped out the code supporting the grid_file property, which allowed
;       variables to be found in a secondary grid file if they were not in the main
;       file. This code was never completely implemented in any case.
;     - Generalised the VectorGet method to handle variables on the rho grid as well
;       as the u & v grid.
;   Mark Hadfield, 2015-12:
;     - The default XI_RANGE and ETA_RANgE for an HsliceGrid are now set to retrieve
;       all data.
;     - Fixed the VectorGet method: the padding required to handle variables on the
;       u & v grid. now works even for single-column or wingle-row retrievals.
;   Mark Hadfield, 2016-03:
;     - HsliceGrid now builds up its output data in a dictionary before conversion to
;       a structure for output. (I tried outputting the dictionary, but this leads to
;       a performance penalty downstream.)
;     - The HsliceData code for the case where DEPTHS is set has been reworked to
;       avoid "illegal operand" messages.
;   Mark Hadfield, 2016-04:
;     - The XsliceGrid method now supports "spline" and "piecewise" types in addition
;       to the previous straight-line/great-circle type, now called "linear". The input
;       keywords have been changed.
;     - The XsliceDefGrid method has been omitted, as it got too complicated with the
;       different keywords.
;   Mark Hadfield, 2016-09:
;     - All occurrences of the constant i (unit imaginary number) ahs been replaced with !const.i
;   Mark Hadfield, 2019-04:
;     - Deleted the Lslice methods, which were introduced some time ago and only
;       partially implemented and documented. These are now supported by the new
;       MGHromsRegrid class.
;-
;+
; METHOD NAME:
;   MGHromsHistory::CsliceGrid
;
; CATEGORY:
;   Ocean models.
;   Object graphics
;
; CALLING SEQUENCE:
;   slice = history->CsliceGrid(arg)
;
; POSITIONAL PARAMETERS:
;   arg
;     This argument can be either a string specifying the name of a
;     variable for which we want grid information or a structure with
;     tags specifying horizontal, vertical and time dimensions.
;
; KEYWORD PARAMETERS:
;   ALONG_RANGE
;     Range of points to be retrieved in the along-slice
;     direction. Default is to retrieve all physically meaningful data
;     in the domain
;
;   DIRECTION
;     Set this to 0 for a slice in the xi direction, 1 for a slice in
;     the eta direction. Default is 0.
;
;   LONLAT
;     This keyword specifies whether to retrieve horizontal position
;     in (lon, lat) or (x, y) form. Default is to retrieve (lon, lat)
;     data (LONLAT=1) if possible, otherwise retrieve (x, y)
;     (LONLAT=0).
;
;   SLICE
;     Position (grid index) of slice in across-slice direction.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2000-12:
;     Written.
;   Mark Hadfield, 2001-08:
;     Overhauled, removing some confusion about location and interpolation.
;-
;+
; METHOD NAME:
;   MGHromsHistory::HsliceData
;
; PURPOSE:
;   This class retrieves a slice through a ROMS 2D or 3D variable on a
;   set of Hslices (bed layer, constant s, constant z or
;   constant-sigma surfaces).
;
; CATEGORY:
;   Ocean models.
;   Object graphics
;
; CALLING SEQUENCE:
;   slice = history->HsliceData(varname)
;
; INPUTS:
;   history
;     A reference to a ROMS history sequence object.
;
;   variable
;     The name of a 2-D or 3-D variable in the netCDF file.
;
; KEYWORD PARAMETERS:
;   DEPTHS
;     Set this keyword to a numeric vector to specify the
;     depths of z surfaces on which data are to be plotted. This
;     keyword should be specified only for variables having a depth
;     coordinate and it cannot be used together with LEVELS or SIGMAS.
;
;   LAYERS
;     Set this keyword to an integer vector to specify the bed layers to
;     be plotted.  This keyword should be specified only for variables
;     having a bed-layer dimension.
;
;   LEVELS
;     Set this keyword to an integer vector to specify the
;     s-coordinate levels to be plotted.  This keyword should be
;     specified only for variables having a depth coordinate and it
;     cannot be used together with DEPTH or SIGMA.
;
;   SIGMAS
;     Set this keyword to a numeric vector to specify the
;     sigma values of constant-sigma surfaces on which data are to be
;     plotted. This keyword should be specified only for variables having
;     a depth coordinate and it cannot be used together with DEPTH or LEVEL.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2000-10:
;     Written based on code in MGH_ROMS_HSLICE_MOVIE
;   Mark Hadfield, 2011-05:
;     Added support for sediment bed layers and dynamic bathymetry.
;-
;+
; METHOD NAME:
;   MGHromsHistory::HsliceGrid
;
; CATEGORY:
;   Ocean models.
;   Object graphics
;
; CALLING SEQUENCE:
;   slice = history->HsliceGrid(arg)
;
; POSITIONAL PARAMETERS:
;   arg
;     This argument can be either a string specifying the name of a
;     variable for which we want grid information or a structure with
;     tags specifying horizontal, vertical and time dimensions. If it
;     is undefined then grid information is returned for variable "h".
;
; KEYWORD PARAMETERS:
;   ETA_RANGE
;     Range of points to be retrieved in the eta direction. Default
;     depends on variable dimensions.
;
;   LONLAT
;     This keyword specifies whether to retrieve horizontal position
;     in (lon, lat) or (x, y) form. Default is to retrieve (lon, lat)
;     data (LONLAT=1) if possible, otherwise retrieve (x, y) (LONLAT=0).
;
;   XI_RANGE
;     Range of points to be retrieved in the xi direction. Default
;     depends on variable dimensions.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2000-10:
;     Written based on code in MGH_ROMS_HSLICE_MOVIE
;   Mark Hadfield, 2011-05:
;     Added support for sediment bed layers and dynamic bathymetry.
;-
;+
; METHOD NAME:
;   MGHromsHistory::VarGet
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
;   MGHromsHistory::VarDims
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
function MGHromsHistory::Init, files, _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   return, self->MGHncSequence::Init(files, _STRICT_EXTRA=_extra)

end

pro MGHromsHistory::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGHncSequence::Cleanup

end

; ** C-slice and C-transect methods *******************************************
;
; MGHromsHistory::CsliceData
;
function MGHromsHistory::CsliceData, variable, $
     ALONG_RANGE=along_range, DIRECTION=direction, $
     GRID=grid, INDEX=index, LONLAT=lonlat, RECORD=record

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Get dimensions of variable

   dims = self->VarDims(variable)

   ;; Get grid info if necessary

   self->CsliceDefGrid, $
      ALONG_RANGE=along_range, DIRECTION=direction, $
      GRID=grid, INDEX=index, LONLAT=lonlat

   ;; Some handy constants

   has_vert = strlen(dims.vertical) gt 0
   has_time = strlen(dims.time) gt 0

   if ~ has_vert then $
      message, 'Variable has no vertical dimension'

   ;; Set defaults

   if n_elements(mask_value) eq 0 then mask_value = !values.f_nan

   ;; Check consistency of function arguments with variable dimensions.

   if has_time then begin
      if n_elements(record) eq 0 then record = 0
      if record lt 0 then $
         record = self->DimInfo(dims.time, /DIMSIZE) + record
   endif else begin
      if n_elements(record) gt 0 then begin
         message, 'The RECORD keyword is not required or allowed when ' + $
            'the variable '+variable+' has no time dimension'
      endif
   endelse

   ;; Get s-coordinate data

   scoord = self->GetScoord(dims.vertical)

   ;; Constants & abbreviations:

   ara0 = grid.along_range[0]
   ara1 = grid.along_range[1]
   aran = ara1-ara0+1

   ;; Specify parameters for getting data...

   ;; ...horizontal dimension

   case grid.direction of
      0: begin
         offset = [ara0,grid.index]
         count = [aran,1]
      end
      1: begin
         offset = [grid.index,ara0]
         count = [1,aran]
      end
   endcase

   delta = [0,0]

   case strjoin(dims.horizontal, ' ') of
      'xi_rho eta_rho':
      'xi_u eta_u': begin
         offset[0] -= 1
         count[0] += 1
         delta[0] -= 1
      end
      'xi_v eta_v': begin
         offset[1] -= 1
         count[1] += 1
         delta[1] -= 1
      end
      'xi_psi eta_psi': begin
         offset -= 1
         count += 1
         delta -= 1
      end
   endcase

   ;; ...vertical dimension

   offset = [offset,0]
   count = [count,0]
   delta = [delta,0]

   ;; ...time dimension

   if has_time then begin
      offset = [offset,record]
      count = [count,1]
      delta = [delta,0]
   endif

   ;; Get data, re-grid if necessary, and return

   result = self->VarGet(variable, OFFSET=offset, COUNT=count, /AUTOSCALE)

   return, reform(mgh_stagger(result, DELTA=delta))

end

; MGHromsHistory::CsliceDefGrid
;
pro MGHromsHistory::CsliceDefGrid, $
     DIRECTION=direction, INDEX=index, ALONG_RANGE=along_range, $
     GRID=grid, LONLAT=lonlat

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; If grid information is supplied, then check it is consistent with the
   ;; current function arguments; otherwise get the grid information.

   if n_elements(grid) gt 0 then begin

      if n_elements(along_range) gt 0 && ~ array_equal(along_range, grid.along_range) then $
         message, 'ALONG_RANGE does not match grid data'
      if n_elements(direction) gt 0 && ~ array_equal(direction, grid.direction) then $
         message, 'DIRECTION does not match grid data'
      if n_elements(index) gt 0 && ~ array_equal(index, grid.index) then $
         message, 'INDEX does not match grid data'
      if n_elements(lonlat) gt 0  && ~ array_equal(lonlat, grid.lonlat) then $
         message, 'LONLAT does not match grid data'

   endif else begin

      grid = self->CsliceGrid(DIRECTION=direction, INDEX=index, ALONG_RANGE=along_range, LONLAT=lonlat)

   endelse

end

; MGHromsHistory::CsliceGrid
;
function MGHromsHistory::CsliceGrid, $
     ALONG_RANGE=along_range, DIRECTION=direction, $
     INDEX=index, LONLAT=lonlat

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  if n_elements(direction) eq 0 then direction = 0

  if n_elements(lonlat) eq 0 then $
       lonlat = self->HasVar('lon_rho') && self->HasVar('lat_rho')

  ;; Get dimensions for RHO points

  dim_rho = self->DimRho()

  ;; Some useful switches

  has_mask = self->HasVar('mask_rho')

  ;; Establish horizontal range. Default is to include all interior RHO
  ;; points on a slice through the centre of the domain.

  if n_elements(index) eq 0 then $
       index = round(0.5*(dim_rho[1-direction]-2))

  if n_elements(along_range) eq 0 then along_range = [1,-2]

  ;; Interpret negative values in INDEX and ALONG_RANGE as offsets
  ;; from the end of the grid.

  if index lt 0 then $
       index += dim_rho[1-direction]
  if along_range[0] lt 0 then $
       along_range[0] += dim_rho[direction]
  if along_range[1] lt 0 then $
       along_range[1] += dim_rho[direction]

  n_along = along_range[1]-along_range[0]+1

  ;; Check values are within bounds

  ivalid = [0,dim_rho[1-direction]-1]
  avalid = [0,dim_rho[direction]-1]

  if index lt ivalid[0] then message, 'Slice is out of bounds'
  if index gt ivalid[1] then message, 'Slice is out of bounds'
  if along_range[0] lt avalid[0] then message, 'Slice is out of bounds'
  if along_range[1] gt avalid[1] then message, 'Slice is out of bounds'

  mgh_undefine, svalid, avalid

  ;; Create arrays of xi & eta location

  case direction of
    0: begin
      xi  = along_range[0]+lindgen(n_along)
      eta = replicate(index, n_along)
    end
    1: begin
      xi  = replicate(index, n_along)
      eta = along_range[0]+lindgen(n_along)
    end
  endcase

  ;; Specify parameters for retrieving rho grid data on the
  ;; slice

  case direction of
    0: begin
      count = [n_along,1]
      offset = [along_range[0],index]
    end
    1: begin
      count = [1,n_along]
      offset = [index,along_range[0]]
    end
  endcase

  ;; Retrieve inverse grid spacing data and calculate arc (along-slice)
  ;; distance

  case direction of
    0: ppvar = 'pm'
    1: ppvar = 'pn'
  endcase

  pp = mgh_stagger(reform(self->VarGet(ppvar, COUNT=count, OFFSET=offset)), DELTA=-1)

  ;; Calculate arc distance from origin

  arc = fltarr(n_along)
  for i=1,n_along-1 do arc[i] = arc[i-1] + 1/pp[i-1]

  ;; Retrieve horizontal position

  case lonlat of
    0: begin
      xvar = 'x_rho'
      yvar = 'y_rho'
    end
    1: begin
      xvar = 'lon_rho'
      yvar = 'lat_rho'
    end
  endcase

  x = reform(self->VarGet(xvar, OFFSET=offset, COUNT=count))
  y = reform(self->VarGet(yvar, OFFSET=offset, COUNT=count))

  ;; Retrieve bathymetry and mask

  h = 0
  h = reform(self->VarGet('h', OFFSET=offset, COUNT=count))

  if self->HasVar('mask_rho') then begin
    mask = reform(self->VarGet('mask_rho', OFFSET=offset, COUNT=count))
  endif else begin
    mask = replicate(1, n_along)
  endelse

  ;; Return result structure

  return, {lonlat: lonlat, direction: direction, $
           index: index, along_range: along_range, $
           x: x, y: y, arc: arc, xi: xi, eta: eta, $
           h: h, mask: mask}
end

; MGHromsHistory::CsliceZ
;
function MGHromsHistory::CsliceZ, variable, $
     ALONG_RANGE=along_range, BATH=bath, DIRECTION=direction, $
     GRID=grid, INDEX=index, LONLAT=lonlat, ZETA=zeta

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  ;; Get grid info if necessary.

  self->CsliceDefGrid, $
       DIRECTION=direction, INDEX=index, ALONG_RANGE=along_range, $
       GRID=grid, LONLAT=lonlat

  ;; Get the variable's vertical dimension

  dims = self->VarDims(variable)
  if strlen(dims.vertical) eq 0 then $
       message, 'Variable has no vertical dimension'

  ;; Get s-coordinate data

  scoord = self->GetScoord(dims.vertical)

  ;; Handy constants

  n_arc = n_elements(grid.arc)
  n_sc = n_elements(scoord.s)

  ;; Set defaults

  if n_elements(bath) eq 0 then bath = grid.h

  if n_elements(zeta) eq 0 then zeta = fltarr(n_arc)

  ;; Create result array & load data into it one column at a time.

  result = fltarr(n_sc, n_arc)

  cs = mgh_roms_s_to_cs(scoord.s, $
                        THETA_S=scoord.theta_s, THETA_B=scoord.theta_b, $
                        VSTRETCH=scoord.vstretch)

  for i=0,n_arc-1 do $
       result[0,i] = mgh_roms_s_to_z(scoord.s, bath[i], $
                                     ZETA=zeta[i], HC=scoord.hc, CS=cs, $
                                     VTRANSFORM=scoord.vtransform)

  return, transpose(result)

end

; MGHromsHistory::CtranData
;
function MGHromsHistory::CtranData, variable, $
     ALONG_RANGE=along_range, DIRECTION=direction, $
     DEPTHS=depths, LEVELS=levels, SIGMAS=sigmas, $
     GRID=grid, INDEX=index, LONLAT=lonlat, RECORD=record, $
     USE_BATH=use_bath, USE_ZETA=use_zeta

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  ;; Set defaults

  if n_elements(use_bath) eq 0 then use_bath = 0
  if n_elements(use_zeta) eq 0 then use_zeta = 0

  ;; Get dimensions of variable

  dims = self->VarDims(variable)

  ;; Get grid info if necessary

  self->CsliceDefGrid, $
       ALONG_RANGE=along_range, DIRECTION=direction, $
       GRID=grid, INDEX=index, LONLAT=lonlat

  ;; Some handy switches

  has_vert = strlen(dims.vertical) gt 0
  has_time = strlen(dims.time) gt 0

  ;; Check consistency of function arguments with variable dimensions.

  if has_vert then begin
    n_key = (n_elements(depths) gt 0) + (n_elements(levels) gt 0) + $
            (n_elements(sigmas) gt 0)
    if n_key gt 1 then $
         message, 'The DEPTHS, LEVELS & SIGMAS keywords cannot be used together'
  endif else begin
    if n_elements(depths) gt 0 then begin
      message, 'The DEPTHS keyword is not required or allowed when ' + $
               'the variable '+variable+' has no vertical dimension'
    endif
    if n_elements(levels) gt 0 then begin
      message, 'The LEVELS keyword is not required or allowed when ' + $
               'the variable '+variable+' has no vertical dimension'
    endif
    if n_elements(levels) gt 0 then begin
      message, 'The SIGMAS keyword is not required or allowed when ' + $
               'the variable '+variable+' has no vertical dimension'
    endif
  endelse

  if has_time then begin
    if n_elements(record) eq 0 then record = 0
    if record lt 0 then $
         record = self->DimInfo(grid.dims.time,/DIMSIZE) + record
  endif else begin
    if n_elements(record) gt 0 then begin
       message, 'The RECORD keyword is not required or allowed when ' + $
                'the variable '+variable+' has no time dimension'
    endif
    if use_bath then begin
       message, 'The USE_BATH option may not be activated when ' + $
                'the variable '+variable+' has no time dimension'
    endif
    if use_zeta then begin
      message, 'The USE_ZETA option may not be activated when ' + $
        'the variable '+variable+' has no time dimension'
    endif
  endelse

  ;; If variable has a vertical dimension, get s-coordinate data

  if has_vert then scoord = self->GetScoord(dims.vertical)

  ;; Establish levels/depths at which data are required

  get_depths = n_elements(depths) gt 0
  get_sigmas = n_elements(sigmas) gt 0
  get_levels = n_elements(levels) gt 0

  case 1B of
    get_depths: n_tran = n_elements(depths)
    get_levels: n_tran = n_elements(levels)
    get_sigmas: n_tran = n_elements(sigmas)
    else: begin
      n_tran = 1
      if has_vert then levels = scoord.n_s-1
    endelse
  endcase

  ;; More handy constants

  ara0 = grid.along_range[0]
  ara1 = grid.along_range[1]
  aran = ara1-ara0+1

  ;; Create array to hold result.

  result = fltarr(aran, n_tran)

  ;; ...horizontal dimension

  case grid.direction of
    0: begin
      offset_h = [ara0,grid.index]
      count_h = [aran,1]
    end
    1: begin
      offset_h = [grid.index,ara0]
      count_h = [1,aran]
    end
  endcase

  delta_h = [0,0]

  case strjoin(dims.horizontal, ' ') of
    'xi_rho eta_rho':
    'xi_u eta_u': begin
      offset_h[0] -= 1
      count_h[0] += 1
      delta_h[0] -= 1
    end
    'xi_v eta_v': begin
      offset_h[1] -= 1
      count_h[1] += 1
      delta_h[1] -= 1
    end
    'xi_psi eta_psi': begin
      offset_h -= 1
      count_h += 1
      delta_h -= 1
    end
  endcase

  ;; Read variable

  if get_depths || get_sigmas then begin

    ;; The DEPTHS or SIGMAS keyword has been specified, so get slices
    ;; and interpolate vertically

    ;; Build up offset & count vectors for VarGet

    offset = [offset_h,0]
    count = [count_h,0]
    delta = [delta_h,0]

    if has_time then begin
      offset = [offset,record]
      count = [count,1]
      delta = [delta,0]
    endif

    ;; Get a slice

    vslice = self->VarGet(variable, OFFSET=offset, COUNT=count)

    vslice = reform(mgh_stagger(vslice, DELTA=delta))

    ;; Get zeta & bathymetry @ variable location

    if keyword_set(use_bath) then begin
      message, "Sorry I don't do USE_bath=1 yet!"
    endif else begin
      bath = grid.h
    endelse

    if keyword_set(use_zeta) then begin
      message, "Sorry I don't do USE_ZETA=1 yet!"
    endif else begin
      zeta = fltarr(aran)
    endelse

    ;; Loop horizontally thru domain interpolating to all depths

    cs = mgh_roms_s_to_cs(grid.s, $
                          THETA_S=grid.theta_s, THETA_B=grid.theta_b, $
                          VSTRETCH=grid.vstretch)

    for i=0,aran-1 do begin
      zz = mgh_roms_s_to_z(grid.s, bath[i], $
                           ZETA=zeta[i], HC=grid.hc, CS=cs, $
                           VTRANSFORM=grid.vtransform)
      if get_sigmas then begin
        zs = sigmas*zeta[i]-(1-sigmas)*bath[i]
        result[i,*] = interpol(vslice[i,*], zz, zs, /QUADRATIC)
      endif else begin
        result[i,*] = interpol(vslice[i,*], zz, -float(depths), /QUADRATIC)
        l_bottom = where(depths gt bath[i], n_bottom)
        if n_bottom gt 0 then result[i,l_bottom] = !values.f_nan
      endelse
    endfor

  endif else begin

    ;; DEPTHS & SIGMAS not specified, so no vertical interpolation
    ;; necessary. Get data for each level separately

    for k=0,n_tran-1 do begin

      ;; Build up offset & count vectors for VarGet

      offset = offset_h
      count = count_h
      delta = delta_h

      if has_vert then begin
        offset = [offset,levels[k]]
        count = [count,1]
        delta = [delta,0]
      endif

      if has_time gt 0 then begin
        offset = [offset,record]
        count = [count,1]
        delta = [delta,0]
      endif

      ;; Get data and load into result array.

      data = self->VarGet(variable, OFFSET=offset, COUNT=count)

      data = reform(mgh_stagger(data, DELTA=delta))

      result[*,k] = data

    endfor

  endelse

  ;; Reset masked values. This is commented out for now because
  ;; there are some tricky issues to do with staggered grids
  ;; and multi-slice averaging that I haven't sorted out yet.

;  if size(grid.mask, /N_DIMENSIONS) gt 0 then begin
;    n_s = grid.average ? 1 : n_slice
;    for s=0,n_s-1 do begin
;      ww = where(grid.mask[*,s] lt 0.01, count)
;      if count gt 0 then begin
;        for k=0,n_tran-1 do begin
;          r = result[*,k,s]
;          r[ww] = mask_value
;          result[0,k,s] = r
;        endfor
;      endif
;    endfor
;  endif

  return, result

end

; MGHromsHistory::DimRho
;
;   Return the dimensions of the rho grid in an array.
;
function MGHromsHistory::DimRho

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   result = [self->DimInfo('xi_rho', /DIMSIZE), $
             self->DimInfo('eta_rho', /DIMSIZE)]

   if self->HasDim('s_rho') then $
        result = [result, self->DimInfo('s_rho', /DIMSIZE)]

   return, result

end

; MGHromsHistory::GetCppOptions
;
;   Return a string array containing the list of CPP options
;
function MGHromsHistory::GetCppOptions, COUNT=count

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   aname = 'CPP_options'

   if self->HasAtt(/GLOBAL, aname) then begin
      s = strtrim(strsplit(self->AttGet(/GLOBAL, aname), ',', /EXTRACT), 2)
      count = n_elements(s)
      return, s
   endif else begin
      count = 0
      return, ''
   endelse

end

; MGHromsHistory::GetScoord
;
function MGHromsHistory::GetScoord, dim_vertical

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

   if strlen(dim_vertical) eq 0 then return, -1

   if ~ self->HasDim('s_rho') then return, -1

   n_s_rho = self->DimInfo('s_rho', /DIMSIZE)

   n = n_s_rho + (dim_vertical eq 's_w')

   theta_s = self->HasVar('theta_s') ? self->VarGet('theta_s') : self->AttGet('theta_s', /GLOBAL)

   theta_b = self->HasVar('theta_b') ? self->VarGet('theta_b') : self->AttGet('theta_b', /GLOBAL)

   hc = self->HasVar('hc') ? self->VarGet('hc') : self->AttGet('hc', /GLOBAL)

   vstretch = self->HasVar('Vstretching') ? self->VarGet('Vstretching'): 1

   vtransform = self->HasVar('Vtransform') ? self->VarGet('Vtransform'): 1

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
   s = (s > (-1)) < 0

   return, {theta_s: theta_s, theta_b: theta_b, hc: hc, $
      vstretch: vstretch, vtransform: vtransform, n: n, s: s}

end

; MGHromsHistory::GetTransportBox
;
function MGHromsHistory::GetTransportBox, $
  ETA_RANGE=eta_range, LONLAT=lonlat, XI_RANGE=xi_range, _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Pass RECORD, VAR_UBAR and VAR_VBAR via inheritance

   ;; Are the grid locations available in (lon,lat) coordinates?

   if n_elements(lonlat) eq 0 then $
      lonlat = self->HasVar('lon_rho') && self->HasVar('lat_rho')

   ;; Specify subset of the grid. Here XI_RANGE and ETA_RANGE define
   ;; a rectangular block of RHO cells; transport is evaluated
   ;; around the *outside* of this block.

   dim_rho = [self->DimInfo('xi_rho', /DIMSIZE),self->DimInfo('eta_rho', /DIMSIZE)]

   if n_elements(xi_range) eq 0 then xi_range = [1,-2]
   if xi_range[0] lt 0 then xi_range[0] += dim_rho[0]
   if xi_range[1] lt 0 then xi_range[1] += dim_rho[0]

   if n_elements(eta_range) eq 0 then eta_range = [1,-2]
   if eta_range[0] lt 0 then eta_range[0] += dim_rho[1]
   if eta_range[1] lt 0 then eta_range[1] += dim_rho[1]

   ;; Get 4 slices, anti-clockwise beginning with southern boundary
   for s=0,3 do begin

      case s of
         0: slice = self->GetTransportSlice(DIRECTION=0, eta_range[0]-1, ALONG_RANGE=xi_range, _STRICT_EXTRA=_extra)
         1: slice = self->GetTransportSlice(DIRECTION=1, xi_range[1], ALONG_RANGE=eta_range, _STRICT_EXTRA=_extra)
         2: slice = self->GetTransportSlice(DIRECTION=0, eta_range[1], ALONG_RANGE=xi_range, _STRICT_EXTRA=_extra)
         3: slice = self->GetTransportSlice(DIRECTION=1, xi_range[0]-1, ALONG_RANGE=eta_range, _STRICT_EXTRA=_extra)
      endcase

      ;; Number of points in slice (incl end points)
      n_slice = n_elements(slice.distance)

      ;; Boundary number for each point in slice
      slice_bb = replicate(s, n_slice-1)

      ;; Difference distance & transport
      slice_dd = slice.distance[1:n_slice-1] - slice.distance[0:n_slice-2]
      slice_dt = slice.transport[1:n_slice-1] - slice.transport[0:n_slice-2]

      ;; Extract depth & horizontal position
      slice_hh = slice.depth
      if lonlat then begin
         slice_xx = slice.lon
         slice_yy = slice.lat
      endif else begin
         slice_xx = slice.x
         slice_yy = slice.y
      endelse

      mgh_undefine, slice

      ;; For slices 1 & 2 positive transport is directed out of the box
      if s eq 1 or s eq 2 then slice_dt = - slice_dt

      ;; For slices 2 & 3 positive direction runs clockwise
      if s eq 2 or s eq 3 then begin
         slice_dd = reverse(slice_dd)
         slice_dt = reverse(slice_dt)
         slice_hh = reverse(slice_hh)
         slice_xx = reverse(slice_xx)
         slice_yy = reverse(slice_yy)
      endif

      if s eq 0 then begin
         box_bb = slice_bb
         box_dd = slice_dd
         box_dt = slice_dt
         box_hh = slice_hh
         box_xx = slice_xx
         box_yy = slice_yy
      endif else begin
         box_bb = [box_bb,slice_bb]
         box_dd = [box_dd,slice_dd]
         box_dt = [box_dt,slice_dt]
         box_hh = [box_hh,slice_hh]
         box_xx = [box_xx,slice_xx]
         box_yy = [box_yy,slice_yy]
      endelse

   endfor

   n_box = n_elements(box_dd)

   if lonlat then begin
      result = {boundary: box_bb, depth: box_hh, lon: box_xx, lat: box_yy, distance: dblarr(n_box+1), transport: dblarr(n_box+1)}
   endif else begin
      result = {boundary: box_bb, depth: box_hh, x: box_xx, y: box_yy, distance: dblarr(n_box+1), transport: dblarr(n_box+1)}
   endelse

   result.distance[0] = 0
   for i=0,n_box-1 do $
      result.distance[i+1] = result.distance[i] + box_dd[i]

   result.transport[0] = 0
   for i=0,n_box-1 do $
      result.transport[i+1] = result.transport[i] + box_dt[i]

   return, result

end

; MGHromsHistory::GetTransportSlice
;
;   Calculate & return transport normal to the specified slice.
;
;   The slices in question are vertical slices along the grid through
;   normal-velocity points. I used to call this method
;   GetTransportCslice, and it does use the Cslice terminology for specifying
;   the slice location. However it actually uses none of the Cslice code
;   and I am thinking of redefining a Cslice to be always along rho
;   points, in which case the use of the term Cslice would be misleading.
;
;   TO DO:
;     - Include zeta in depth?
;
function MGHromsHistory::GetTransportSlice, slice , $
   ALONG_RANGE=along_range, DIRECTION=direction, LONLAT=lonlat, RECORD=record, $
   VAR_UBAR=var_ubar, VAR_VBAR=var_vbar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(var_ubar) eq 0 then var_ubar = 'ubar'
   if n_elements(var_vbar) eq 0 then var_vbar = 'vbar'

   if n_elements(record) eq 0 then record = 0

   ;; Are the grid locations available in (lon,lat) coordinates?
   if n_elements(lonlat) eq 0 then $
      lonlat = self->HasVar('lon_rho') && self->HasVar('lat_rho')

   ;; Default slice orientation is east-west
   if n_elements(direction) eq 0 then direction = 0

   ;; Default slice is through the PSI points on the southern boundary
   ;; (western if DIRECTION is 1)
   if n_elements(slice) eq 0 then slice = 0

   ;; Default along_range is from western/southern to eastern/northern
   ;; boundary
   if n_elements(along_range) eq 0 then begin
      n_rho = self->DimRho()
      case direction of
         0: along_range = [1,n_rho[0]-2]
         1: along_range = [1,n_rho[1]-2]
      endcase
   endif

   ;; Abbreviations for along_range
   ar0 = along_range[0]
   ar1 = along_range[1]
   arn = ar1 - ar0 + 1

   ;; Create a structure to hold the result

   if lonlat then begin
      result = {distance: dblarr(arn+1), transport: dblarr(arn+1), depth: fltarr(arn+1), lon: dblarr(arn), lat: dblarr(arn)}
   endif else begin
      result = {distance: dblarr(arn+1), transport: dblarr(arn+1), depth: fltarr(arn+1), x: dblarr(arn), y: dblarr(arn)}
   endelse

   ;; Retrieve depth at cell centres. Average over the rho slices on
   ;; either side
   case direction of
      0: hh = mgh_avg(self->VarGet('h', OFFSET=[ar0,slice], COUNT=[arn,2]),2)
      1: hh = mgh_avg(self->VarGet('h', OFFSET=[slice,ar0], COUNT=[2,arn]),1)
   endcase

   ;; Retrieve the relevant metric coefficients at cell
   ;; centres. Average over the rho slices on either side
   case direction of
      0: pp = mgh_avg(self->VarGet('pm', OFFSET=[ar0,slice], COUNT=[arn,2]),2)
      1: pp = mgh_avg(self->VarGet('pn', OFFSET=[slice,ar0], COUNT=[2,arn]),1)
   endcase

   ;; Retrieve horizontal position at cell centres. Average over the rho
   ;; slices on either side
   case direction of
      0: begin
         xx = mgh_avg(self->VarGet((lonlat ? 'lon_rho' : 'x_rho'), OFFSET=[ar0,slice], COUNT=[arn,2]),2)
         yy = mgh_avg(self->VarGet((lonlat ? 'lat_rho' : 'y_rho'), OFFSET=[ar0,slice], COUNT=[arn,2]),2)
      end
      1: begin
         xx = mgh_avg(self->VarGet((lonlat ? 'lon_rho' : 'x_rho'), OFFSET=[slice,ar0], COUNT=[2,arn]),1)
         yy = mgh_avg(self->VarGet((lonlat ? 'lat_rho' : 'y_rho'), OFFSET=[slice,ar0], COUNT=[2,arn]),1)
      end
   endcase

   ;; Retrieve the velocity mask at cell centres
   mm = 1
   case direction of
      0: begin
         if self->HasVar('mask_v') then $
            mm = reform(self->VarGet('mask_v', OFFSET=[ar0,slice], COUNT=[arn,1]))
      end
      1: begin
         if self->HasVar('mask_u') then $
            mm = reform(self->VarGet('mask_u', OFFSET=[slice,ar0], COUNT=[1,arn]))
      end
   endcase

   ;; Retrieve and mask the slice-normal velocity at cell centres
   case direction of
      0: begin
         vv = mm * reform(self->VarGet(var_vbar, OFFSET=[ar0,slice,record], COUNT=[arn,1,1], /AUTOSCALE))
      end
      1: begin
         vv = mm * reform(self->VarGet(var_ubar, OFFSET=[slice,ar0,record], COUNT=[1,arn,1], /AUTOSCALE))
      end
   endcase

   ;; Calculate along-slice distance at cell faces from origin (by default
   ;; the psi point at the southwest corner of the domain).
   result.distance[0] = 0
   for i=0,arn-1 do $
      result.distance[i+1] = result.distance[i] + 1/pp[i]

   ;; Calculate cumulative transport at cell faces
   result.transport[0] = 0
   for i=0,arn-1 do $
      result.transport[i+1] = result.transport[i] + vv[i]*hh[i]/pp[i]

   ;; Mask & store the bathymetry
   ww = where(mm eq 0, n_masked)
   if n_masked gt 0 then hh[ww] = !values.f_nan
   result.depth = hh

   ;; Store horizontal position at cell centres
   if lonlat then begin
      result.lon = xx
      result.lat = yy
   endif else begin
      result.x = xx
      result.y = yy
   endelse

   return, result

end

; MGHromsHistory::GetTransportPslice
;
;   Calculate & return transport normal to the specified P-slice.
;
;   The transport increases along the Pslice where the velocity is
;   from right to left.
;
function MGHromsHistory::GetTransportPslice, $
   GRID=grid, RECORD=record, VAR_UBAR=var_ubar, VAR_VBAR=var_vbar, VAR_ZETA=var_zeta, USE_ZETA=use_zeta

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(grid) eq 0 then grid = self->PsliceGrid()

   if n_elements(record) eq 0 then record = 0

   if n_elements(var_ubar) eq 0 then var_ubar = 'ubar'
   if n_elements(var_vbar) eq 0 then var_vbar = 'vbar'
   if n_elements(var_zeta) eq 0 then var_zeta = 'zeta'

   if n_elements(use_zeta) eq 0 then use_zeta = !false

   ;; Get ubar & vbar data for a rectangular region spanning the
   ;; X-slice

   xr = [floor(min(grid.point_xi_rho)),ceil(max(grid.point_xi_rho))]
   er = [floor(min(grid.point_eta_rho)),ceil(max(grid.point_eta_rho))]

   xr0 = xr[0]  &  xrn = xr[1]-xr[0]+1
   er0 = er[0]  &  ern = er[1]-er[0]+1

   uv = self->VectorGet([var_ubar,var_vbar], OFFSET=[xr0,er0,record], COUNT=[xrn,ern,1])

   message, 'OK, time to start thinking!'

   ;; Interpolate to slice interval mid-points

   xx = mgh_stagger(grid.xi, DELTA=[-1])-xr0
   ee = mgh_stagger(grid.eta, DELTA=-1)-er0

   uv = interpolate(temporary(uv), xx, ee)

   ;; Rotate into slice-relative coordinates.

   uv = temporary(uv)*exp(-!const.i*grid.angle)

   ;; Calculate depth at interval mid-points

   h = mgh_stagger(grid.h, DELTA=[-1])

   if keyword_set(use_zeta) then begin
      zeta = self->VarGet(var_zeta, OFFSET=[xr0,er0,record], COUNT=[xrn,ern,1])
      zeta = mgh_fill2d(zeta)
      xx = mgh_stagger(grid.xi, DELTA=[-1])-xr0
      ee = mgh_stagger(grid.eta, DELTA=-1)-er0
      h += interpolate(temporary(zeta), xx, ee)
   endif

   ;; The normal velocity is the imaginary part of the complex velocity

   result = dblarr(grid.n_points)
   for i=1,grid.n_points-1 do begin
      result[i] = result[i-1] + imaginary(uv[i-1])*h[i-1]*(grid.arc[i]-grid.arc[i-1])
   endfor

   return, result

end

; MGHromsHistory::GetTransportXslice
;
;   Calculate & return transport normal to the specified X-slice.
;
;   The transport increases along the Xslice where the velocity is
;   from right to left.
;
function MGHromsHistory::GetTransportXslice, $
  GRID=grid, RECORD=record, VAR_UBAR=var_ubar, VAR_VBAR=var_vbar, VAR_ZETA=var_zeta, USE_ZETA=use_zeta

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(grid) eq 0 then grid = self->XsliceGrid()

   if n_elements(record) eq 0 then record = 0

   if n_elements(var_ubar) eq 0 then var_ubar = 'ubar'
   if n_elements(var_vbar) eq 0 then var_vbar = 'vbar'
   if n_elements(var_zeta) eq 0 then var_zeta = 'zeta'

   if n_elements(use_zeta) eq 0 then use_zeta = !false

   ;; Get ubar & vbar data for a rectangular region spanning the
   ;; X-slice

   xr = [floor(min(grid.xi)),ceil(max(grid.xi))]
   er = [floor(min(grid.eta)),ceil(max(grid.eta))]

   xr0 = xr[0]  &  xrn = xr[1]-xr[0]+1
   er0 = er[0]  &  ern = er[1]-er[0]+1

   uv = self->VectorGet([var_ubar,var_vbar], OFFSET=[xr0,er0,record], COUNT=[xrn,ern,1])

   ;; Interpolate to slice interval mid-points

   xx = mgh_stagger(grid.xi, DELTA=[-1])-xr0
   ee = mgh_stagger(grid.eta, DELTA=-1)-er0

   uv = interpolate(temporary(uv), xx, ee)

   ;; Rotate into slice-relative coordinates.

   uv = temporary(uv)*exp(-!const.i*grid.angle)

   ;; Calculate depth at interval mid-points

   h = mgh_stagger(grid.h, DELTA=[-1])

   if keyword_set(use_zeta) then begin
      zeta = self->VarGet(var_zeta, OFFSET=[xr0,er0,record], COUNT=[xrn,ern,1])
      zeta = mgh_fill2d(zeta)
      xx = mgh_stagger(grid.xi, DELTA=[-1])-xr0
      ee = mgh_stagger(grid.eta, DELTA=-1)-er0
      h += interpolate(temporary(zeta), xx, ee)
   endif

   ;; The normal velocity is the imaginary part of the complex velocity

   result = dblarr(grid.n_points)
   for i=1,grid.n_points-1 do begin
      result[i] = result[i-1] + imaginary(uv[i-1])*h[i-1]*(grid.arc[i]-grid.arc[i-1])
   endfor

   return, result

end

function MGHromsHistory::GetWalls

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  vars = self->VarNames()

  opts = self->GetCppOptions()

  dum = where(vars eq 'wall1',count)

  if count gt 0 then begin

    ;; The SCRUM style of specifying walls (may not work for all
    ;; versions) Assume that if 'wall1' exists then so do 'wall2',
    ;; 'wall3', 'wall4'

    walls = strarr(4)
    for w=0,n_elements(walls)-1 do begin
      ww = 'wall'+strtrim(w+1,2)
      if max(vars eq ww) gt 0 then walls[w] = string(self->VarGet(ww))
    endfor

    return, walls eq 'T' or walls eq 't'

  endif else begin

    ;; The ROMS style of specifying walls

    return,  [max(strcmp(opts,'WESTERN_WALL' )) gt 0, $
      max(strcmp(opts,'SOUTHERN_WALL')) gt 0, $
      max(strcmp(opts,'EASTERN_WALL' )) gt 0, $
      max(strcmp(opts,'NORTHERN_WALL')) gt 0]

  endelse

end

; MGHromsHistory::GetZGrid
;
;   Return a 3-D array of Z values on the model grid
;
function MGHromsHistory::GetZGrid, $
  HORIZONTAL=horizontal, VERTICAL=vertical, ZETA=zeta

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  compile_opt OBSOLETE

  if n_elements(horizontal) eq 0 then horizontal = 'rho'

  if n_elements(vertical) eq 0 then vertical = 'rho'

  ;; Retrieve bathymetry and interpolate to output grid

  h = self->VarGet('h')

  h = mgh_roms_stagger(h, TO=horizontal)

  dim_hor = size(h, /DIMENSIONS)

  ;; Process ZETA argument and interpolate to output grid

  if n_elements(zeta) eq 0 then begin
    dim_rho = self->DimRho()
    zeta = make_array(DIMENSION=dim_rho[0:1])
  endif

  zzz = mgh_roms_stagger(zeta, TO=horizontal)

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

  result = make_array(DIMENSION=[dim_hor, n_s])

  cs = mgh_roms_s_to_cs(s, THETA_S=theta_s, THETA_B=theta_b, VSTRETCH=vstretch)

  for i=0,dim_hor[0]-1 do begin
    for j=0,dim_hor[1]-1 do begin
      result[i,j,*] = mgh_roms_s_to_z(s, h[i,j], ZETA=zzz[i,j], CS=cs, HC=hc, VTRANSFORM=vtransform)
    endfor
  endfor

  return, result

end

; MGHromsHistory::HasVar
;
function MGHromsHistory::HasVar, var

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

      isa(var, 'STRING') && (var eq 'Dsbl'): begin
         result = self->HasVar('Hsbl')
      end

      isa(var, 'STRING') && (var eq 'Dbbl'): begin
         result = self->HasVar('Hbbl') && self->HasVar('h')
      end

      else: begin
         result = self->MGHncSequence::HasVar(var)
      end

   endcase

   return, result

end

; ** H-slice methods *******************************************
;
; MGHromsHistory::HsliceData
;
function MGHromsHistory::HsliceData, var, $
     DEPTHS=depths, LAYERS=layers, LEVELS=levels, SIGMAS=sigmas, $
     GRID=grid, LONLAT=lonlat, MASK_VALUE=mask_value, RECORD=record, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
     USE_BATH=use_bath, USE_ZETA=use_zeta

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Check for problems with inputs

   if n_elements(var) ne 1 then $
        message, 'The name of a variable must be supplied'

   if ~ (isa(var, 'STRING') || isa(var, 'STRUCT'))  then $
        message, 'A variable identifier must be supplied'

   ;; Set defaults

   if n_elements(mask_value) eq 0 then mask_value = !values.f_nan

   if n_elements(depths) gt 0 || n_elements(sigmas) gt 0 then begin
      if n_elements(use_bath) eq 0 then use_bath = self->HasVar('bath')
      if n_elements(use_zeta) eq 0 then use_zeta = self->HasVar('zeta')
   endif

   ;; If no grid information is supplied, then get it. Otherwise check
   ;; grid is consistent with the current variable name & function
   ;; arguments.

   if n_elements(grid) gt 0 then begin
      dims = self->VarDims(var)
      if ~ array_equal(dims.horizontal, grid.dims.horizontal) gt 0 then $
           message, 'Horizontal dimensions of variable do not match GRID data'
      if dims.bed ne grid.dims.bed then $
           message, 'Vertical dimension of variable does not match GRID data'
      if dims.vertical ne grid.dims.vertical then $
           message, 'Vertical dimension of variable does not match GRID data'
      if dims.time ne grid.dims.time then $
           message, 'Time dimension of variable does not match GRID data'
      mgh_undefine, dims
      if n_elements(xi_range) gt 0 && ~ array_equal(xi_range,grid.xi_range) then $
           message, 'XI_RANGE does not match grid data'
      if n_elements(eta_range) gt 0 && ~ array_equal(eta_range,grid.eta_range) then $
           message, 'ETA_RANGE does not match grid data'
   endif else begin
      grid = self->HsliceGrid(var, ETA_RANGE=eta_range, LONLAT=lonlat, XI_RANGE=xi_range)
   endelse

   ;; Check consistency of function arguments with variable dimensions

   if grid.dims.vertical then begin
      n_key = (n_elements(depths) gt 0) + (n_elements(levels) gt 0) + $
              (n_elements(sigmas) gt 0)
      if n_key gt 1 then $
           message, 'The DEPTHS, LEVELS & SIGMAS keywords cannot be used together'
   endif else begin
      fmt = '(%"The %s keyword is not required or allowed when the ' + $
            'variable %s has no vertical dimension")'
      if n_elements(levels) gt 0 then message, string(FORMAT=fmt, 'LEVELS', var)
      if n_elements(depths) gt 0 then message, string(FORMAT=fmt, 'DEPTHS', var)
      if n_elements(sigmas) gt 0 then message, string(FORMAT=fmt, 'SIGMAS', var)
   endelse

   if grid.dims.bed then begin
      ;; Actually, I haven't thought of anything to test here
   endif else begin
      fmt = '(%"The %s keyword is not required or allowed when the ' + $
            'variable %s has no vertical dimension")'
      if n_elements(layers) gt 0 then message, string(FORMAT=fmt, 'LAYERS', var)
   endelse

   if grid.dims.time then begin
      if n_elements(record) eq 0 then record = 0
      if record lt 0 then $
           record = self->DimInfo(grid.dims.time, /DIMSIZE) + record
   endif else begin
      if n_elements(record) gt 0 then begin
         message, 'The RECORD keyword is not required or allowed when ' + $
                  'the variable '+var+' has no time dimension'
      endif
      if keyword_set(use_bath) then begin
         message, 'The USE_BATH option may not be activated when ' + $
                  'the variable '+var+' has no time dimension'
      endif
      if keyword_set(use_zeta) then begin
         message, 'The USE_ZETA option may not be activated when ' + $
                  'the variable '+var+' has no time dimension'
      endif
   endelse

   ;; Abbreviations for xi and eta range

   xra0 = grid.xi_range[0]
   xra1 = grid.xi_range[1]
   xran = xra1-xra0+1
   era0 = grid.eta_range[0]
   era1 = grid.eta_range[1]
   eran = era1-era0+1

   ;; If vertical interpolation is required, get bath & zeta at the
   ;; variable location

   if n_elements(depths) gt 0 || n_elements(sigmas) gt 0 then begin

      if use_bath || use_zeta then begin
         offset = [xra0,era0,record]
         count = [xran,eran,1]
         delta = [0,0]
         hdims = strjoin(grid.dims.horizontal, ' ')
         case !true of
            strmatch(hdims, 'xi_rho eta_rho'):
            strmatch(hdims, 'xi_u eta_u') || strmatch(hdims, 'xi_u eta_rho'): begin
               count += [1,0,0]
               delta -= [1,0]
            end
            strmatch(hdims, 'xi_v eta_v') || strmatch(hdims, 'xi_rho eta_v'): begin
               count += [0,1,0]
               delta -= [0,1]
            end
            strmatch(hdims, 'xi_psi eta_psi'): begin
               count += [1,1,0]
               delta -= [1,1]
            end
         endcase
      endif

      if use_bath then begin
         bath = reform(self->VarGet('bath', OFFSET=offset, COUNT=count))
         bath = mgh_stagger(bath, DELTA=delta)
      endif else begin
         if size(grid.h, /N_DIMENSIONS) ne 2 then $
              message, 'Static bathymetry data needed but not found'
         bath = grid.h
      endelse

      if use_zeta then begin
         zeta = reform(self->VarGet('zeta', OFFSET=offset, COUNT=count))
         zeta = mgh_stagger(zeta, DELTA=delta)
      endif else begin
         zeta = fltarr(xran, eran)
      endelse

   endif

   ;; Read & unpack variable, interpolating if necessary

   case !true of

      n_elements(depths) gt 0: begin

         ;; DEPTHS keyword has been set so interpolate to constant-z levels

         n_slice = n_elements(depths)

         ;; DEPTH specified, so get 3D data & interpolate vertically.
         ;; Since vertical interpolation is required we make a result
         ;; array of floating point type immediately.

         result = replicate(!values.f_nan, xran, eran, n_slice)

         ;; Build up OFFSET & COUNT vectors for the netCDF get
         ;; operation

         offset = [xra0,era0,0]
         count = [xran,eran,grid.n_level]

         if strlen(grid.dims.time) gt 0 then begin
            offset = [offset, record]
            count = [count, 1]
         endif

         ;; Get 3D data.

         var3d = self->VarGet(var, OFFSET=offset, COUNT=count, /AUTOSCALE)

         ;; Loop horizontally thru domain interpolating to z
         ;; Clip s to avoid slightly, out-of-bounds values generated
         ;; by netCDF file packing.

         cs = mgh_roms_s_to_cs(grid.s, $
                               THETA_S=grid.theta_s, THETA_B=grid.theta_b, $
                               VSTRETCH=grid.vstretch)

         for i=0,xran-1 do begin
            for j=0,eran-1 do begin
               zz = mgh_roms_s_to_z(grid.s, bath[i,j], $
                                    CS=cs, ZETA=zeta[i,j], HC=grid.hc, $
                                    VTRANSFORM=grid.vtransform)
               varss = reform(var3d[i,j,*])
               if min(finite(varss)) gt 0 then begin
                  varzz = interpol(varss, zz, -float(depths), /SPLINE)
                  varzz[where(depths gt bath[i,j], /NULL)] = !values.f_nan
                  result[i,j,*] = varzz
               endif
            endfor
         endfor

         !null = !null

      end

      n_elements(sigmas) gt 0: begin

         ;; SIGMAS keyword has been set, so interpolate to constant-sigma levels

         n_slice = n_elements(sigmas)

         ;; SIGMAS specified, so get 3D data & interpolate vertically.
         ;; Since vertical interpolation is required we make a result
         ;; array of floating point type immediately.

         result = fltarr(xran, eran, n_slice)

         ;; Build up OFFSET & COUNT vectors for the netCDF get
         ;; operation

         offset = [xra0,era0,0]
         count = [xran,eran,n_elements(grid.s)]

         if strlen(grid.dims.time) gt 0 then begin
            offset = [offset, record]
            count = [count, 1]
         endif

         ;; Get 3D data.

         var3d = self->VarGet(var, OFFSET=offset, COUNT=count, /AUTOSCALE)

         ;; Loop horizontally thru domain interpolating to specified sigma:
         ;; Clip s to avoid slightly out-of-bounds values generated
         ;; by netCDF file packing.

         cs = mgh_roms_s_to_cs(grid.s, $
                               THETA_S=grid.theta_s, THETA_B=grid.theta_b, $
                               VSTRETCH=grid.vstretch)

         for i=0,xran-1 do begin
            for j=0,eran-1 do begin
               zz = mgh_roms_s_to_z(grid.s, bath[i,j], $
                                    ZETA=zeta[i,j], HC=grid.hc, CS=cs, $
                                    VTRANSFORM=grid.vtransform)
               zsij = sigmas*zeta[i,j]-(1-sigmas)*bath[i,j]
               varss = interpol(var3d[i,j,*], zz, zsij, /SPLINE)
               result[i,j,*] = varss
               !null = check_math()
            endfor
         endfor

      end

      else: begin

         ;; Neither DEPTHS nor SIGMAS specified, so vertical interpolation
         ;; is not necessary

         case 1B of
            strlen(grid.dims.vertical) gt 0: begin
               slice = n_elements(levels) gt 0 ? levels : grid.n_level-1
               n_slice = n_elements(slice)
            end
            strlen(grid.dims.bed) gt 0: begin
               ;; Note that the bed layer index increases downward
               slice = n_elements(layers) gt 0 ? layers : 0
               n_slice = n_elements(slice)
            end
            else: begin
               n_slice = 1
            end
         endcase

         result = fltarr(xran, eran, n_slice)

         for k=0,n_slice-1 do begin

            offset = [xra0,era0]
            count = [xran,eran]

            if grid.dims.vertical || grid.dims.bed then begin
               offset = [offset, slice[k]]
               count = [count, 1]
            endif

            if grid.dims.time then begin
               offset = [offset, record]
               count = [count, 1]
            endif

            var2d = self->VarGet(var, OFFSET=offset, COUNT=count, /AUTOSCALE)

            result[*,*,k] = temporary(var2d)

         endfor

      end

   endcase

   ;; Reset masked values

   if size(grid.mask, /N_DIMENSIONS) eq 2 then begin

      l_masked = where(grid.mask lt 0.5, n_masked)

      if n_masked gt 0 then begin
         for k=0,n_slice-1 do begin
            r = result[*,*,k]
            r[l_masked] = mask_value
            result[0,0,k] = r
         endfor
      endif

   end

   return, result

end

; MGHromsHistory::HsliceGrid
;
function MGHromsHistory::HsliceGrid, arg, $
     ETA_RANGE=eta_range, LONLAT=lonlat, XI_RANGE=xi_range

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(arg) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'arg'

   dims = self->VarDims(arg)

   ;; Get dimensions for RHO points

   dim_rho = self->DimRho()

   ;; Establish horizontal range. Default is to retrieve all data.

   hdims = strjoin(dims.horizontal, ' ')
   case !true of
      strmatch(hdims, 'xi_rho eta_rho'): begin
         if n_elements(xi_range) eq 0 then xi_range = [0,dim_rho[0]-1]
         if n_elements(eta_range) eq 0 then eta_range = [0,dim_rho[1]-1]
      end
      strmatch(hdims, 'xi_u eta_u') || strmatch(hdims, 'xi_u eta_rho'): begin
         if n_elements(xi_range) eq 0 then xi_range = [0,dim_rho[0]-2]
         if n_elements(eta_range) eq 0 then eta_range = [0,dim_rho[1]-1]
      end
      strmatch(hdims, 'xi_v eta_v') || strmatch(hdims, 'xi_rho eta_v'): begin
         if n_elements(xi_range) eq 0 then xi_range = [0,dim_rho[0]-1]
         if n_elements(eta_range) eq 0 then eta_range = [0,dim_rho[1]-2]
      end
      strmatch(hdims, 'xi_psi eta_psi'): begin
         if n_elements(xi_range) eq 0 then xi_range = [0,dim_rho[0]-2]
         if n_elements(eta_range) eq 0 then eta_range = [0,dim_rho[1]-2]
      end
   endcase

   ;; Interpret negative values in XI_RANGE and ETA_RANGE as offsets
   ;; from the end of the grid.

   if xi_range[0] lt 0 then xi_range[0] += dim_rho[0]
   if xi_range[1] lt 0 then xi_range[1] += dim_rho[0]
   if eta_range[0] lt 0 then eta_range[0] += dim_rho[1]
   if eta_range[1] lt 0 then eta_range[1] += dim_rho[1]

   ;; Are the grid locations available in (lon,lat) coordinates?

   if n_elements(lonlat) eq 0 then $
        lonlat = self->HasVar('lon_rho') && self->HasVar('lat_rho')

   ;; Calculate parameters for retrieving & processing horizontal grid
   ;; data, which are defined on rho points.

   offset = [xi_range[0],eta_range[0]]
   count = [xi_range[1]-xi_range[0]+1,eta_range[1]-eta_range[0]+1]
   delta = [0,0]

   hdims = strjoin(dims.horizontal, ' ')
   case !true of
      strmatch(hdims, 'xi_rho eta_rho'):
      strmatch(hdims, 'xi_u eta_u') || strmatch(hdims, 'xi_u eta_rho'): begin
         count += [1,0]
         delta -= [1,0]
      end
      strmatch(hdims, 'xi_v eta_v') || strmatch(hdims, 'xi_rho eta_v'): begin
         count += [0,1]
         delta -= [0,1]
      end
      strmatch(hdims, 'xi_psi eta_psi'): begin
         count += [1,1]
         delta -= [1,1]
      end
   endcase

   ;; Retrieve & interpolate horizontal position

   x = lonlat $
       ? self->VarGet('lon_rho', OFFSET=offset, COUNT=count) $
       : self->VarGet('x_rho', OFFSET=offset, COUNT=count)
   y = lonlat $
       ? self->VarGet('lat_rho', OFFSET=offset, COUNT=count) $
       : self->VarGet('y_rho', OFFSET=offset, COUNT=count)

   x = mgh_stagger(x, DELTA=delta)
   y = mgh_stagger(y, DELTA=delta)

   ;; Retrieve & interpolate mask data

   mask = -1
   if self->HasVar('mask_rho') then begin
      mask = self->VarGet('mask_rho', OFFSET=offset, COUNT=count)
      if delta[0] eq -1 then $
           mask = mask[0:count[0]-2,*] * mask[1:count[0]-1,*]
      if delta[1] eq -1 then $
           mask = mask[*,0:count[1]-2] * mask[*,1:count[1]-1]
   endif

   ;; Retrieve & interpolate angle data

   angle = -1
   if self->HasVar('angle') then begin
      angle = self->VarGet('angle', OFFSET=offset, COUNT=count)
      angle = mgh_stagger(angle, DELTA=delta)
   endif

   ;; Horizontal grid specification is now complete.

   result = dictionary()

   result.dims = dims
   result.xi_range = xi_range
   result.eta_range = eta_range
   result.lonlat = lonlat
   result.x = x
   result.y = y
   result.mask = mask
   result.angle = angle

   case !true of

      strlen(dims.vertical) gt 0: begin

         ;; For a variable with a vertical dimension we return
         ;; s-coordinate data

         ;; Retrieve & interpolate static bathymetry (if available)

         if self->HasVar('h') then begin
            h = self->VarGet('h', OFFSET=offset, COUNT=count)
            result.h = mgh_stagger(h, DELTA=delta)
         endif

         ;; Retrieve s-coordinate parameters and s values for the variable. Note that
         ;; there is similar, but not identical, code in the GetScoord method.

         result.theta_s = self->HasVar('theta_s') ? self->VarGet('theta_s') : self->AttGet('theta_s', /GLOBAL)
         result.theta_b = self->HasVar('theta_b') ? self->VarGet('theta_b') : self->AttGet('theta_b', /GLOBAL)
         result.hc = self->HasVar('hc') ? self->VarGet('hc') : self->AttGet('hc', /GLOBAL)
         result.vstretch = self->HasVar('Vstretching') ? self->VarGet('Vstretching'): 1
         result.vtransform = self->HasVar('Vtransform') ? self->VarGet('Vtransform'): 1

         ;; The following code is complicated by a couple of changes in ROMS
         ;; output files:
         ;;  - In older files, the bottom value of s_w is omitted.
         ;;  - In ROMS 2.1 and earlier, the names of the s-coordinate variables did
         ;;    not match the corresponding dimensions--heavens knows why. In ROMS 2.2
         ;;    this was fixed.

         if self->HasVar(dims.vertical) then begin
            s = (self->VarGet(dims.vertical) < 0) > (-1)
         endif else begin
            n_s_rho = self->DimInfo('s_rho', /DIMSIZE)
            case dims.vertical of
               's_rho': begin
                  if self->HasVar('sc_r') then begin
                     s =  self->VarGet('sc_r')
                  endif else begin
                     s = mgh_stagger(mgh_range(-1, 0, STRIDE=1.0D/n_s_rho), DELTA=-1)
                  endelse
               end
               's_w': begin
                  if self->HasVar('sc_w') then begin
                     s = (self->VarGet('sc_w') < 0) > (-1)
                  endif else begin
                     s = mgh_range(-1, 0, STRIDE=1.0D/n_s_rho)
                  endelse
                  if n_elements(s) eq n_s_rho then s = [-1,s]
               end
            endcase
         endelse

         result.s = s
         result.n_level = n_elements(s)

      end

      strlen(dims.bed) gt 0: begin

         ;; For a variable with a bed-layer dimension we return
         ;; the number of layers

         result.n_layer = self->DimInfo(dims.bed, /DIMSIZE)

      end

      else:

   endcase

   return, result->ToStruct(/RECURSIVE, /NO_COPY)

end

; MGHromsHistory::HsliceMean
;
function MGHromsHistory::HsliceMean, var, $
     ETA_RANGE=eta_range, GRID=grid, LONLAT=lonlat, $
     RECORD_RANGE=record_range, TIME_RANGE=time_range, XI_RANGE=xi_range, $
     _REF_eXTRA=_extra

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  ;; Check for problems with inputs

  if n_elements(var) ne 1 then $
     message, 'The name of a variable must be supplied'

  if size(var, /TNAME) ne 'STRING' then $
     message, 'The name of a variable must be supplied'

  ;; If no grid information is supplied, then get it. Otherwise check
  ;; grid is consistent with the current variable name & function
  ;; arguments.

  if n_elements(grid) gt 0 then begin
     dims = self->VarDims(var)
     if ~ array_equal(dims.horizontal, grid.dims.horizontal) gt 0 then $
        message, 'Horizontal dimensions of variable do not match GRID data'
     if dims.vertical ne grid.dims.vertical then $
        message, 'Vertical dimension of variable does not match GRID data'
     if dims.time ne grid.dims.time then $
        message, 'Time dimension of variable does not match GRID data'
     mgh_undefine, dims
     if n_elements(xi_range) gt 0 && ~ array_equal(xi_range,grid.xi_range) then $
        message, 'XI_RANGE does not match grid data'
     if n_elements(eta_range) gt 0 && ~ array_equal(eta_range,grid.eta_range) then $
        message, 'ETA_RANGE does not match grid data'
  endif else begin
     grid = self->HsliceGrid(var, ETA_RANGE=eta_range, LONLAT=lonlat, XI_RANGE=xi_range)
  endelse

  ;; Select records over which average is to be taken

  ;; Establish records to be processed

  has_time = strlen(grid.dims.time) gt 0

  if has_time then begin
     case !true of
        self->HasVar(grid.dims.time): $
           time_var = grid.dims.time
        self->HasVar('ocean_time'): $
           time_var = 'ocean_time'
        self->HasVar('scrum_time'): $
           time_var = 'scrum_time'
        else: $
           message, 'Time variable not found'
     endcase
     time = self->VarGet(time_var)
     if n_elements(time_range) gt 0 then begin
        record_range = mgh_subset(time, time_range)
     endif
     if n_elements(record_range) eq 0 then begin
        n_time = self->DimInfo(grid.dims.time, /DIMSIZE)
        record_range = [0,n_time-1]
     endif
     if n_elements(time_range) eq 0 then $
        time_range = time[record_range]
     msg = ['Getting', var, 'data between records', strtrim(record_range,2), $
        'times', mgh_format_float(time[record_range])]
     message, /INFORM, strjoin(temporary(msg), ' ')
     rra0 = record_range[0]
     rra1 = record_range[1]
     rran = rra1-rra0+1
  endif else begin
     rran = 1
     msg = ['Variable', self.var, 'does not vary with time']
     message, /INFORM, strjoin(temporary(msg), ' ')
  endelse

  ;; Get first slice

  result = self->HsliceData(var, ETA_RANGE=eta_range, GRID=grid, $
     LONLAT=lonlat, RECORD=rra0, XI_RANGE=xi_range, $
     _STRICT_EXTRA=_extra)

  if rran gt 0 then begin
     for r=rra0+1,rra1 do begin
        result += self->HsliceData(var, ETA_RANGE=eta_range, GRID=grid, $
           LONLAT=lonlat, RECORD=r, XI_RANGE=xi_range, $
           _STRICT_EXTRA=_extra)
     endfor
     result /= float(rran)
  endif

  return, result

end

; MGHromsHistory::LocateXY
;
function MGHromsHistory::LocateXY, x, y, LONLAT=lonlat

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

; ** P-slice and P-transect methods *******************************************

; MGHromsHistory::PsliceDefGrid
;
pro MGHromsHistory::PsliceDefGrid, $ $
   GRID=grid, LONLAT=lonlat, VERT_XI=vert_xi, VERT_ETA=vert_eta

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; If grid information is supplied, then check it is consistent with the
   ;; current function arguments; otherwise get the grid information.

   if n_elements(grid) gt 0 then begin

      if n_elements(lonlat) gt 0  && ~ array_equal(lonlat, grid.lonlat) then $
         message, 'LONLAT does not match grid data'
      if n_elements(vertx) gt 0  && ~ array_equal(vertx, grid.vertx) then $
         message, 'VERT_XI does not match grid data'
      if n_elements(verty) gt 0  && ~ array_equal(verty, grid.verty) then $
         message, 'VERT_ETA does not match grid data'

   endif else begin

      grid = self->PsliceGrid(LONLAT=lonlat, VERT_XI=vert_xi, VERT_ETA=vert_eta)

   endelse

end

; MGHromsHistory::PsliceGrid
;
function MGHromsHistory::PsliceGrid, $
   LONLAT=lonlat, VERT_XI=vert_xi, VERT_ETA=vert_eta

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(lonlat) eq 0 then $
      lonlat = self->HasVar('lon_rho') && self->HasVar('lat_rho')

   ;; The default P-slice is the "southern" normal velocity model boundary,
   ;; as for the GetTransportSlice method.

   dim_rho = self->DimRho()

   if n_elements(vert_xi) eq 0 then $
      vert_xi = [0,dim_rho[0]-2]
   if n_elements(vert_eta) eq 0 then $
      vert_eta = replicate(0, n_elements(vertx))

   ;; Check slice consistency

   n_vert = n_elements(vert_xi)

   if n_elements(vert_eta) ne n_vert then message, 'Inconsistent vertex lengths'

   vx = round(vert_xi)
   ve = round(vert_eta)

   if ~ array_equal(vert_xi, vx) then message, 'Non-integral vertices'
   if ~ array_equal(vert_eta, ve) then message, 'Non-integral vertices'

   if min(vx) lt 0 then message, 'Slice is out of bounds'
   if max(vx) gt dim_rho[0]-2 then message, 'Slice is out of bounds'

   if min(ve) lt 0 then message, 'Slice is out of bounds'
   if max(ve) gt dim_rho[1]-2 then message, 'Slice is out of bounds'

   dx = mgh_diff(vx)
   de = mgh_diff(ve)

   if max(dx*de) gt 0 then message, 'X & Y vertices both vary over the same segment'
   if min(abs(dx)+abs(de)) eq 0 then message, 'X & Y vertices both unchanged over the same segment'

   ;; Identify all the psi points spanned by the P-slice segments. A bit tricky, this.

   px = list()
   pe = list()

   px->Add, [vx[0]]
   pe->Add, [ve[0]]

   for i=0,n_vert-2 do begin
      if abs(dx[i]) gt 0 then begin
         px->Add, ([vx[i]:vx[i+1]])[1:*]
         pe->Add, replicate(ve[i], abs(dx[i]))
      endif else begin
         px->Add, replicate(vx[i], abs(de[i]))
         pe->Add, ([ve[i]:ve[i+1]])[1:*]
      endelse
   endfor

   mgh_undefine, vx, ve, dx, de

   point_xi = px->ToArray(DIMENSION=1)
   point_eta = pe->ToArray(DIMENSION=1)

   mgh_undefine, px, pe

   ;; Various properties now hold by construction and should not need to be checked:
   ;;
   ;;  - The point_xi and point_eta arrays have the same number of elements
   ;;  - Consecutive Pslice points never coincide
   ;;  - Consecutive Pslice points always shift by -1 or +1 in xi or eta, but not both,
   ;;    i.e. the segments between points always have length 1.
   ;;
   ;; The following code calculates the direction of each segment and relies on these
   ;; properties. Directions follow the GetTransportBox convention:
   ;;   0: positive xi, i.e. parallel with the southern boundary
   ;;   1: positive eta,  i.e. parallel with the eastern boundary
   ;;   2: negative xi, i.e. parallel with the northern boundary
   ;;   3: negative eta,  i.e. parallel with the western boundary

   n_point = n_elements(point_xi)

   direction = replicate(255B, n_point-1)

   dx = mgh_diff(point_xi)
   de = mgh_diff(point_eta)

   direction[where(dx gt 0, /NULL)] = 0
   direction[where(de gt 0, /NULL)] = 1
   direction[where(dx lt 0, /NULL)] = 2
   direction[where(de lt 0, /NULL)] = 3

   if max(direction) gt 3 then message, 'Uh oh, a segment has not been assigned a direction'

   ;; On the rho grid, locate the psi points spanned by the Pslice and the segment centres.

   point_xi_rho = point_xi + 0.5D0
   point_eta_rho = point_eta + 0.5D0

   centre_xi_rho = mgh_stagger(point_xi_rho, DELTA=-1)
   centre_eta_rho = mgh_stagger(point_eta_rho, DELTA=-1)

   ;; Interpolate various variables to the Pslice points and centres

   h = self->VarGet('h')
   point_h = interpolate(h, point_xi_rho, point_eta_rho)
   centre_h = interpolate(h, centre_xi_rho, centre_eta_rho)
   mgh_undefine, h

   if lonlat then begin
      lon = self->VarGet('lon_rho')
      lat = self->VarGet('lat_rho')
      point_lon = interpolate(lon, point_xi_rho, point_eta_rho)
      point_lat = interpolate(lat, point_xi_rho, point_eta_rho)
      mgh_undefine, lon, lat
   endif else begin
      x = self->VarGet('x_rho')
      y = self->VarGet('y_rho')
      point_x = interpolate(x, point_xi_rho, point_eta_rho)
      point_y = interpolate(y, point_xi_rho, point_eta_rho)
      mgh_undefine, x, y
   endelse

   pm = self->VarGet('pm')
   pn = self->VarGet('pn')
   centre_pm = interpolate(pm, centre_xi_rho, centre_eta_rho)
   centre_pn = interpolate(pn, centre_xi_rho, centre_eta_rho)
   mgh_undefine, pm, pn

   ;; For each segment, determine whether it has land on either side

   has_mask = self->HasVar('mask_rho')

   if has_mask then begin
      mask = self->VarGet('mask_rho')
      centre_mask = bytarr(n_point-1)
      for i=0,n_point-2 do begin
         xx = centre_xi_rho[i]
         ee = centre_eta_rho[i]
         if direction[i] mod 2 eq 0 then begin
            centre_mask[i] = mask[round(xx),floor(ee)] and mask[round(xx),ceil(ee)]
         endif else begin
            centre_mask[i] = mask[floor(xx),round(ee)] and mask[ceil(xx),round(ee)]
         endelse
      endfor
   endif

   ;; Calculate arc distance along the Pslice. If the Pslice zigzags through the
   ;; grid near a diagonal, this will be much larger than the straight line distance

   point_arc = dblarr(n_point)
   for i=0,n_point-2 do begin
      pp = (direction[i] mod 2 eq 0) ? centre_pm[i] : centre_pn[i]
      point_arc[i+1] = point_arc[i] + 1.0D/pp
   endfor

   ;; Compile and return result

   result = dictionary()

   result.lonlat = lonlat
   result.vert_xi = vert_xi
   result.vert_eta = vert_eta
   result.n_point = n_point
   result.n_segment = n_segment
   result.direction = direction
   result.point_xi = point_xi
   result.point_eta = point_eta
   result.point_xi_rho = point_xi_rho
   result.point_eta_rho = point_eta_rho
   result.centre_xi_rho = centre_xi_rho
   result.centre_eta_rho = centre_eta_rho
   result.point_h = point_h
   result.centre_h = centre_h
   if lonlat then begin
      result.point_lon = point_lon
      result.point_lat = point_lat
   endif else begin
      result.point_x = point_x
      result.point_y = point_y
   endelse
   result.centre_mask = centre_mask
   result.point_arc = point_arc

   return, result->ToStruct(/RECURSIVE, /NO_COPY)

end

; MGHromsHistory::PsliceZ
;
function MGHromsHistory::PsliceZ, variable, $
   ALONG_RANGE=along_range, BATH=bath, DIRECTION=direction, $
   GRID=grid, INDEX=index, LONLAT=lonlat, ZETA=zeta

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Get grid info if necessary.

   self->PsliceDefGrid, $
      DIRECTION=direction, INDEX=index, ALONG_RANGE=along_range, $
      GRID=grid, LONLAT=lonlat

   ;; Get the variable's vertical dimension

   dims = self->VarDims(variable)
   if strlen(dims.vertical) eq 0 then $
      message, 'Variable has no vertical dimension'

   ;; Get s-coordinate data

   scoord = self->GetScoord(dims.vertical)

   ;; Handy constants

   n_arc = n_elements(grid.arc)
   n_sc = n_elements(scoord.s)

   ;; Set defaults

   if n_elements(bath) eq 0 then bath = grid.h

   if n_elements(zeta) eq 0 then zeta = fltarr(n_arc)

   ;; Create result array & load data into it one column at a time.

   result = fltarr(n_sc, n_arc)

   cs = mgh_roms_s_to_cs(scoord.s, $
      THETA_S=scoord.theta_s, THETA_B=scoord.theta_b, $
      VSTRETCH=scoord.vstretch)

   for i=0,n_arc-1 do $
      result[0,i] = mgh_roms_s_to_z(scoord.s, bath[i], $
      ZETA=zeta[i], HC=scoord.hc, CS=cs, $
      VTRANSFORM=scoord.vtransform)

   return, transpose(result)

end

; MGHromsHistory::RsliceData
;
function MGHromsHistory::RsliceData, var

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   message, 'This function has not been fully implemented yet'

   ;; Check for problems with inputs

   if n_elements(var) ne 1 then $
      message, 'The name of a variable must be supplied'

   if ~ (isa(var, 'STRING') || isa(var, 'STRUCT'))  then $
      message, 'A variable identifier must be supplied'

   ;; Set defaults

   if n_elements(mask_value) eq 0 then mask_value = !values.f_nan

   if n_elements(depths) gt 0 || n_elements(sigmas) gt 0 then begin
      if n_elements(use_bath) eq 0 then use_bath = self->HasVar('bath')
      if n_elements(use_zeta) eq 0 then use_zeta = self->HasVar('zeta')
   endif

   ;; If no grid information is supplied, then get it. Otherwise check
   ;; grid is consistent with the current variable name & function
   ;; arguments.

   if n_elements(grid) gt 0 then begin
      dims = self->VarDims(var)
      if ~ array_equal(dims.horizontal, grid.dims.horizontal) gt 0 then $
         message, 'Horizontal dimensions of variable do not match GRID data'
      if dims.bed ne grid.dims.bed then $
         message, 'Vertical dimension of variable does not match GRID data'
      if dims.vertical ne grid.dims.vertical then $
         message, 'Vertical dimension of variable does not match GRID data'
      if dims.time ne grid.dims.time then $
         message, 'Time dimension of variable does not match GRID data'
      mgh_undefine, dims
      if n_elements(xi_range) gt 0 && ~ array_equal(xi_range,grid.xi_range) then $
         message, 'XI_RANGE does not match grid data'
      if n_elements(eta_range) gt 0 && ~ array_equal(eta_range,grid.eta_range) then $
         message, 'ETA_RANGE does not match grid data'
   endif else begin
      grid = self->HsliceGrid(var, ETA_RANGE=eta_range, LONLAT=lonlat, XI_RANGE=xi_range)
   endelse

   ;; Check consistency of function arguments with variable dimensions

   if grid.dims.vertical then begin
      n_key = (n_elements(depths) gt 0) + (n_elements(levels) gt 0) + $
         (n_elements(sigmas) gt 0)
      if n_key gt 1 then $
         message, 'The DEPTHS, LEVELS & SIGMAS keywords cannot be used together'
   endif else begin
      fmt = '(%"The %s keyword is not required or allowed when the ' + $
         'variable %s has no vertical dimension")'
      if n_elements(levels) gt 0 then message, string(FORMAT=fmt, 'LEVELS', var)
      if n_elements(depths) gt 0 then message, string(FORMAT=fmt, 'DEPTHS', var)
      if n_elements(sigmas) gt 0 then message, string(FORMAT=fmt, 'SIGMAS', var)
   endelse

   if grid.dims.bed then begin
      ;; Actually, I haven't thought of anything to test here
   endif else begin
      fmt = '(%"The %s keyword is not required or allowed when the ' + $
         'variable %s has no vertical dimension")'
      if n_elements(layers) gt 0 then message, string(FORMAT=fmt, 'LAYERS', var)
   endelse

   if grid.dims.time then begin
      if n_elements(record) eq 0 then record = 0
      if record lt 0 then $
         record = self->DimInfo(grid.dims.time, /DIMSIZE) + record
   endif else begin
      if n_elements(record) gt 0 then begin
         message, 'The RECORD keyword is not required or allowed when ' + $
            'the variable '+var+' has no time dimension'
      endif
      if keyword_set(use_bath) then begin
         message, 'The USE_BATH option may not be activated when ' + $
            'the variable '+var+' has no time dimension'
      endif
      if keyword_set(use_zeta) then begin
         message, 'The USE_ZETA option may not be activated when ' + $
            'the variable '+var+' has no time dimension'
      endif
   endelse

   ;; Abbreviations for xi and eta range

   xra0 = grid.xi_range[0]
   xra1 = grid.xi_range[1]
   xran = xra1-xra0+1
   era0 = grid.eta_range[0]
   era1 = grid.eta_range[1]
   eran = era1-era0+1

   ;; If vertical interpolation is required, get bath & zeta at the
   ;; variable location

   if n_elements(depths) gt 0 || n_elements(sigmas) gt 0 then begin

      if use_bath || use_zeta then begin
         offset = [xra0,era0,record]
         count = [xran,eran,1]
         delta = [0,0]
         hdims = strjoin(grid.dims.horizontal, ' ')
         case !true of
            strmatch(hdims, 'xi_rho eta_rho'):
            strmatch(hdims, 'xi_u eta_u') || strmatch(hdims, 'xi_u eta_rho'): begin
               count += [1,0,0]
               delta -= [1,0]
            end
            strmatch(hdims, 'xi_v eta_v') || strmatch(hdims, 'xi_rho eta_v'): begin
               count += [0,1,0]
               delta -= [0,1]
            end
            strmatch(hdims, 'xi_psi eta_psi'): begin
               count += [1,1,0]
               delta -= [1,1]
            end
         endcase
      endif

      if use_bath then begin
         bath = reform(self->VarGet('bath', OFFSET=offset, COUNT=count))
         bath = mgh_stagger(bath, DELTA=delta)
      endif else begin
         if size(grid.h, /N_DIMENSIONS) ne 2 then $
            message, 'Static bathymetry data needed but not found'
         bath = grid.h
      endelse

      if use_zeta then begin
         zeta = reform(self->VarGet('zeta', OFFSET=offset, COUNT=count))
         zeta = mgh_stagger(zeta, DELTA=delta)
      endif else begin
         zeta = fltarr(xran, eran)
      endelse

   endif

   ;; Read & unpack variable, interpolating if necessary

   case !true of

      n_elements(depths) gt 0: begin

         ;; DEPTHS keyword has been set so interpolate to constant-z levels

         n_slice = n_elements(depths)

         ;; DEPTH specified, so get 3D data & interpolate vertically.
         ;; Since vertical interpolation is required we make a result
         ;; array of floating point type immediately.

         result = replicate(!values.f_nan, xran, eran, n_slice)

         ;; Build up OFFSET & COUNT vectors for the netCDF get
         ;; operation

         offset = [xra0,era0,0]
         count = [xran,eran,grid.n_level]

         if strlen(grid.dims.time) gt 0 then begin
            offset = [offset, record]
            count = [count, 1]
         endif

         ;; Get 3D data.

         var3d = self->VarGet(var, OFFSET=offset, COUNT=count, /AUTOSCALE)

         ;; Loop horizontally thru domain interpolating to z
         ;; Clip s to avoid slightly, out-of-bounds values generated
         ;; by netCDF file packing.

         cs = mgh_roms_s_to_cs(grid.s, $
            THETA_S=grid.theta_s, THETA_B=grid.theta_b, $
            VSTRETCH=grid.vstretch)

         for i=0,xran-1 do begin
            for j=0,eran-1 do begin
               zz = mgh_roms_s_to_z(grid.s, bath[i,j], $
                  CS=cs, ZETA=zeta[i,j], HC=grid.hc, $
                  VTRANSFORM=grid.vtransform)
               varss = reform(var3d[i,j,*])
               if min(finite(varss)) gt 0 then begin
                  varzz = interpol(varss, zz, -float(depths), /SPLINE)
                  varzz[where(depths gt bath[i,j], /NULL)] = !values.f_nan
                  result[i,j,*] = varzz
               endif
            endfor
         endfor

         !null = !null

      end

      n_elements(sigmas) gt 0: begin

         ;; SIGMAS keyword has been set, so interpolate to constant-sigma levels

         n_slice = n_elements(sigmas)

         ;; SIGMAS specified, so get 3D data & interpolate vertically.
         ;; Since vertical interpolation is required we make a result
         ;; array of floating point type immediately.

         result = fltarr(xran, eran, n_slice)

         ;; Build up OFFSET & COUNT vectors for the netCDF get
         ;; operation

         offset = [xra0,era0,0]
         count = [xran,eran,n_elements(grid.s)]

         if strlen(grid.dims.time) gt 0 then begin
            offset = [offset, record]
            count = [count, 1]
         endif

         ;; Get 3D data.

         var3d = self->VarGet(var, OFFSET=offset, COUNT=count, /AUTOSCALE)

         ;; Loop horizontally thru domain interpolating to specified sigma:
         ;; Clip s to avoid slightly out-of-bounds values generated
         ;; by netCDF file packing.

         cs = mgh_roms_s_to_cs(grid.s, $
            THETA_S=grid.theta_s, THETA_B=grid.theta_b, $
            VSTRETCH=grid.vstretch)

         for i=0,xran-1 do begin
            for j=0,eran-1 do begin
               zz = mgh_roms_s_to_z(grid.s, bath[i,j], $
                  ZETA=zeta[i,j], HC=grid.hc, CS=cs, $
                  VTRANSFORM=grid.vtransform)
               zsij = sigmas*zeta[i,j]-(1-sigmas)*bath[i,j]
               varss = interpol(var3d[i,j,*], zz, zsij, /SPLINE)
               result[i,j,*] = varss
               !null = check_math()
            endfor
         endfor

      end

      else: begin

         ;; Neither DEPTHS nor SIGMAS specified, so vertical interpolation
         ;; is not necessary

         case 1B of
            strlen(grid.dims.vertical) gt 0: begin
               slice = n_elements(levels) gt 0 ? levels : grid.n_level-1
               n_slice = n_elements(slice)
            end
            strlen(grid.dims.bed) gt 0: begin
               ;; Note that the bed layer index increases downward
               slice = n_elements(layers) gt 0 ? layers : 0
               n_slice = n_elements(slice)
            end
            else: begin
               n_slice = 1
            end
         endcase

         result = fltarr(xran, eran, n_slice)

         for k=0,n_slice-1 do begin

            offset = [xra0,era0]
            count = [xran,eran]

            if grid.dims.vertical || grid.dims.bed then begin
               offset = [offset, slice[k]]
               count = [count, 1]
            endif

            if grid.dims.time then begin
               offset = [offset, record]
               count = [count, 1]
            endif

            var2d = self->VarGet(var, OFFSET=offset, COUNT=count, /AUTOSCALE)

            result[*,*,k] = temporary(var2d)

         endfor

      end

   endcase

   ;; Reset masked values

   if size(grid.mask, /N_DIMENSIONS) eq 2 then begin

      l_masked = where(grid.mask lt 0.5, n_masked)

      if n_masked gt 0 then begin
         for k=0,n_slice-1 do begin
            r = result[*,*,k]
            r[l_masked] = mask_value
            result[0,0,k] = r
         endfor
      endif

   end

   return, result

end

; MGHromsHistory::RsliceGrid
;
function MGHromsHistory::RsliceGrid

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   message, 'This function has not been fully implemented yet'

   ;; Horizontal grid specification is now complete.

   result = dictionary()

   result.dims = dims
   result.lon = self->VarGet('lon')
   result.lat = self->VarGet('lat')
   result.mask = mask
   result.angle = angle

   case !true of

      strlen(dims.vertical) gt 0: begin

         ;; For a variable with a vertical dimension we return
         ;; s-coordinate data

         ;; Retrieve & interpolate static bathymetry (if available)

         if self->HasVar('h') then begin
            h = self->VarGet('h', OFFSET=offset, COUNT=count)
            result.h = mgh_stagger(h, DELTA=delta)
         endif

         ;; Retrieve s-coordinate parameters and s values for the variable. Note that
         ;; there is similar, but not identical, code in the GetScoord method.

         result.theta_s = self->HasVar('theta_s') ? self->VarGet('theta_s') : self->AttGet('theta_s', /GLOBAL)
         result.theta_b = self->HasVar('theta_b') ? self->VarGet('theta_b') : self->AttGet('theta_b', /GLOBAL)
         result.hc = self->HasVar('hc') ? self->VarGet('hc') : self->AttGet('hc', /GLOBAL)
         result.vstretch = self->HasVar('Vstretching') ? self->VarGet('Vstretching'): 1
         result.vtransform = self->HasVar('Vtransform') ? self->VarGet('Vtransform'): 1

         ;; The following code is complicated by a couple of changes in ROMS
         ;; output files:
         ;;  - In older files, the bottom value of s_w is omitted.
         ;;  - In ROMS 2.1 and earlier, the names of the s-coordinate variables did
         ;;    not match the corresponding dimensions--heavens knows why. In ROMS 2.2
         ;;    this was fixed.

         if self->HasVar(dims.vertical) then begin
            s = (self->VarGet(dims.vertical) < 0) > (-1)
         endif else begin
            n_s_rho = self->DimInfo('s_rho', /DIMSIZE)
            case dims.vertical of
               's_rho': begin
                  if self->HasVar('sc_r') then begin
                     s =  self->VarGet('sc_r')
                  endif else begin
                     s = mgh_stagger(mgh_range(-1, 0, STRIDE=1.0D/n_s_rho), DELTA=-1)
                  endelse
               end
               's_w': begin
                  if self->HasVar('sc_w') then begin
                     s = (self->VarGet('sc_w') < 0) > (-1)
                  endif else begin
                     s = mgh_range(-1, 0, STRIDE=1.0D/n_s_rho)
                  endelse
                  if n_elements(s) eq n_s_rho then s = [-1,s]
               end
            endcase
         endelse

         result.s = s
         result.n_level = n_elements(s)

      end

      strlen(dims.bed) gt 0: begin

         ;; For a variable with a bed-layer dimension we return
         ;; the number of layers

         result.n_layer = self->DimInfo(dims.bed, /DIMSIZE)

      end

      else:

   endcase

   return, result->ToStruct(/RECURSIVE, /NO_COPY)

end

; MGHromsHistory::TimeVarName
;
function MGHromsHistory::TimeVarName, var

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

; ** X-slice and X-transect methods *******************************************

; MGHromsHistory::XsliceData
;
function MGHromsHistory::XsliceData, variable, $
     GRID=grid, MASK_VALUE=mask_value, $
     RECORD=record, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Get dimensions of variable

   dims = self->VarDims(variable)

   ;; Get grid info if necessary

   if n_elements(grid) eq 0 then grid = self->XsliceGrid(_STRICT_EXTRA=extra)

   ;; Some handy constants

   has_vert = strlen(dims.vertical) gt 0
   has_time = strlen(dims.time) gt 0

   if ~ has_vert then $
      message, 'Variable has no vertical dimension'

   ;; Set defaults

   if n_elements(mask_value) eq 0 then mask_value = !values.f_nan

   ;; Check consistency of function arguments with variable dimensions.

   if has_time then begin
      if n_elements(record) eq 0 then record = 0
      if record lt 0 then $
         record = self->DimInfo(dims.time, /DIMSIZE) + record
   endif else begin
      if n_elements(record) gt 0 then begin
         message, 'The RECORD keyword is not required or allowed when ' + $
            'the variable '+variable+' has no time dimension'
      endif
   endelse

   ;; Get s-coordinate data

   scoord = self->GetScoord(dims.vertical)

   ;; Define a rectangular region spanning the slice.

   case strjoin(dims.horizontal,' ') of
      'xi_rho eta_rho': begin
         xr = [floor(min(grid.xi)),ceil(max(grid.xi))]
         er = [floor(min(grid.eta)),ceil(max(grid.eta))]
      end
      else: begin
         message, 'Sorry, I don''t do non-rho variables'
      end
   endcase

   xr0 = xr[0]  &  xrn = xr[1]-xr[0]+1
   er0 = er[0]  &  ern = er[1]-er[0]+1

   ;; Specify parameters for getting data

   offset = [xr0,er0,0]
   count = [xrn,ern,0]

   if has_time gt 0 then begin
      offset = [offset, record]
      count = [count, 1]
   endif

   ;; Get data and interpolate onto transect

   result = fltarr(grid.n_points, scoord.n)

   data = self->VarGet(variable, OFFSET=offset, COUNT=count, /AUTOSCALE)
   for k=0,scoord.n-1 do begin
      result[*,k] = interpolate(data[*,*,k], grid.xi-xr0, grid.eta-er0)
   endfor

   return, result

end

; MGHromsHistory::XsliceGrid
;
function MGHromsHistory::XsliceGrid, $
     INTERVAL=interval, LONLAT=lonlat, N_INTERMEDIATE=n_intermediate, VERTX=vertx, VERTY=verty, TYPE=type

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(type) eq 0 then type = 'linear'

   ;; Are the grid locations available in (lon,lat) coordinates?

   if n_elements(lonlat) eq 0 then $
      lonlat = self->HasVar('lon_rho') && self->HasVar('lat_rho')

   ;; Retrieve (x,y) or (lon,lat) positions for RHO grid (these
   ;; include exterior points)

   x_rho = lonlat ? self->VarGet('lon_rho') : self->VarGet('x_rho')
   y_rho = lonlat ? self->VarGet('lat_rho') : self->VarGet('y_rho')

   ;; Get dimensions of RHO grid

   dim_rho = self->DimRho()

   ;; Process the supplied VERTX & VERTY data according to the Xslice type

   case strlowcase(type) of

      'linear': begin

         ;; A linear Xslice is a straight line (a great circle if LONLAT = 1)
         ;; defined by two end points.

         ;; The default vertices are at the centre of the west (minimum xi) and
         ;; east (maximum xi) normal-velocity boundaries.

         if n_elements(vertx) eq 0 then $
            vertx = interpolate(x_rho, [0.5,dim_rho[0]-1.5], replicate(0.5*(dim_rho[1]-1), 2))
         if n_elements(verty) eq 0 then $
            verty = interpolate(y_rho, [0.5,dim_rho[0]-1.5], replicate(0.5*(dim_rho[1]-1), 2))

         ;; Check arguments

         if n_elements(vertx) ne 2 then $
            message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'vertx'
         if n_elements(verty) ne 2 then $
            message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'verty'

         if n_elements(interval) gt 0 then $
            message, 'The INTERVAL keyword must not be specified for a linear Xslice'

         if n_elements(n_intermediate) eq 0 then n_intermediate = 100

         ;; Specify the path

         if lonlat then begin
            xy = map_2points(vertx[0], verty[0], vertx[1], verty[1], NPATH=n_intermediate+2)
            x = reform(xy[0,*])
            y = reform(xy[1,*])
            mgh_undefine, xy
         endif else begin
            x = mgh_range(vertx, N_ELEMENTS=n_intermediate+2)
            y = mgh_range(verty, N_ELEMENTS=n_intermediate+2)
         endelse

      end

      'spline': begin

         ;; A spline Xslice is generated by the spline_p procedure and  defined by two
         ;; or more end points.

         ;; The default vertices are at the centre of the west (minimum xi) and
         ;; east (maximum xi) normal-velocity boundaries, as with the linear Xslice

         if n_elements(vertx) eq 0 then $
            vertx = interpolate(x_rho, [0.5,dim_rho[0]-1.5], replicate(0.5*(dim_rho[1]-1), 2))

         if n_elements(verty) eq 0 then $
            verty = interpolate(y_rho, [0.5,dim_rho[0]-1.5], replicate(0.5*(dim_rho[1]-1), 2))

         if n_elements(verty) ne n_elements(vertx) then $
            message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'vertx', 'verty'

         if n_elements(n_intermediate) gt 0 then $
            message, 'The N_INTERMEDIaTE keyword must not be specified for a spline Xslice'

         ;; Specify the path

         spline_p, vertx, verty, x, y, INTERVAL=interval

      end

      'piecewise': begin

         ;; A piecewise Xslice is a series of straight lines (great circles if LONLAT = 1), each
         ;; defined by two end points.

         ;; The default vertices are at the centre of the west (minimum xi) and
         ;; east (maximum xi) normal-velocity boundaries.

         if n_elements(vertx) eq 0 then $
            vertx = interpolate(x_rho, [0.5,dim_rho[0]-1.5], replicate(0.5*(dim_rho[1]-1), 2))
         if n_elements(verty) eq 0 then $
            verty = interpolate(y_rho, [0.5,dim_rho[0]-1.5], replicate(0.5*(dim_rho[1]-1), 2))

         ;; Check arguments

         if n_elements(interval) gt 0 then $
            message, 'The INTERVAL keyword must not be specified for a piecewise linear Xslice'

         if n_elements(verty) ne n_elements(vertx) then $
            message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'vertx', 'verty'

         n_vert = n_elements(vertx)
         if n_vert lt 2 then $
            message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'vertx'

         if n_elements(n_intermediate) eq 0 then n_intermediate = round(100/float(n_vert-1))

         n_points = n_vert + (n_vert-1)*n_intermediate

         x = dblarr(n_points)
         y = dblarr(n_points)

         ;; Specify the path

         for i=0,n_vert-2 do begin
            n0 = i*(n_intermediate+1)
            if lonlat then begin
               xy = map_2points(vertx[i], verty[i], vertx[i+1], verty[i+1], NPATH=n_intermediate+2)
               x[n0:n0+n_intermediate+1] = reform(xy[0,*])
               y[n0:n0+n_intermediate+1] = reform(xy[1,*])
               mgh_undefine, xy
            endif else begin
               x[n0:n0+n_intermediate+1] = mgh_range(vertx[i:i+1], N_ELEMENTS=n_intermediate+2)
               y[n0:n0+n_intermediate+1] = mgh_range(verty[i:i+1], N_ELEMENTS=n_intermediate+2)
            endelse
         endfor


      end

   endcase

   ;; For each point on the slice get (xi,eta) positions relative
   ;; to the RHO grid.

   loc = self->LocateXY(x, y, LONLAT=lonlat)

   xi  = reform(loc[0,*])
   eta = reform(loc[1,*])

   n_points = n_elements(xi)

   ;; Calculate spacing and angle of the path (spacing is in metres;
   ;; angle has same convention as the ROMS angle variable, i.e. positive
   ;; is anti-clockwise from E in radians). If lonlat is defined, assume a
   ;; spherical Earth.

   if lonlat then begin
      radius = 6371.D3
      dist = dblarr(n_points-1)
      angle = dblarr(n_points-1)
      for i=0,n_points-2 do begin
         ra = map_2points(x[i], y[i], x[i+1], y[i+1])
         dist[i] = ra[0]*radius*!dpi/180
         angle[i] = (90-ra[1])*!dpi/180
      endfor
   endif else begin
      cc = mgh_diff(x + !const.i*y)
      dist = abs(cc)
      angle = atan(cc, /PHASE)
      mgh_undefine, cc
   endelse

   ;; Calculate distance along the path

   arc = fltarr(n_points)
   for i=0,n_points-2 do arc[i+1] = arc[i] + dist[i]

   ;; We want to interpolate bathymetry, grid angle, and mask info to
   ;; the slice. First define a rectangular region spanning the slice.

   xr = [floor(min(xi)),ceil(max(xi))]
   er = [floor(min(eta)),ceil(max(eta))]

   xr0 = xr[0]  &  xrn = xr[1]-xr[0]+1
   er0 = er[0]  &  ern = er[1]-er[0]+1

   ;; Retrieve and interpolate bathymetry

   h = self->VarGet('h', OFFSET=[xr0,er0], COUNT=[xrn,ern])
   h = interpolate(temporary(h), xi-xr0, eta-er0)

   ;; Retrieve the model grid's rho mask. Use nearest-neighbour
   ;; interpolation, so that the resulting value matches the usual
   ;; square-block interpretation of the mask.

   mask = 0
   if self->HasVar('mask_rho') then begin
      mask = self->VarGet('mask_rho', OFFSET=[xr0,er0], COUNT=[xrn,ern])
      mask = mask[round(xi-xr0),round(eta-er0)]
   endif

   ;; Adjust angle so it is relative to ROMS grid. Keep
   ;; angle relative to geographic coordinates as true_angle

   true_angle = angle

   if self->HasVar('angle') then begin
      roms_angle = self->VarGet('angle', OFFSET=[xr0,er0], COUNT=[xrn,ern])
      angle = angle - interpolate(temporary(roms_angle), mgh_stagger(xi, DELTA=-1)-xr0, mgh_stagger(eta, DELTA=-1)-er0)
   endif

   ;; Return result structure

   return, {lonlat: lonlat, n_points: n_points, x: x, y: y, arc: arc, $
      xi: xi, eta: eta, angle: angle, true_angle: true_angle, h: h, mask: mask}

end

; MGHromsHistory::XtranData
;
function MGHromsHistory::XtranData, variable, $
     DEPTHS=depths, LEVELS=levels, SIGMAS=sigmas, $
     GRID=grid, RECORD=record, USE_BATH=use_bath, USE_ZETA=use_zeta, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Set defaults

   if n_elements(use_bath) eq 0 then use_bath = !false
   if n_elements(use_zeta) eq 0 then use_zeta = !false

   ;; Get dimensions of variable

   dims = self->VarDims(variable)

   ;; Get grid info if necessary

   if n_elements(grid) eq 0 then grid = self->XsliceGrid(_STRICT_EXTRA=extra)

   ;; Some handy switches

   has_vert = strlen(dims.vertical) gt 0
   has_time = strlen(dims.time) gt 0

   ;; Check consistency of function arguments with variable dimensions.

   if has_vert then begin
      n_key = (n_elements(depths) gt 0) + (n_elements(levels) gt 0) + $
         (n_elements(sigmas) gt 0)
      if n_key gt 1 then $
         message, 'The DEPTHS, LEVELS & SIGMAS keywords cannot be used together'
   endif else begin
      if n_elements(depths) gt 0 then begin
         message, 'The DEPTHS keyword is not required or allowed when ' + $
            'the variable '+variable+' has no vertical dimension'
      endif
      if n_elements(levels) gt 0 then begin
         message, 'The LEVELS keyword is not required or allowed when ' + $
            'the variable '+variable+' has no vertical dimension'
      endif
      if n_elements(levels) gt 0 then begin
         message, 'The SIGMAS keyword is not required or allowed when ' + $
            'the variable '+variable+' has no vertical dimension'
      endif
   endelse

   if has_time then begin
      if n_elements(record) eq 0 then record = 0
      if record lt 0 then $
         record = self->DimInfo(dims.time, /DIMSIZE) + record
   endif else begin
      if n_elements(record) gt 0 then begin
         message, 'The RECORD keyword is not required or allowed when ' + $
            'the variable '+variable+' has no time dimension'
      endif
      if use_bath then begin
         message, 'The USE_BATH option may not be activated when ' + $
            'the variable '+variable+' has no time dimension'
      endif
      if use_zeta then begin
         message, 'The USE_ZETA option may not be activated when ' + $
            'the variable '+variable+' has no time dimension'
      endif
   endelse

   ;; If variable has a vertical dimension, get s-coordinate data

   if has_vert then scoord = self->GetScoord(dims.vertical)

   ;; Establish levels/depths at which data are required

   get_depths = n_elements(depths) gt 0
   get_sigmas = n_elements(sigmas) gt 0
   get_levels = n_elements(levels) gt 0

   case !true of
      get_depths: n_tran = n_elements(depths)
      get_levels: n_tran = n_elements(levels)
      get_sigmas: n_tran = n_elements(sigmas)
      else: begin
         n_tran = 1
         if has_vert then levels = scoord.n-1
      endelse
   endcase

   ;; Create array to hold result.

   result = fltarr(grid.n_points, n_tran)

   ;; Define a rectangular region spanning the slice.

   case strjoin(dims.horizontal,' ') of
      'xi_rho eta_rho': begin
         xr = [floor(min(grid.xi)),ceil(max(grid.xi))]
         er = [floor(min(grid.eta)),ceil(max(grid.eta))]
      end
      else: begin
         message, 'Sorry, I don''t do non-rho variables'
      end
   endcase

   xr0 = xr[0]  &  xrn = xr[1]-xr[0]+1
   er0 = er[0]  &  ern = er[1]-er[0]+1

   offset_h = [xr0,er0]
   count_h = [xrn,ern]

   ;; Read variable

   if get_depths || get_sigmas then begin

      ;; The DEPTHS or SIGMAS keyword has been specified, so get a
      ;; complete slice and interpolate vertically

      message, "Sorry, I don't do DEPTHS or SIGMAS yet!"

      ;; Build up offset & count vectors for VarGet

      offset = [offset_h,0]
      count = [count_h,n_elements(grid.s)]
      if has_time then begin
         offset = [offset, record]
         count = [count, 1]
      endif

      ;; Get a slice

      vslice = reform(self->VarGet(variable, OFFSET=offset, COUNT=count))

      ;; Get zeta & bathymetry @ variable location

      if keyword_set(use_bath) then begin
         message, "Sorry I don't do USE_bath=1 yet!"
      endif else begin
         bath = grid.h
      endelse

      if keyword_set(use_zeta) then begin
         message, "Sorry I don't do USE_ZETA=1 yet!"
      endif else begin
         zeta = fltarr(n_points)
      endelse

      ;; Loop horizontally thru domain interpolating to all depths

      cs = mgh_roms_s_to_cs(grid.s, $
         THETA_S=grid.theta_s, THETA_B=grid.theta_b, $
         VSTRETCH=grid.vstretch)

      for i=0,aran-1 do begin
         zz = mgh_roms_s_to_z(grid.s, bath[i], $
            ZETA=zeta[i], HC=grid.hc, CS=cs, $
            VTRANSFORM=grid.vtransform)
         if get_sigmas then begin
            zs = sigmas*zeta[i]-(1-sigmas)*bath[i]
            data = interpol(vslice[i,*], zz, zs, /QUADRATIC)
         endif else begin
            data = interpol(vslice[i,*], zz, -float(depths), /QUADRATIC)
            bottom = where(depths gt bath[i], count)
            if count gt 0 then data[bottom] = !values.f_nan
         endelse
         if grid.average then begin
            result[i,*] = result[i,*] + data/n_slice
         endif else begin
            result[i,*,s] = data
         endelse
      endfor

   endif else begin

      ;; DEPTHS & SIGMAS not specified, so no vertical interpolation
      ;; necessary. Get data for each level separately

      ;; Build up offset & count vectors for VarGet

      offset = offset_h
      count = count_h

      if has_vert then begin
         offset = [offset,levels]
         count = [count,1]
      endif

      if has_time gt 0 then begin
         offset = [offset, record]
         count = [count, 1]
      endif

      ;; Get data and interpolate onto transect

      result = self->VarGet(variable, OFFSET=offset, COUNT=count)
      result = interpolate(reform(result), grid.xi-xr0, grid.eta-er0)

   endelse

   return, result

end

; MGHromsHistory::XsliceUVbar
;
;   Calculate & return depth-averaged velocity normal and tangential to
;   an Xslice.  The tangential (normal) velocity is the real
;   (imaginary) part.
;
function MGHromsHistory::XsliceUVbar, $
     GRID=grid, RECORD=record, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Get grid info if necessary

   if n_elements(grid) eq 0 then grid = self->XsliceGrid(_STRICT_EXTRA=extra)

   if n_elements(record) eq 0 then record = 0

   ;; Get ubar & vbar data for a rectangular region spanning the
   ;; X-slice

   xr = [floor(min(grid.xi)),ceil(max(grid.xi))]
   er = [floor(min(grid.eta)),ceil(max(grid.eta))]

   xr0 = xr[0]  &  xrn = xr[1]-xr[0]+1
   er0 = er[0]  &  ern = er[1]-er[0]+1

   uv = self->VectorGet(['ubar','vbar'], OFFSET=[xr0,er0,record], COUNT=[xrn,ern,1])

   ;; Interpolate to slice interval mid-points

   xx = mgh_stagger(grid.xi, DELTA=[-1])-xr0
   ee = mgh_stagger(grid.eta, DELTA=-1)-er0

   uv = interpolate(temporary(uv), xx, ee)

   ;; Rotate into slice-relative coordinates and return.

   return, temporary(uv)*exp(!const.i*grid.grid_angle)*exp(-!const.i*grid.angle)

end

; MGHromsHistory::XsliceZ
;
function MGHromsHistory::XsliceZ, variable, $
     GRID=grid, BATH=bath, ZETA=zeta, $
     _REF_EXTRA=extra

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  ;; Get grid info if necessary

  if n_elements(grid) eq 0 then grid = self->XsliceGrid(_STRICT_EXTRA=extra)

  ;; Get the variable's vertical dimension

  dims = self->VarDims(variable)
  if strlen(dims.vertical) eq 0 then $
       message, 'Variable has no vertical dimension'

  ;; Get s-coordinate data

  scoord = self->GetScoord(dims.vertical)

  ;; Handy constants

  n_arc = n_elements(grid.arc)
  n_sc = n_elements(scoord.s)

  ;; Set defaults

  if n_elements(bath) eq 0 then bath = grid.h

  if n_elements(zeta) eq 0 then zeta = fltarr(n_arc)

  ;; Create result array & load data into it one column at a time.

  result = fltarr(n_sc, n_arc)

  cs = mgh_roms_s_to_cs(scoord.s, $
                        THETA_S=scoord.theta_s, THETA_B=scoord.theta_b, $
                        VSTRETCH=scoord.vstretch)

  for i=0,n_arc-1 do $
       result[0,i] = mgh_roms_s_to_z(scoord.s, bath[i], $
                                     ZETA=zeta[i], HC=scoord.hc, CS=cs, $
                                     VTRANSFORM=scoord.vtransform)

  return, transpose(result)

end

; MGHromsHistory::PatchGrid
;
function MGHromsHistory::PatchGrid, $
  LONLAT=lonlat, VERTX=vertx, VERTY=verty

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Are the grid locations available in (lon,lat) coordinates?

   if n_elements(lonlat) eq 0 then $
      lonlat = self->HasVar('lon_rho') && self->HasVar('lat_rho')

   ;; Get rho grid dimensions

   dim = self->DimRho()

   ;; Retrieve (x,y) or (lon,lat) positions for RHO grid (these
   ;; include exterior points)

   x_rho = lonlat ? self->VarGet('lon_rho') : self->VarGet('x_rho')
   y_rho = lonlat ? self->VarGet('lat_rho') : self->VarGet('y_rho')

   ;; The psi grid defines the corners of the grid cells & will be needed later.

   x_psi = mgh_stagger(x_rho, DELTA=[-1,-1])
   y_psi = mgh_stagger(y_rho, DELTA=[-1,-1])

   mgh_undefine, x_rho, y_rho

   ;; Get the land mask, if any.

   has_mask = self->HasVar('mask_rho')
   if has_mask then begin
      mask = round(self->VarGet('mask_rho'))
      mask = mask[1:-2,1:-2]
   endif

   ;; Construct a [2,n] polygon containing the vertices

   if n_elements(vertx) ne n_elements(verty) then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'vertx', 'verty'

   if size(vertx, /N_DIMENSIONS) gt 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'vertx'
   if size(verty, /N_DIMENSIONS) gt 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'verty'

   if n_elements(vertx) eq 0 then begin
      pol = mgh_perim(x_psi, y_psi)
   endif else begin
      pol = [transpose(vertx),transpose(verty)]
   endelse

   ;; Calculate the fractional area of each grid (rho) cell that is contained within the
   ;; perimeter polygon. Consider interior cells only

   frac = fltarr(dim[0:1]-2)

   for j=0,dim[1]-3 do begin
      for i=0,dim[0]-3 do begin
         if has_mask && mask[i,j] eq 0 then continue
         pc = pol
         lx = [x_psi[i,j],x_psi[i+1,j],x_psi[i+1,j+1],x_psi[i,j+1]]
         ly = [y_psi[i,j],y_psi[i+1,j],y_psi[i+1,j+1],y_psi[i,j+1]]
         pc = mgh_polyclip2(pc, mgh_line_coeff(lx[0], ly[0], lx[1], ly[1]), COUNT=count)
         if count eq 0 then continue
         pc = mgh_polyclip2(pc, mgh_line_coeff(lx[1], ly[1], lx[2], ly[2]), COUNT=count)
         if count eq 0 then continue
         pc = mgh_polyclip2(pc, mgh_line_coeff(lx[2], ly[2], lx[3], ly[3]), COUNT=count)
         if count eq 0 then continue
         pc = mgh_polyclip2(pc, mgh_line_coeff(lx[3], ly[3], lx[0], ly[0]), COUNT=count)
         if count eq 0 then continue
         cx = reform(pc[0,*])
         cy = reform(pc[1,*])
         frac[i,j] = abs(total(cx*shift(cy, -1) - cy*shift(cx, -1))) / $
            abs(total(lx*shift(ly, -1) - ly*shift(lx, -1)))
      endfor
   endfor

   ;; Establish the bounds of a rectangle enclosing the non-empty cells.
   ;; Express the corners of this rectangle as indices into the entire
   ;; rho grid

   if max(frac) eq 0 then message, 'No non-empty cells found'

   for i=0,dim[0]-3 do begin
      if max(frac[i,*]) gt 0 then break
   endfor
   xr0 = i+1
   for i=dim[0]-3,0,-1 do begin
      if max(frac[i,*]) gt 0 then break
   endfor
   xr1 = i+1

   for j=0,dim[1]-3 do begin
      if max(frac[*,j]) gt 0 then break
   endfor
   er0 = j+1
   for j=dim[1]-3,0,-1 do begin
      if max(frac[*,j]) gt 0 then break
   endfor
   er1 = j+1

   xrn = xr1-xr0+1
   ern = er1-er0+1

   ;; Trim off the cells outside the bounding rectangle

   frac = frac[xr0-1:xr1-1,er0-1:er1-1]

   ;; Return result structure

   return, {lonlat: lonlat, xi_range: [xr0,xr1], eta_range: [er0,er1], frac: frac}

end

; MGHromsHistory::VarDims

function MGHromsHistory::VarDims, var

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
      if strmatch(dim[d], 'xi_*') then $
         result.horizontal[0] = dim[d]
      if strmatch(dim[d], 'eta_*') then $
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

; MGHromsHistory::VarDimNames
;
function MGHromsHistory::VarDimNames, var, COUNT=count

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

      isa(var, 'STRING') && (var eq 'Dsbl'): begin
         result = self->MGHncSequence::VarDimNames('Hsbl', COUNT=count)
      end

      isa(var, 'STRING') && (var eq 'Dbbl'): begin
         result = self->MGHncSequence::VarDimNames('Hbbl', COUNT=count)
      end

      isa(var, 'STRING') && strmatch(var, 'bed_sum:*'): begin
         result = self->VarDimNames(strmid(var, 8), COUNT=count)
         result = result[[0,1,3]]
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

; MGHromsHistory::VarGet
;
function MGHromsHistory::VarGet, var, $
     AUTOSCALE=autoscale, COUNT=count, OFFSET=offset

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(var) eq 0 then $
      message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', var

   if n_elements(autoscale) eq 0 then autoscale = 1B

   case !true of

      isa(var, 'STRING') && (var eq 'psi'): begin
         result = self->VarGet('psi(ubar,vbar,zeta)', AUTOSCALE=autoscale, COUNT=count, OFFSET=offset)
      end

      isa(var, 'STRING') && strmatch(var, 'psi(*,*,*)'): begin
         pp = stregex(var, '(^psi\()(.+)(,)(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4,6]], ll[[2,4,6]])
         dim_rho = self->DimRho()
         my_count = n_elements(count) gt 0 ? count : [dim_rho+[-1,-1],0]
         if my_count[0] eq 0 then my_count[0] = dim_rho[0]-1
         if my_count[1] eq 0 then my_count[1] = dim_rho[1]-1
         my_offset = n_elements(offset) gt 0 ? offset : [0,0,0]
         ;; Get grid data
         pm = self->MGHncSequence::VarGet('pm', COUNT=my_count[0:1]+[1,1], OFFSET=my_offset[0:1], /AUTOSCALE)
         pn = self->MGHncSequence::VarGet('pn', COUNT=my_count[0:1]+[1,1], OFFSET=my_offset[0:1], /AUTOSCALE)
         if self->MGHncSequence::HasVar('mask_rho') then begin
            mask = self->MGHncSequence::VarGet('mask_rho', COUNT=my_count[0:1]+[1,1], OFFSET=my_offset[0:1], /AUTOSCALE)
         endif
         h = self->MGHncSequence::VarGet('h', COUNT=my_count[0:1]+[1,1], OFFSET=my_offset[0:1], /AUTOSCALE)
         h_u = mgh_stagger(h, DELTA=[-1,0])
         h_v = mgh_stagger(h, DELTA=[0,-1])
         mgh_undefine, h
         ;; Get time_varying data
         ubar = self->MGHncSequence::VarGet(vv[0], COUNT=my_count+[0,1,0], OFFSET=my_offset, /AUTOSCALE)
         vbar = self->MGHncSequence::VarGet(vv[1], COUNT=my_count+[1,0,0], OFFSET=my_offset, /AUTOSCALE)
         l_miss = where(~ finite(ubar), n_miss)
         if n_miss gt 0 then ubar[l_miss] = 0
         l_miss = where(~ finite(vbar), n_miss)
         if n_miss gt 0 then vbar[l_miss] = 0
         n_dim = size(ubar, /N_DIMENSIONS)
         n_rec = n_dim eq 3 ? (size(ubar, /DIMENSIONS))[2] : 1
         zeta = self->MGHncSequence::VarGet(vv[2], COUNT=my_count+[1,1,0], OFFSET=my_offset, /AUTOSCALE)
         result = dblarr([dim_rho[0:1]-1,n_rec])
         for r=0,n_rec-1 do begin
            zb = mgh_fill2d(zeta[*,*,r])
            zeta_u = mgh_stagger(zb, DELTA=[-1,0])
            zeta_v = mgh_stagger(zb, DELTA=[0,-1])
            mgh_undefine, zb
            result[*,*,r] = mgh_roms_psi(ubar[*,*,r]*(h_u+zeta_u), vbar[*,*,r]*(h_v+zeta_v), pm, pn, MASK=mask)
         endfor
      end

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

      isa(var, 'STRING') && (var eq 'Dsbl'): begin
         result = - self->MGHncSequence::VarGet('Hsbl', AUTOSCALE=autoscale, COUNT=count, OFFSET=offset)
      end

      isa(var, 'STRING') && (var eq 'Dbbl'): begin
         dim_rho = self->DimRho()
         my_offset = n_elements(offset) gt 0 ? offset : [0,0,0]
         my_count = n_elements(count) gt 0 ? count : [dim_rho[0:1],0]
         result = self->MGHncSequence::VarGet('Hbbl', AUTOSCALE=autoscale, COUNT=my_count, OFFSET=my_offset)
         h = self->MGHncSequence::VarGet('h', AUTOSCALE=autoscale, COUNT=my_count[0:1], OFFSET=my_offset[0:1])
         if my_count[2] ne 1 then begin
            h = rebin(reform(h, [size(h, /DIMENSIONS),1]), [size(h, /DIMENSIONS),my_count[2]])
         endif
         result += h
      end

      isa(var, 'STRING') && (var eq 'bed_total'): begin
         if n_elements(offset) gt 0 then begin
            my_offset = [offset[0:1],0,offset[2]]
         endif else begin
            my_offset = [0,0,0,0]
         endelse
         if n_elements(count) gt 0 then begin
            my_count = [count[0:1],0,count[2]]
         endif else begin
            my_count = [0,0,0,0]
         endelse
         result = self->VarGet('bed_thickness', $
            COUNT=my_count, OFFSET=my_offset)
         return, total(result, 3)
      end

      isa(var, 'STRING') && strmatch(var, 'bed_sum:*'): begin
         if n_elements(offset) gt 0 then begin
            my_offset = [offset[0:1],0,offset[2]]
         endif else begin
            my_offset = [0,0,0,0]
         endelse
         if n_elements(count) gt 0 then begin
            my_count = [count[0:1],0,count[2]]
         endif else begin
            my_count = [0,0,0,0]
         endelse
         result = self->VarGet(strmid(var, 8), $
            COUNT=my_count, OFFSET=my_offset)
         if size(result, /N_DIMENSIONS) ge 3 then $
            result = total(result, 3)
         return, result
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

pro MGHromsHistory::VarInfo, var, $
   ALL=all, ATT_NAMES=att_names, DATATYPE=datatype, DIM_NAMES=dim_names, $
   DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_ATTS=n_atts, N_DIMS=n_dims

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case !true of

      isa(var, 'STRING') && (var eq 'psi'): begin
         self->VarInfo, 'psi(ubar,vbar,zeta)', $
            ATT_NAMES=att_names, DATATYPE=datatype, DIM_NAMES=dim_names, $
            DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_ATTS=n_atts, N_DIMS=n_dims
      end

      isa(var, 'STRING') && strmatch(var, 'psi(*,*,*)'): begin
         pp = stregex(var, '(^psi\()(.+)(,)(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4,6]], ll[[2,4,6]])
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, vv[0], $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
         dim_names[0] = 'xi_psi'
         dim_names[1] = 'eta_psi'
         dimensions[0] = self->MGHncSequence::DimInfo('xi_psi', /DIMSIZE)
         dimensions[1] = self->MGHncSequence::DimInfo('eta_psi', /DIMSIZE)
      end

      isa(var, 'STRING') && strmatch(var, 'abs(*,*)'): begin
         pp = stregex(var, '(^abs\()(.+)(,)(.+)(\))', LENGTH=ll, /SUBEXPR)
         vv = strmid(var, pp[[2,4]], ll[[2,4]])
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, vv[0], $
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

      isa(var, 'STRING') && (var eq 'Dsbl'): begin
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, 'Hsbl', $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
      end

      isa(var, 'STRING') && (var eq 'Dbbl'): begin
         att_names = ''
         datatype = 'FLOAT'
         n_atts = 0
         self->VarInfo, 'Hbbl', $
            DIM_NAMES=dim_names, DIMENSIONS=dimensions, FILL_VALUE=fill_value, N_DIMS=n_dims
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
; MGHromsHistory::VectorGet
;
function MGHromsHistory::VectorGet, var, $
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

pro MGHromsHistory__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, {MGHromsHistory, inherits MGHncSequence}

end
