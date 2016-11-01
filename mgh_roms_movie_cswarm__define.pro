;+
; CLASS NAME:
;   MGH_ROMS_Movie_Cswarm
;
; PURPOSE:
;   This procedure generates and displays an animated graph showing
;   ROMS float locations on vertical slice (C-slice).
;
; CALLING SEQUENCE:
;   mgh_new, 'MGH_ROMS_Movie_Cswarm', ffile, hfile
;
; SUPERCLASS:
;   MGH_Datamator
;
; POSITIONAL PARAMETERS:
;  ffile
;    Synonym for FLOAT_FILE property.
;
;  hfile
;    Synonym for HISTORY_FILE property.
;
; PROPERTIES:
;   ALONG_RANGE
;     A 2-element integer vector specifying the grid subset to be plotted
;    in the along-slice direction. The range is specified relative to
;    the rho grid. Negative values are taken to be offsets from the
;    "northern" edge of the grid. The default is [0,-1] which is
;    equivalent to [0,dim(eta_rho)-1].
;
;  FLOAT_FILE (Init, Get)
;    A reference to a ROMS float-file object.
;
;  FLOAT_RANGE (Init)
;  FLOAT_STRIDE (Init)
;    Range and stride specifying a subset of floats to be
;    plotted. These are ignored if FLOATS is specified.
;
;  FLOATS (Init)
;    An integer vector specifying the floats to be plotted
;
;  GRAPH_PROPERTIES (Init)
;    A structure containing keywords to be passed to the graph.
;
;  HISTORY_FILE (Init, Get)
;    A reference to a ROMS history-file object.
;
;  RECORD_RANGE (Init)
;  RECORD_STRIDE (Init)
;    Range and stride specifying a subset of records to be
;    plotted. These are ignored if RECORDS is specified.
;
;  RECORDS (Init)
;    An integer vector specifying the records to be plotted
;
;  SYMBOL_PROPERTIES (Init)
;    A structure containing keywords to be passed to the plotting
;    symbol.
;
;  XAXIS_PROPERTIES (Init)
;  YAXIS_PROPERTIES (Init)
;    Structures containing keywords to be passed to the X & Y axis.
;
;  XI_RANGE (Init)
;    A 2-element integer vector specifying the grid points to be plotted
;    in the xi direction. The range is specified relative to the rho grid.
;    Negative values are taken to be offsets from the "eastern" edge of
;    the grid. The default is [0,-1] which is equivalent to [0,dim(xi_rho)-1].
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2001-03:
;     Written.
;   Mark Hadfield, 2001-07:
;     Now a subclass of MGH_Datamator. Updated for IDL 5.5.
;   Mark Hadfield, 2010-12:
;     In Init method, scalar float-position data in the Init method is
;     now forced to array.
;   Mark Hadfield, 2011-07:
;     Removed the float_destroy and history_destroy fields:
;     unnecessary with automatic gargarbage collection.
;   Mark Hadfield, 2012-04:
;     Fix bug: missing endfor statement.
;-
function MGH_ROMS_Movie_Cswarm::Init, ffile, hfile, $
     ALONG_RANGE=along_range, $
     DEPTH_RANGE=depth_range, $
     DIRECTION=direction, $
     EXAGGERATE_ZETA=exaggerate_zeta, $
     FLOAT_FILE=float_file, $
     FLOAT_RANGE=float_range, $
     FLOAT_STRIDE=float_stride, $
     FLOATS=floats, $
     GRAPH_PROPERTIES=graph_properties, $
     HISTORY_FILE=history_file, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     RECORDS=records, $
     SYMBOL_PROPERTIES=symbol_properties, $
     USE_ZETA=use_zeta, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
     _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process float-file argument

   if n_elements(float_file) eq 0 && n_elements(ffile) gt 0 then $
        float_file = ffile

   if n_elements(float_file) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'float_file'

   case size(float_file, /TNAME) of
      'STRING': begin
         self.float_file = obj_new('MGHromsFloat', float_file)
      end
      'OBJREF': begin
         self.float_file = float_file
      end
   endcase
   oflt = self.float_file

   ;; Process history-file argument

   if n_elements(history_file) eq 0 then $
        if n_elements(hfile) gt 0 then history_file = hfile

   if n_elements(history_file) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'history_file'

   case size(history_file, /TNAME) of
      'STRING': begin
         self.history_file = obj_new('MGHromsHistory', history_file)
      end
      'OBJREF': begin
         self.history_file = history_file
      end
   endcase
   ohis = self.history_file

   ;; Defaults

   if n_elements(direction) eq 0 then direction = 0

   if n_elements(exaggerate_zeta) eq 0 then exaggerate_zeta = 1.

   if n_elements(use_zeta) eq 0 then use_zeta = 0

   ;; Get number of floats and resolve float-selection arguments

   n_float = oflt->DimInfo('drifter', /DIMSIZE)

   mgh_resolve_indices, n_float, float_range, float_stride, floats

   ;; Establish variable and dimension name associated with float data

   case !true of
      oflt->HasVar('Xgrid') && oflt->HasAtt('Xgrid', 'time'): $
           time_var = oflt->AttGet('Xgrid','time')
      oflt->HasVar('Ygrid') && oflt->HasAtt('Ygrid', 'time'): $
           time_var = oflt->AttGet('Ygrid','time')
      else: $
           time_var = 'ocean_time'
   endcase

   time_dim = (oflt->VarDimNames(time_var))[0]

   ;; Get number of records in float file and resolve record-selection
   ;; arguments

   n_time = oflt->DimInfo(time_dim, /DIMSIZE)

   mgh_resolve_indices, n_time, record_range, record_stride, records

   ;; Get grid dimensions

   dim_rho = ohis->DimRho()

   ;; Set up a C-slice defining the plane onto which the float
   ;; locations will be projected.

   if n_elements(slice) eq 0 then $
        slice = round(0.5*(dim_rho[1-direction]-2))

   if n_elements(along_range) eq 0 then along_range = [0,-1]

   if slice lt 0 then $
        slice = dim_rho[1-direction] + slice

   if along_range[0] lt 0 then $
        along_range[0] = dim_rho[direction] + along_range[0]
   if along_range[1] lt 0 then $
        along_range[1] = dim_rho[direction] + along_range[1]

   n_along = along_range[1]-along_range[0]+1

   dims = {horizontal: ['xi_rho','eta_rho'], vertical:'s_w', time:''}

   grid = ohis->CsliceGrid(dims, ALONG_RANGE=along_range, DIRECTION=direction, $
                           SLICE=slice)
   self.grid = ptr_new(grid)

   ;; Calculate Z for all points on the C-slice assuming zeta=0.
   ;; Include values for the bottom (s = 0).  If USE_ZETA is set, this
   ;; will be recalculated at each frame.

   grid_z = ohis->CSliceZ(GRID=grid)

   ;; If time-varying zeta data from the history file are to be used,
   ;; then retrieve a list of times from this file. Plus other stuff

   if use_zeta then begin

      ;; Grid for zeta retrievals.

      zgrid = ohis->CsliceGrid('zeta', ALONG_RANGE=grid.along_range, $
                               DIRECTION=grid.direction, SLICE=grid.slice)

      ;; Time from history file, needed for temporal interpolation

      ztime = ohis->VarGet(ohis->TimeVarName())

      ;; Array to hold two frames of grid Z data

      ntemp = replicate(-1,2)
      gztemp = fltarr(n_elements(grid.arc), n_elements(grid.s), 2)

   end

   ;; Specify the float-file variable representing along-slice
   ;; position of the float.

   case grid.direction of
      0: avar = 'Xgrid'
      1: avar = 'Ygrid'
   endcase

   ;; Set default for DEPTH_RANGE

   if n_elements(depth_range) eq 0 then $
        depth_range = [-1.05,0.05]*max(grid.h)

   ;; Create base graph

   mgh_new, 'MGHgrGraph2D', ASPECT=0.7, $
            NAME='ROMS float x-z animation', RESULT=ograph, $
            _STRICT_EXTRA=graph_properties

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fsize

   ograph->NewFont, SIZE=fsize
   ograph->NewFont, SIZE=0.9*fsize

   ;; Draw axes

   ograph->NewAxis, 0, RANGE=mgh_minmax(grid.arc), /EXACT, $
        TITLE='Distance (km)', TICKFORMAT='mgh_tf_linear', $
        TICKFRMTDATA={scale:1.E-3, format:'(F10.1)'}
   ograph->NewAxis, 1, RANGE=depth_range, /EXACT, $
        TITLE='Depth (m)', TICKFORMAT='mgh_tf_negative'

   ;; Background to act as selection target.

   ograph->NewBackground

   ;; Plot bathymetry

   xb = grid.arc#replicate(1.,2)
   yb = [[mgh_reproduce(-1.05*max(grid.h), grid.h)],[-grid.h]]
   zb = make_array(DIMENSION=size(xb, /DIMENSIONS))
   ograph->NewAtom, 'IDLgrSurface', $
        DATAX=temporary(xb), DATAY=temporary(yb), DATAZ=temporary(zb), $
        STYLE=2, COLOR=[127,127,127], NAME='Bathymetry'

   ;; Plot walls

   walls = ohis->GetWalls()

   ara0 = grid.along_range[0]
   ara1 = grid.along_range[1]

   if walls[direction] and (ara0 eq 0) then begin
      arc = [grid.arc[0], mgh_avg(grid.arc[0:1])]
      ograph->NewAtom, 'IDLgrPolygon', arc[[0,1,1,0]], depth_range[[0,0,1,1]], $
           COLOR=[127,127,127], NAME='Wall 0'
   endif

   if walls[direction+2] and (ara1 eq dim_rho[direction]-1) then begin
      arc = [mgh_avg(grid.arc[ara1-1:ara1]), grid.arc[ara1]]
      ograph->NewAtom, 'IDLgrPolygon', arc[[0,1,1,0]], depth_range[[0,0,1,1]], $
           COLOR=[127,127,127], NAME='Wall 1'
   endif

   ;; Add various graphics objects that will (or might) be animated...

   ;; ...title

   otitle = ograph->NewTitle('')

   ;; ...invisible plot object with symbols indicating float
   ;; positions. Make the symbols reasonably simple to avoid slowing
   ;; the animation too much.

   if n_elements(symbol_properties) gt 0 then begin
      n_sym = n_elements(symbol_properties)
      osym = objarr(n_sym)
      for i=0,n_sym-1 do begin
         osym[i] = ograph->NewSymbol(0, N_VERTICES=6, COLOR=mgh_color('red'), $
                                     NORM_SIZE=0.01, $
                                     _STRICT_EXTRA=symbol_properties[i])
      endfor
   endif else begin
      osym = ograph->NewSymbol(0, N_VERTICES=6, COLOR=mgh_color('red'), NORM_SIZE=0.01)
   endelse

   ograph->NewAtom, 'IDLgrPlot', LINESTYLE=6, SYMBOL=osym, RESULT=oplt

   ;; ...mesh or ruled surface showing levels.

   ograph->NewAtom, 'IDLgrSurface', STYLE=3, COLOR=mgh_color('blue'), $
        DATAZ=fltarr(n_elements(grid.arc), n_elements(grid.s)), $
        DATAX=cmreplicate(grid.arc, n_elements(grid.s)), $
        DATAY=grid_z, RESULT=omesh

   ;; Create an animator window to display and manage the movie.

   ok = self->MGH_Datamator::Init(CHANGEABLE=0, GRAPHICS_TREE=ograph, $
                                  _STRICT_EXTRA=_extra)
   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Datamator'

   ;; Step through times in the netCDF file, generating new frames &
   ;; plotting data

   oframe = objarr(2+use_zeta)

   for r=0,n_elements(records)-1 do begin

      rec = records[r]

      ;; Check to see if the user has selected the "Finish Loading" menu item

      if self->Finished() then break

      ;; Update title

      t = oflt->VarGet('ocean_time', OFFSET=rec, COUNT=[1])
      t = t[0]
      oframe[0] = obj_new('MGH_Command', OBJECT=otitle, 'SetProperty', $
                          STRINGS=string(t/(24*3600), FORMAT='(F0.3)')+' days')

      ;; If necessary, recalculate surface height and grid z-values

      if use_zeta then begin

         if (t lt ztime[0]) || (t gt ztime[n_elements(ztime)-1]) then begin
            message, 'Time '+mgh_format_float(t)+' not found in ' + $
                     'history file'
         endif

         n0 = (value_locate(ztime, t))[0]
         n1 = (n0 + 1) < n0

         ;; Ensure grid Z data for record n0 is loaded into the first
         ;; location in the temporary array

         case 1 of
            n0 eq ntemp[0]:     ; No action necessary
            n0 eq ntemp[1]: begin
               gztemp[0,0,0] = gztemp[*,*,1]
               ntemp[0] = n0
            end
            else: begin
               zzzz = ohis->CtranData('zeta', GRID=zgrid, RECORD=n0) * exaggerate_zeta
               gztemp[0,0,0] = ohis->CsliceZ(GRID=grid, ZETA=temporary(zzzz))
               ntemp[0] = n0
            endelse
         endcase

         ;; Ditto for record n1

         case 1 of
            n1 eq ntemp[1]:     ; No action necessary
            else: begin
               zzzz = ohis->CtranData('zeta', GRID=zgrid, RECORD=n1) * exaggerate_zeta
               gztemp[0,0,1] = ohis->CsliceZ(GRID=grid, ZETA=temporary(zzzz))
               ntemp[1] = n1
            endelse
         endcase

         ;; Interpolate in time

         a0 = (n0 eq n1) ? 0. : (ztime[n1]-t)/(ztime[n1]-ztime[n0])

         grid_z = reform(a0*gztemp[*,*,0] + (1-a0)*gztemp[*,*,1])

      endif

      ;; Retrieve float horizontal grid position and project it onto
      ;; the slice,

      agg = oflt->VarGetFloat(avar, FLOATS=floats, RECORDS=rec) - grid.along_range[0]

      arc = mgh_interpolate(grid.arc, agg, MISSING=!values.f_nan)

      ;; Retrieve & project float vertical position

      zgg = oflt->VarGetFloat('Zgrid', FLOATS=floats, RECORDS=rec)

      dep = mgh_interpolate(grid_z, agg, zgg, MISSING=!values.f_nan)

      ;; Load float positions into plot object

      if isa(arc, /SCALAR) then arc = [arc]
      if isa(dep, /SCALAR) then dep = [dep]

      oframe[1] = obj_new('MGH_Command', OBJECT=oplt, 'SetProperty', $
                          DATAX=arc, DATAY=dep)

      ;; Update mesh if necessary

      if use_zeta then $
           oframe[2] = obj_new('MGH_Command', OBJECT=omesh, 'SetProperty', $
                               DATAY=grid_z)

      ;; Add frame to animator & display

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; MGH_ROMS_Movie_Cswarm::Cleanup
;
pro MGH_ROMS_Movie_Cswarm::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::Cleanup

end

; MGH_ROMS_Movie_Cswarm::GetProperty
;
pro MGH_ROMS_Movie_Cswarm::GetProperty, $
     ALL=all, FLOAT_FILE=float_file, GRID=grid, $
     HISTORY_FILE=history_file, _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::GetProperty, ALL=all, _STRICT_EXTRA=_extra

   float_file = self.float_file
   history_file = self.history_file

   if ptr_valid(self.grid) then grid = *self.grid

   if arg_present(all) then begin
      all = create_struct(all, 'float_file', float_file, $
                          'history_file', history_file)
   endif

end

; MGH_ROMS_Movie_Cswarm::SetProperty
;
pro MGH_ROMS_Movie_Cswarm::SetProperty, _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=_extra

end

; MGH_ROMS_Movie_Cswarm::About
;
;   Print information about the window and its contents
;
pro MGH_ROMS_Movie_Cswarm::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::About, lun

   self->GetProperty, $
        FLOAT_FILE=float_file, GRID=grid, HISTORY_FILE=history_file

   if obj_valid(float_file) then begin
      printf, lun, FORMAT='(%"%s: the float file sequence is %s. Its files are:")', $
              mgh_obj_string(self), mgh_obj_string(float_file)
      self.float_file->GetProperty, FILE_NAME=file_name
      print, file_name
   endif

   if obj_valid(history_file) then begin
      printf, lun, self, ': the history file sequence is '+ $
              mgh_obj_string(history_file)+'. Its files are:'
      self.history_file->GetProperty, FILE_NAME=file_name
      print, file_name
   endif

   printf, lun, FORMAT='(%"%s: the slice direction is %d")', $
           mgh_obj_string(self), grid.direction

   printf, lun, FORMAT='(%"%s: the slice index & range are %d %d %d")', $
           mgh_obj_string(self), grid.slice, grid.along_range

end

; MGH_ROMS_Movie_Cswarm::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro MGH_ROMS_Movie_Cswarm::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::BuildMenuBar

end


; MGH_ROMS_Movie_Cswarm::Event
;
function MGH_ROMS_Movie_Cswarm::Event, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   return, self->MGH_Datamator::Event(event)

end

; MGH_ROMS_Movie_Cswarm::ExportData
;
pro MGH_ROMS_Movie_Cswarm::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::ExportData, values, labels

   self->GetProperty, $
        ANIMATION=animation, FLOAT_FILE=float_file, $
        HISTORY_FILE=history_file, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[1]->GetProperty, KEYWORDS=keywords

   labels = [labels, 'Float Object', 'History Object', 'X,Z Data']
   values = [values, ptr_new(float_file), ptr_new(history_file), $
             ptr_new(complex(keywords.datax,keywords.datay))]

end

pro MGH_ROMS_Movie_Cswarm__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {MGH_ROMS_Movie_Cswarm, inherits MGH_Datamator, $
         float_file: obj_new(), history_file: obj_new(), $
         grid: ptr_new()}

end
