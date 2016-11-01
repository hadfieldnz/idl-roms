;+
; CLASS NAME:
;   MGH_ROMS_Movie_Swarm3d
;
; PURPOSE:
;   This procedure generates and displays an animated graph showing
;   ROMS float locations in 3D space
;
; CALLING SEQUENCE:
;   mgh_new, 'MGH_ROMS_Movie_Swarm3d', file
;
; INIT ARGUMENTS:
;  ffile
;    Synonym for FLOAT_FILE property.
;
;  hfile
;    Synonym for HISTORY_FILE property.
;
; PROPERTIES
;  ETA_RANGE (Init)
;    A 2-element integer vector specifying the grid points to be plotted
;    in the eta direction. The range is specified relative to the rho grid.
;    Negative values are taken to be offsets from the "northern" edge of
;    the grid. The default is [0,-1] which is equivalent to [0,dim(eta_rho)-1].
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
;   Mark Hadfield, 2011-08:
;     - Added support for map projections via the MAP_STRUCTURE
;       keyword.
;     - Removed FLOAT_DESTROY and HISTORY_DESTROY keywords.
;-
function MGH_ROMS_Movie_Swarm3d::Init, ffile, hfile, $
     DEPTH_RANGE=depth_range, $
     FLOAT_FILE=float_file, $
     FLOAT_RANGE=float_range, $
     FLOAT_STRIDE=float_stride, $
     FLOATS=floats, $
     GRAPH_PROPERTIES=graph_properties, $
     HISTORY_FILE=history_file, $
     MAP_STRUCTURE=map_structure, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     RECORDS=records, $
     SYMBOL_PROPERTIES=symbol_properties, $
     XI_RANGE=xi_range, $
     ETA_RANGE=eta_range, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
     ZAXIS_PROPERTIES=zaxis_properties, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process float-file argument

   if n_elements(float_file) eq 0 then $
        if n_elements(ffile) gt 0 then float_file = ffile

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

   if n_elements(history_file) eq 0 && n_elements(hfile) gt 0 then $
        history_file = hfile

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

   ;; Get number of floats and resolve float-selection arguments

   n_float = oflt->DimInfo('drifter', /DIMSIZE)

   mgh_resolve_indices, n_float, float_range, float_stride, floats

   ;; Get grid size and resolve grid-subset arguments

   dim_rho = ohis->DimRho()

   if n_elements(xi_range) eq 0 then xi_range = [0,-1]
   if xi_range[0] lt 0 then xi_range[0] += dim_rho[0]
   if xi_range[1] lt 0 then xi_range[1] += dim_rho[0]

   if n_elements(eta_range) eq 0 then eta_range = [0,-1]
   if eta_range[0] lt 0 then eta_range[0] += dim_rho[1]
   if eta_range[1] lt 0 then eta_range[1] += dim_rho[1]

   xra0 = xi_range[0]
   xra1 = xi_range[1]
   xran = xra1-xra0+1
   era0 = eta_range[0]
   era1 = eta_range[1]
   eran = era1-era0+1

   ;; Are longitude & latitude data available?

   lonlat = ohis->HasVar('lon_rho') && ohis->HasVar('lat_rho')

   ;; Are we using a map projection?

   if n_elements(map_structure) gt 0 then begin
      if ~ lonlat then $
           message, 'Cannot use map structure without lon, lat data'
      self.map_structure = ptr_new(map_structure)
   endif

   use_map_structure = ptr_valid(self.map_structure)

   ;; Get x & y positions @ RHO points. This is used for float
   ;; positions, walls and land mask.

   o = [xra0,era0]  &  c = [xran,eran]

   x_rho = ohis->VarGet(lonlat ? 'lon_rho' : 'x_rho', OFFSET=o, COUNT=c)
   y_rho = ohis->VarGet(lonlat ? 'lat_rho' : 'y_rho', OFFSET=o, COUNT=c)

   if use_map_structure then begin
      xy = map_proj_forward(x_rho, y_rho, MAP_STRUCTURE=*self.map_structure)
      x_rho = reform(xy[0,*], size(x_rho, /DIMENSIONS))
      y_rho = reform(xy[1,*], size(y_rho, /DIMENSIONS))
      mgh_undefine, xy
   endif

   ;; Establish variable and dimension name associated with float data

   case 1B of
      oflt->HasVar('Xgrid') && oflt->HasAtt('Xgrid','time'): $
           time_var = oflt->AttGet('Xgrid','time')
      oflt->HasVar('Ygrid') && oflt->HasAtt('Ygrid','time'): $
           time_var = oflt->AttGet('Ygrid','time')
      else: $
           time_var = 'ocean_time'
   endcase

   time_dim = (oflt->VarDimNames(time_var))[0]

   ;; Get bathymetry

   h = ohis->VarGet('h', OFFSET=[xra0,era0], COUNT=[xran,eran])

   if ohis->HasVar('mask_rho') then $
        h = h * ohis->VarGet('mask_rho', OFFSET=[xra0,era0], $
                             COUNT=[xran,eran])

   if n_elements(depth_range) eq 0 then depth_range = max(h)*[-0.05, 1.05]

   ;; Get number of records in float file and resolve record-selection
   ;; arguments

   n_time = oflt->DimInfo(time_dim, /DIMSIZE)

   mgh_resolve_indices, n_time, record_range, record_stride, records

   ;; Create base graph

   x_range = mgh_minmax(x_rho)
   y_range = mgh_minmax(y_rho)
   z_range = -reverse(depth_range)

   mgh_new, 'MGHgrGraph3D', PLOT_BOX=[-0.5,-0.5,-0.2,1.0,1.0,0.4], $
            NAME='ROMS float location animation', RESULT=ograph, $
            _STRICT_EXTRA=graph_properties

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   olmodel = ograph->Get(POSITION=2)

   ograph->NewAtom, 'IDLgrLight', MODEL=olmodel, LOCATION=[0.5,0.5,0.8], $
        TYPE=1, INTENSITY=0.7, NAME='Positional'
   ograph->NewAtom, 'IDLgrLight', MODEL=olmodel, TYPE=0, INTENSITY=0.5, $
        NAME='Ambient'

   ;; Draw axes

   if lonlat && (~ use_map_structure) then begin
      xap = {tickformat: 'mgh_tf_longitude', tickfrmtdata: {format:'(F10.1)'}}
      yap = {tickformat: 'mgh_tf_latitude', tickfrmtdata: {format:'(F10.1)'} }
   endif else begin
      xap = {title: 'X (km)', tickformat: 'mgh_tf_linear', $
             tickfrmtdata: {scale:1.E-3, format:'(F10.1)'}}
      yap = {title: 'Y (km)', tickformat: 'mgh_tf_linear', $
             tickfrmtdata: {scale:1.E-3, format:'(F10.1)'}}
   endelse

   ograph->NewAxis, 0, $
        RANGE=x_range, /EXACT, /EXTEND, $
        _STRICT_EXTRA=mgh_struct_merge(xap, xaxis_properties)
   ograph->NewAxis, 1, $
        RANGE=y_range, /EXACT, /EXTEND, $
        _STRICT_EXTRA=mgh_struct_merge(yap, yaxis_properties)

   ograph->NewAxis, 2, RANGE=z_range, /EXACT, TITLE='Depth (m)', $
        TICKFORMAT='mgh_tf_negative'

   ;; Draw land mask?

   ;; Draw bathymetry

   ograph->NewAtom, 'IDLgrSurface', STYLE=2, COLOR=mgh_color('light blue'), $
        DATAX=x_rho, DATAY=y_rho, DATAZ=-h, NAME='Bathymetry', RESULT=osurf
   self.surface = osurf

   ;; Create a symbol to represent each float. Make the symbol
   ;; reasonably simple to avoid slowing the animation too much.

   ograph->NewSymbol, 3, /FILL, N_VERTICES=6, $
        COLOR=mgh_color('red'), NORM_SIZE=0.008, RESULT=osym, $
        _STRICT_EXTRA=symbol_properties
   self.symbol = osym

   ;; Create an invisible polyline. Symbols will be displayed at the vertices

   ograph->NewAtom, 'IDLgrPolyline', LINESTYLE=6, SYMBOL=osym, RESULT=opoly

   ;; Create a title object

   otitle = ograph->NewTitle('')

   ;; Create an animator window to display and manage the movie.

   mouse_action = ['Rotate','Magnify 3D','Context']
   ok = self->MGH_Datamator::Init(CHANGEABLE=0, GRAPHICS_TREE=ograph, $
                                  MOUSE_ACTION=mouse_action, _STRICT_EXTRA=extra)
   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Datamator'

   ;; Step through times in the netCDF file, generating new frames &
   ;; plotting data

   oframe = objarr(2)

   for r=0,n_elements(records)-1 do begin

      rec = records[r]

      ;; Check to see if the user has selected the "Finish Loading"
      ;; menu item

      if self->Finished() then break

      ;; Title. TO DO: scale result according to "units" property

      t = oflt->VarGet(time_var, OFFSET=rec, COUNT=[1])
      tstring = string(t/(24*3600), FORMAT='(F0.3)')+' days'
      oframe[0] = obj_new('MGH_Command', OBJECT=otitle, 'SetProperty', STRINGS=tstring)

      ;; Get float position data

      xgrid = oflt->VarGetFloat('Xgrid', FLOATS=floats, RECORDS=rec, AUTOSCALE=0)
      ygrid = oflt->VarGetFloat('Ygrid', FLOATS=floats, RECORDS=rec, AUTOSCALE=0)

      ;; Float grid-relative positions should be in the range
      ;; [0.5,dim_rho-1.5]; values outside this range imply the float
      ;; has not been released yet or has become unbounded. (The xgrid
      ;; & ygrid values in this case are usually 1.E35, but for an MPI
      ;; run they may be zero.)

      l_bound = where(xgrid ge 0.5 and xgrid le dim_rho[0]-1.5 and $
                      ygrid ge 0.5 and ygrid le dim_rho[1]-1.5, n_bound, $
                      COMPLEMENT=l_unbound, NCOMPLEMENT=n_unbound)

      if n_unbound gt 0 then begin
         xgrid[l_unbound] = !values.f_nan
         ygrid[l_unbound] = !values.f_nan
      endif

      if n_bound gt 0 then begin
         xgrid[l_bound] -= xra0
         ygrid[l_bound] -= era0
      endif

      x = mgh_interpolate(x_rho, xgrid, ygrid, MISSING=!values.f_nan)
      y = mgh_interpolate(y_rho, xgrid, ygrid, MISSING=!values.f_nan)

      ;; Get vertical position data

      z = oflt->VarGetFloat('depth', FLOATS=floats, RECORDS=rec, AUTOSCALE=1)

      ;; Display

      if n_bound gt 0 then begin
         oframe[1] = $
              obj_new('MGH_Command', OBJECT=opoly, 'SetProperty', HIDE=0, $
                      DATA=transpose([[x[l_bound]],[y[l_bound]],[z[l_bound]]]))
      endif else begin
         oframe[1] = $
              obj_new('MGH_Command', OBJECT=opoly, 'SetProperty', HIDE=1)
      endelse

      ;; Add frame to animator & display

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; MGH_ROMS_Movie_Swarm3d::Cleanup
;
pro MGH_ROMS_Movie_Swarm3d::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Datamator::Cleanup

end

; MGH_ROMS_Movie_Swarm3d::GetProperty
;
pro MGH_ROMS_Movie_Swarm3d::GetProperty, $
     ALL=all, SURFACE_STYLE=surface_style, SYMBOL_SIZE=symbol_size, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::GetProperty, All=all, _STRICT_EXTRA=extra

   self.surface->GetProperty, STYLE=surface_style

   self.symbol->GetProperty, SIZE=symbol_size

   if arg_present(all) then $
        all = create_struct(all, 'surface_style', surface_style, $
                            'symbol_size', symbol_size)

end

; MGH_ROMS_Movie_Swarm3d::SetProperty
;
pro MGH_ROMS_Movie_Swarm3d::SetProperty, $
     SURFACE_STYLE=surface_style, SYMBOL_SIZE=symbol_size, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=extra

   if n_elements(surface_style) gt 0 then begin
      if surface_style gt 4 then message, 'Sorry I don''t do lego surfaces'
      self.surface->SetProperty, STYLE=surface_style
   endif

   self.symbol->SetProperty, SIZE=symbol_size

end

; MGH_ROMS_Movie_Swarm3d::About
;
;   Print information about the window and its contents
;
pro MGH_ROMS_Movie_Swarm3d::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::About, lun

   if obj_valid(self.float_file) then begin
      printf, lun, self, ': the float file sequence is ', $
              mgh_obj_string(self.float_file),'. Its files are:'
      self.float_file->GetProperty, FILE_NAME=file_name
      print, file_name
   endif

   if obj_valid(self.history_file) then begin
      printf, lun, self, ': the history file sequence is ', $
              mgh_obj_string(self.history_file),'. Its files are:'
      self.history_file->GetProperty, FILE_NAME=file_name
      print, file_name
   endif

end

; MGH_ROMS_Movie_Swarm3d::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro MGH_ROMS_Movie_Swarm3d::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::BuildMenuBar

   obar = mgh_widget_self(self.menu_bar)

   if obj_valid(obar) then begin

      obar->NewItem, PARENT='Tools', SEPARATOR=[1,0,0], $
        ['Set Surface Style...','Set Symbol Size...','Arrange Lights...']

   endif

end

; MGH_ROMS_Movie_Swarm3d::EventMenuBar
;
function MGH_ROMS_Movie_Swarm3d::EventMenuBar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      'TOOLS.SET SURFACE STYLE': begin
         mgh_new, 'MGH_GUI_SetList', CAPTION='Style', CLIENT=self, /FLOATING, $
                  GROUP_LEADER=self.base, /IMMEDIATE, $
                  ITEM_STRING=['Points','Mesh','Filled','Ruled XZ', $
                               'Ruled YZ'], $
                  PROPERTY_NAME='SURFACE_STYLE'
         return, 1
      end

      'TOOLS.SET SYMBOL SIZE': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Size', CLIENT=self $
                  , /FLOATING, GROUP_LEADER=self.base, /IMMEDIATE $
                  , N_ELEMENTS=3, PROPERTY_NAME='SYMBOL_SIZE'
         return, 1
      end

      'TOOLS.ARRANGE LIGHTS': begin
         self->GetProperty, GRAPHICS_TREE=graphics_tree
         olights = graphics_tree->Get(POSITION=2)
         mgh_new, 'MGH_GUI_LightEditor', CLIENT=self $
                  , LIGHT=olights->Get(/ALL) $
                  , /FLOATING, GROUP_LEADER=self.base, /IMMEDIATE
         return, 1
      end

      else: return, self->MGH_Datamator::EventMenuBar(event)

   endcase

end

pro MGH_ROMS_Movie_Swarm3d__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {MGH_ROMS_Movie_Swarm3d, inherits MGH_Datamator, $
         float_file: obj_new(), history_file: obj_new(), $
         surface: obj_new(), symbol: obj_new(), $
         map_structure: ptr_new()}

end
