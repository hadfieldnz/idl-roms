;+
; CLASS NAME:
;   MGH_ROMS_Movie_Hsurf
;
; PURPOSE:
;   This class generates and displays a an animated sequence of graphs showing
;   a slice through a ROMS 2D or 3D output field on an x-y surface in the form
;   of an IDLgrSurface.
;
; CATEGORY:
;   Ocean models.
;   Object graphics
;
; CALLING SEQUENCE:
;   mgh_new, 'mgh_roms_movie_hsurf', history, varname
;
; POSITIONAL PARAMETERS:
;   history
;     A reference to a ROMS history sequence object or a string arrary
;     specifying a list of ROMS history files osr a single string with
;     wildcards specifying a list of ROMS history files.
;
;   variable
;     The name of a 2-D or 3-D variable in the netCDF file.
;
; KEYWORD PARAMETERS:
;   DATA_MULTIPLIER (input, scalar numeric)
;     Number by which data values are multiplied before they are loaded into
;     the density surface. Default depends on the variable and is calculated by
;     MGH_ROMS_RESOLVE_DATA
;
;   DATA_RANGE (input, 2-element numeric)
;     Data range for the density surface. Default depends on the variable and
;     is calculated by MGH_ROMS_RESOLVE_DATA
;
;   DEPTH
;     Set this keyword to a scalar numeric value to specify the
;     depth of a z surface on which data are to be plotted. This
;     keyword should be specified only for variables having a depth
;     coordinate and it cannot be used together with LEVEL or SIGMA.
;
;   LAYER
;     Set this keyword to a scalar integer to specify the bed layer to
;     be plotted.  This keyword should be specified only for variables
;     having a bed-layer dimension.
;
;   LEVEL
;     Set this keyword to a scalar integer to specify the
;     s-coordinate level to be plotted.  This keyword should be
;     specified only for variables having a depth coordinate and it
;     cannot be used together with DEPTH or SIGMA.
;
;   MAP_STRUCTURE
;     If a map structure is supplied via this keyword, and if there
;     are longitude and latitude data in the file, then the longitudes
;     and latitudes are converted with MAP_PROJ_FORWARD before display.
;
;   SIGMA
;     Set this keyword to a scalar numeric value to specify the
;     sigma values of a copnstant-sigma surface on which data are to be
;     plotted. This keyword should be specified only for variables having
;     a depth coordinate and it cannot be used together with DEPTH or LEVEL.
;
;   XI_RANGE
;   ETA_RANGE
;     Use these keywords to display a subset of the domain. They are
;     interpreted relative to the rho grid.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, Apr 1999:
;     Written as ROMS_HSURF_MOVIE.
;   Mark Hadfield, Aug 1999:
;     Miscellaneous improvements including sequences of files & LEVEL
;     keyword.
;   Mark Hadfield, Jan 2000:
;     Converted from a procedure to a class with a front-end procedure
;     (MGH_ROMS_HSLICE_MOVIE) that can be used to create the object.
;   Mark Hadfield, Aug 2000:
;     MGH_ROMS_HSURF_MOVIE is now the class name; an object can be
;     conveniently created with the MGH_NEW procedure.
;   Mark Hadfield, Nov 2000:
;     Changed default ETA_RANGE & XI_RANGE so that only physically
;     realistic points are shown. Added USE_ZETA keyword. Increased
;     flexibility of the code handling the variable's dimensions, so
;     that variables without a time dimension can be plotted.
;   Mark Hadfield, Dec 2000:
;     Rewritten to take advantage of MGHromsHistory's new HsliceData &
;     HsliceGet methods.
;   Mark Hadfield, Jun 2001:
;     Rewritten as a subclass of MGH_Datamator.
;   Mark Hadfield, 2009-03:
;     Now MGH_ROMS_MOVIE_HSURF.
;   Mark Hadfield, 2012-02:
;     Brought over changes from MGH_ROMS_MOVIE_HSLICE.
;-
function MGH_ROMS_Movie_Hsurf::Init, $
     history, variable, $
     ANIMATION_PROPERTIES=animation_properties, $
     BYTE_RANGE=byte_range, $
     DATA_MULTIPLIER=data_multiplier, $
     DATA_RANGE=data_range, $
     DEPTH=depth, LAYER=layer, LEVEL=level, SIGMA=sigma, $
     GRAPH_PROPERTIES=graph_properties, $
     MAP_STRUCTURE=map_structure, $
     RECORD_AVERAGE=record_average, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     RECORDS=records, STYLE=style, $
     SHOW_TIME=show_time, $
     TITLE=title, USE_BATH=use_bath, USE_ZETA=use_zeta, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
     X_RANGE=x_range, Y_RANGE=y_range, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
     ZAXIS_PROPERTIES=zaxis_properties, $
     _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process history argument

   if n_elements(history) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'history'

   case size(history, /TNAME) of
      'STRING': begin
         ohis = obj_new('MGHromsHistory', history)
         self.history_file = ohis
      end
      'OBJREF': begin
         ohis = history
         self.history_file = history
      end
      else: $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'history'
   endcase

   ;; Check variable name

   if n_elements(variable) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'variable'

   if n_elements(variable) gt 1 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'variable'

   if size(variable, /TNAME) ne 'STRING' then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'variable'

   self.variable = variable

   ;; Other defaults

   mgh_roms_resolve_data, self.variable, $
        DATA_RANGE=data_range, DATA_MULTIPLIER=data_multiplier

   self.data_range = data_range

   if n_elements(map_structure) gt 0 then $
        self.map_structure = ptr_new(map_structure)

   use_map_structure = ptr_valid(self.map_structure)

   if n_elements(title) eq 0 then title = ''
   if n_elements(show_time) eq 0 then show_time = 1B

   ;; Default style is filled
   self.style = n_elements(style) gt 0 ? style : 2
   if self.style gt 4 then message, 'Sorry I don''t do lego surfaces'

   ;; Set ETA_RANGE and XI_RANGE. These specify subsets of the rho
   ;; grid.

   dim_rho = [ohis->DimInfo('xi_rho', /DIMSIZE), $
              ohis->DimInfo('eta_rho', /DIMSIZE)]

   if n_elements(xi_range) eq 0 then xi_range = [1,dim_rho[0]-2]
   if n_elements(eta_range) eq 0 then eta_range = [1,dim_rho[1]-2]

   if xi_range[0] lt 0 then xi_range[0] += dim_rho[0]
   if xi_range[1] lt 0 then xi_range[1] += dim_rho[0]

   if eta_range[0] lt 0 then eta_range[0] += dim_rho[1]
   if eta_range[1] lt 0 then eta_range[1] += dim_rho[1]

   ;; Get grid data required for horizontal slice retrievals

   var_xi_range = xi_range
   var_eta_range = eta_range

   var_dims = ohis->VarDims(variable)

   vdh0 = var_dims.horizontal[0]
   if vdh0 eq 'xi_u' || vdh0 eq 'xi_psi' then $
        var_xi_range += [-1,0]

   vdh1 = var_dims.horizontal[1]
   if vdh1 eq 'eta_v' || vdh1 eq 'eta_psi' then $
        var_eta_range += [-1,0]

   grid = ohis->HsliceGrid(variable, ETA_RANGE=var_eta_range, XI_RANGE=var_xi_range)

   self.lonlat = grid.lonlat

   ;; Get x & y positions for all RHO points (used for plotting walls
   ;; and land mask)

   x_rho = self.lonlat ? ohis->VarGet('lon_rho') : ohis->VarGet('x_rho')
   y_rho = self.lonlat ? ohis->VarGet('lat_rho') : ohis->VarGet('y_rho')

   ;; Convert all position data to map projection if appropriate

   if use_map_structure then begin
      if ~ grid.lonlat then $
           message, 'Cannot use map structure without lon, lat data'
      xy = map_proj_forward(grid.x, grid.y, MAP_STRUCTURE=*self.map_structure)
      grid.x = reform(xy[0,*], size(grid.x, /DIMENSIONS))
      grid.y = reform(xy[1,*], size(grid.y, /DIMENSIONS))
      mgh_undefine, xy
      xy = map_proj_forward(x_rho, y_rho, MAP_STRUCTURE=*self.map_structure)
      x_rho = reform(xy[0,*], size(x_rho, /DIMENSIONS))
      y_rho = reform(xy[1,*], size(y_rho, /DIMENSIONS))
      mgh_undefine, xy
   endif

   ;; The variable will be plotted on a density plane; specify its
   ;; style and vertex positions, also the value to use at masked points.

   case strjoin(grid.dims.horizontal, ' ') of
      'xi_rho eta_rho': begin
         if n_elements(style) eq 0 then style = 0
         if n_elements(mask_value) eq 0 then mask_value = !values.f_nan
      end
      'xi_u eta_u': begin
         if n_elements(style) eq 0 then style = 0
         if n_elements(mask_value) eq 0 then mask_value = 0
      end
      'xi_v eta_v': begin
         if n_elements(style) eq 0 then style = 0
         if n_elements(mask_value) eq 0 then mask_value = 0
      end
      'xi_psi eta_psi': begin
         if n_elements(style) eq 0 then style = 0
         if n_elements(mask_value) eq 0 then mask_value = !values.f_nan
      end
   endcase

   x_var = style eq 0 ? mgh_stagger(grid.x, DELTA=[1,1]) : grid.x
   y_var = style eq 0 ? mgh_stagger(grid.y, DELTA=[1,1]) : grid.y

   ;; Establish records to be plotted (if applicable). Establish name of
   ;; time variable--the logic here is a bit shaky, as conventions
   ;; for defining time dimensions and variables in ROMS files are
   ;; not entirely consistent.

   if n_elements(record_average) eq 0 then record_average = 1

   has_time = strlen(grid.dims.time) gt 0

   if has_time then begin
      n_time = ohis->DimInfo(grid.dims.time, /DIMSIZE)
      mgh_resolve_indices, n_time, record_range, record_stride, records
      n_records = n_elements(records)
      time_var = ohis->TimeVarName(grid.dims.time)
      if isa(time_var, /NULL) then message, 'Time variable not found'
   endif else begin
      n_records = 1
      ;; Check other record keywords for validity here?
      if record_average gt 1 then $
           message, 'Cannot average over records for time-independent variable'
   endelse

   n_frames = long(n_records)/long(record_average)

   ;; Get time units from the history file

   if has_time then begin
      if ohis->HasAtt(time_var, 'units') then begin
         time_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
      endif else begin
         time_units = {scale: 1}
      endelse
   endif

   ;; Default graph aspect ratio, can be overridden via GRAPH_PROPERTIES

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(x_var)
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(y_var)

   ;; Create base graph

   ograph = obj_new('MGHgrGraph3D', NAME='ROMS surface animation', $
                    _STRICT_EXTRA=graph_properties)

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Draw axes

   if grid.lonlat && (~ use_map_structure) then begin
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
   ograph->NewAxis, DIRECTION=2, RANGE=self.data_range, /EXACT, /EXTEND, $
        _STRICT_EXTRA=zaxis_properties

   ;; Add a model containing some lights This is necessary only for a
   ;; filled surface, but seems to have no disadvantages in other
   ;; cases, so do it always.

   olmodel = ograph->Get(POSITION=2)

   ograph->NewAtom, 'IDLgrLight', MODEL=olmodel, LOCATION=[0.5,0.5,0.8], $
        TYPE=1, INTENSITY=0.7, NAME='Positional'
   ograph->NewAtom, 'IDLgrLight', MODEL=olmodel, TYPE=0, INTENSITY=0.5, $
        NAME='Ambient'

   ;; Add various graphics objects that will (or might) be animated...

   ;; ...title

   otitle = ograph->NewTitle(title)

   ;; ...surface showing data values.

   ograph->NewAtom, 'IDLgrSurface', STYLE=self.style, $
        COLOR=[127,127,255], BOTTOM=[127,255,127], DATAX=grid.x, DATAY=grid.y, $
        DATAZ=mgh_reproduce(!values.f_nan,grid.x), /HIDDEN_LINES, RESULT=osurf
   self.surface = osurf

   ;; Create an animator window to display and manage the movie.

   ma = ['Rotate','Pick','Context']

   ok = self->MGH_Datamator::Init(CHANGEABLE=0, GRAPHICS_TREE=ograph, $
                                  MOUSE_ACTION=ma,  _STRICT_EXTRA=extra)
   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Datamator'

   ;; Step through the netCDF file, generating new frames &
   ;; plotting data

   oframe = objarr(2)

   ra = record_average
   dm = data_multiplier

   for f=0,n_frames-1 do begin

      if self->Finished() then break

      rec0 = ra*f

      time = 0  &  slice = 0

      for r=rec0,rec0+ra-1 do begin

         if has_time then begin
            t = ohis->VarGet(time_var, OFFSET=records[r], $
                             COUNT=[1], AUTOSCALE=0)*time_units.scale
            time += temporary(t)/double(ra)
            slice += ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                                      RECORD=records[r], $
                                      DEPTHS=depth, LAYERS=layer, LEVELS=level, $
                                      SIGMAS=sigma, USE_BATH=use_bath, USE_ZETA=use_zeta)
         endif else begin
            slice += ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                                      DEPTHS=depth, LAYERS=layer, LEVELS=level, $
                                      SIGMAS=sigma, USE_BATH=use_bath, USE_ZETA=use_zeta)
         endelse

      endfor

      oframe[0] = obj_new('MGH_Command', OBJECT=self.surface, 'SetProperty', $
                          DATAZ=dm*slice/float(ra))

      ;; Update title with time or date, if applicable

      if has_time && show_time gt 0 then begin
         case show_time of
            1: ttt = string(FORMAT='(%"%0.3f days")', time)
            2: ttt = mgh_dt_string(time+time_units.offset)
         endcase
         if strlen(title) gt 0 then $
              ttt = string(FORMAT='(%"%s: %s")', title, ttt)
         oframe[1] = obj_new('MGH_Command', OBJECT=otitle, $
                             'SetProperty', STRINGS=temporary(ttt))
      endif

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; MGH_ROMS_Movie_Hsurf::Cleanup
;
pro MGH_ROMS_Movie_Hsurf::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Datamator::Cleanup

end

; MGH_ROMS_Movie_Hsurf::GetProperty
;
pro MGH_ROMS_Movie_Hsurf::GetProperty, $
     ALL=all, DATA_RANGE=data_range, $
     HISTORY_FILE=history_file, LONLAT=lonlat, MAP_STRUCTURE=map_structure, $
     STYLE=style, VARIABLE=variable, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::GetProperty, ALL=all, _STRICT_EXTRA=extra

   data_range = self.data_range

   history_file = self.history_file

   lonlat = self.lonlat

   if arg_present(all) || arg_present(map_structure) then $
        map_structure = ptr_valid(self.map_structure) ? *self.map_structure : -1

   style = self.style

   if arg_present(all) then $
        all = create_struct(all, $
                            'data_range', data_range, $
                            'history_file', history_file, $
                            'lonlat', lonlat, 'map_structure', map_structure, $
                            'variable', variable)

end

; MGH_ROMS_Movie_Hsurf::SetProperty
;
pro MGH_ROMS_Movie_Hsurf::SetProperty, $
     DATA_RANGE=data_range, STYLE=style, _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(data_range) gt 0 then begin

      self.data_range = data_range

      self->GetProperty, GRAPHICS_TREE=ograph

      ;; Change data range of vertical axis/axes. Since these
      ;; are master-slave axes, this resets the scaling for all
      ;; slave atoms as well.

      zaxis = ograph->GetAxis(DIRECTION=2, /ALL, COUNT=n_zaxes)
      for i=0,n_zaxes-1 do $
           zaxis[i]->SetProperty, RANGE=self.data_range

   endif

   if n_elements(style) gt 0 then begin

      if style gt 4 then message, 'Sorry I don''t do lego surfaces'

      self.style = style

      self.surface->SetProperty, STYLE=self.style

   endif

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=_extra

end

; MGH_ROMS_Movie_Hsurf::About
;
pro MGH_ROMS_Movie_Hsurf::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::About, lun

   if obj_valid(self.history_file) then begin
      printf, lun, self, ': the history file sequence is ', self.history_file
      self.history_file->GetProperty, FILE_NAME=file_name
      printf, lun, self, ': the files are:', file_name
   endif

end

; MGH_ROMS_Movie_Hsurf::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro MGH_ROMS_Movie_Hsurf::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::BuildMenuBar

   obar = mgh_widget_self(self.menu_bar)

   obar->NewItem, PARENT='Tools', MENU=[1,0], $
        ['Data Range','View Data Values...']

   obar->NewItem, PARENT='Tools.Data Range', ['Set...','Fit this Frame']

end


; MGH_ROMS_Movie_Hsurf::EventMenuBar
;
function MGH_ROMS_Movie_Hsurf::EventMenuBar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      'TOOLS.DATA RANGE.SET': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Range', CLIENT=self, /FLOATING, $
                  GROUP_LEADER=self.base, IMMEDIATE=0, N_ELEMENTS=2, $
                  PROPERTY_NAME='DATA_RANGE'
         return, 0
      end

      'TOOLS.DATA RANGE.FIT THIS FRAME': begin
         self->GetProperty, POSITION=position
         oframe = self.animation->GetFrame(POSITION=position)
         oframe[0]->GetProperty, KEYWORDS=keywords
         data_range = mgh_minmax(keywords.dataz, /NAN)
         if data_range[0] eq data_range[1] then data_range += [-1,1]
         self->SetProperty, DATA_RANGE=data_range
         self->Update
         return, 0
      end

      'TOOLS.VIEW DATA VALUES': begin
         self->GetProperty, POSITION=position
         oframe = self.animation->GetFrame(POSITION=position)
         oframe[0]->GetProperty, KEYWORDS=keywords
         data_dims = size(keywords.dataz, /DIMENSIONS)
         ;; Call REFORM so that XVAREDIT cannot modify values
         xvaredit, reform(keywords.dataz), GROUP=self.base, $
                   X_SCROLL_SIZE=(data_dims[0] < 8), $
                   Y_SCROLL_SIZE=(data_dims[1] < 30)
         return, 0
      end

      else: return, self->MGH_Datamator::EventMenuBar(event)

   endcase

end

; MGH_ROMS_Movie_Hsurf::ExportData
;
pro MGH_ROMS_Movie_Hsurf::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::ExportData, values, labels

   self->GetProperty, $
        ANIMATION=animation, HISTORY_FILE=history_file, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[0]->GetProperty, KEYWORDS=keywords

   labels = [labels,'History Object','Surface Data']
   values = [values,ptr_new(history_file),ptr_new(keywords.dataz)]

end


pro MGH_ROMS_Movie_Hsurf__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {MGH_ROMS_Movie_Hsurf, inherits MGH_Datamator, $
         history_file: obj_new(), $
         variable: '', lonlat: 0B, map_structure: ptr_new(), $
         data_range: fltarr(2), style: 0B, surface: obj_new()}

end
