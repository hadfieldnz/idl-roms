;+
; CLASS NAME:
;   MGH_Roms_Plot_Hstats
;
; PURPOSE:
;   This class generates and displays a graph showing statistical
;   quantities derived from a series of Hslices through a ROMS 2D or
;   3D output field.
;
; CATEGORY:
;   Ocean models.
;   Object graphics
;
; CALLING SEQUENCE:
;   mgh_new, 'mgh_roms_plot_hstats', history, variable
;
; POSITIONAL ARGUMENTS:
;   history
;     A reference to a ROMS history sequence object or a string arrary
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
;   variable
;     The name of a 2-D or 3-D variable in the netCDF file.
;
; KEYWORD ARGUMENTS:
;   DATA_MULTIPLIER (input, numeric scalar)
;     Number by which data values are multiplied before statistical processing.
;     The default depends on the variable and is calculated by MGH_ROMS_RESOLVE_DATA
;
;   DATA_RANGE (input, numeric 2-element vector)
;     Data range for the density surface. Default depends on the variable and
;     is calculated by MGH_ROMS_RESOLVE_DATA
;
;   DATA_TRANSFORMATION (input, numeric string)
;     Function applied to data values before statistical processing.
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
;   RECORD_RANGE (input, numeric 2-element vector)
;     Range of record numbers (0-based indices in the time dimension) over
;     which statistics are to be calculated. The record range may also be
;     specified indirectly via the TIME_RANGE keyword.
;
;   SIGMA
;     Set this keyword to a scalar numeric value to specify the
;     sigma values of a constant-sigma surface on which data are to be
;     plotted. This keyword should be specified only for variables having
;     a depth coordinate and it cannot be used together with DEPTH or LEVEL.
;
;   TIME_RANGE (input, 2-element numeric)
;     Time range (days) over which statistics are to be calculated. The
;     record range may also be specified directly via the RECORD_RANgE keyword.
;
;   XI_RANGE
;   ETA_RANGE
;     Use these keywords to display a subset of the domain. They are
;     interpreted relative to the rho grid.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2004-??:
;     Written.
;   Mark Hadfield, 2005-10:
;     Added support for RECORD_RANGE keyword.
;   Mark Hadfield, 2009-03:
;     Added support for sigma slices via the SIGMA keyword.
;   Mark Hadfield, 2009-10:
;     Removed calls to widget_event(/NOWAIT).
;   Mark Hadfield, 2010-07:
;     Modified code for the mean parameter to reduce memory
;     requirements.
;   Mark Hadfield, 2011-03:
;     Default colour table is now Matlab Jet.
;   Mark Hadfield, 2011-08:
;     - No longer destroys the history-file object.
;     - Modified code for the min & max parameters along the same lines as
;       the mean parameter to reduce memory requirements.
;     - Now handles STRUCT variables reasonably gracefully.
;   Mark Hadfield, 2012-??:
;     - Added TMPFILE property, allowing the caller to specify a file for
;       temporary storage of statistics.
;   Mark Hadfield, 2012-10:
;     - Added DATA_TRANSFORMATION and LOGARITHMIC properties.
;   Mark Hadfield, 2014-08:
;     - Added POST_TRANSFORMATION property.
;     - Added SHOW_TITLE keyword.
;     - Land mask default colour is now a lighter grey.
;   Mark Hadfield, 2015-12:
;     - The code to calculate percentiles now uses Craig Markwardt's function
;       CMAPPLY.
;     - Note that this routine is being largely superseded by the routine
;       MGH_ROMS_HISTORY (I may change that name) which pre-calculates
;       statistics and stores them in a netCDF file.
;   Mark Hadfield, 2016-02:
;     - Temporary file name now calculated in the Init method so does not need to
;       be (and cannot be) passed from the caller.
;   Mark Hadfield, 2016-03:
;     Fixed bug: forgot to include the parameter in the temporary file name.
;   Mark Hadfield, 2016-03:
;     Major change: statistics now calculated by a new function mgh_roms_hslice_statistics.
;-
function Mgh_Roms_Plot_Hstats::Init, $
     history, variable, $
     BYTE_RANGE=byte_range, $
     COLORBAR_PROPERTIES=colorbar_properties, $
     CONTOUR_PROPERTIES=contour_properties, $
     DATA_RANGE=data_range, LOGARITHMIC=logarithmic, $
     DATA_MULTIPLIER=data_multiplier, $
     POST_MULTIPLIER=post_multiplier, $
     DATA_TRANSFORMATION=data_transformation, $
     POST_TRANSFORMATION=post_transformation, $
     DEPTH=depth, LEVEL=level, SIGMA=sigma, $
     MASK_VALUE=mask_value, $
     MAP_STRUCTURE=map_structure, $
     PARAMETER=parameter, STYLE=style, $
     RECORD_RANGE=record_range, TIME_RANGE=time_range, $
     SHOW_CONTOUR=show_contour, SHOW_PLANE=show_plane, $
     SHOW_TITLE=show_title, SHOW_TIME=show_time, TITLE=title, DT_FORMAT=dt_format, $
     RECALC=recalc, USE_ZETA=use_zeta, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
     X_RANGE=x_range, Y_RANGE=y_range, $
     GRAPH_PROPERTIES=graph_properties, $
     LAND_PROPERTIES=land_properties, $
     PALETTE_PROPERTIES=palette_properties, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process variable name argument

   if n_elements(variable) eq 0 then variable = 'zeta'

   if ~ (isa(variable, 'STRING') || isa(variable, 'STRUCT'))  then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'variable'

   if n_elements(variable) ne 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'variable'

   self.variable = variable

   ;; Specify statistical parameter to be calculated

   if n_elements(parameter) eq 0 then parameter = 'mean'

   self.parameter = parameter

   ;; Process history argument.

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

   ;; Other defaults

   mgh_roms_resolve_data, variable, $
      DATA_RANGE=data_range, DATA_MULTIPLIER=data_multiplier

   if n_elements(show_contour) eq 0 then show_contour = !false

   if n_elements(show_plane) eq 0 then show_plane = !true

   if n_elements(map_structure) gt 0 then $
      self.map_structure = ptr_new(map_structure)

   use_map_structure = ptr_valid(self.map_structure)

   if n_elements(show_title) eq 0 then show_title = !true

   if keyword_set(show_title) then begin
      if n_elements(title) eq 0 then title = ''
      if n_elements(show_time) eq 0 then show_time = 2
      if n_elements(dt_format) eq 0 then $
         dt_format = '(C(CYI4.4,"-",CMOI2.2,"-",CDI2.2))'
   endif

   ;; Set ETA_RANGE and XI_RANGE relative to the rho grid.
   ;; The default is to show all  points.

   dim_rho = [ohis->DimInfo('xi_rho', /DIMSIZE),ohis->DimInfo('eta_rho', /DIMSIZE)]

   if n_elements(xi_range) eq 0 then xi_range = [0,dim_rho[0]-1]
   if n_elements(eta_range) eq 0 then eta_range = [0,dim_rho[1]-1]

   if xi_range[0] lt 0 then xi_range[0] += dim_rho[0]
   if xi_range[1] lt 0 then xi_range[1] += dim_rho[0]

   if eta_range[0] lt 0 then eta_range[0] += dim_rho[1]
   if eta_range[1] lt 0 then eta_range[1] += dim_rho[1]

   xra0 = xi_range[0]
   xra1 = xi_range[1]
   xran = xra1-xra0+1
   era0 = eta_range[0]
   era1 = eta_range[1]
   eran = era1-era0+1

   ;; Get the data

   data = mgh_roms_hslice_statistics(ohis, variable, RECALC=recalc, $
      DATA_MULTIPLIER=data_multiplier, DATA_TRANSFORMATION=data_transformation, $
      DEPTH=depth, LEVEL=level, SIGMA=sigma, PARAMETER=parameter, $
      RECORD_RANGE=record_range, TIME_RANGE=time_range, $
      XI_RANGE=xi_range, ETA_RANGE=eta_range)

   ;; Get x & y positions for all RHO points (used for plotting walls
   ;; and land mask)

   self.lonlat = data.grid.lonlat

   x_rho = self.lonlat ? ohis->VarGet('lon_rho') : ohis->VarGet('x_rho')
   y_rho = self.lonlat ? ohis->VarGet('lat_rho') : ohis->VarGet('y_rho')

   ;; Convert all position data to map projection if appropriate

   if use_map_structure then begin
      if ~ self.lonlat then $
         message, 'Cannot use map structure without lon, lat data'
      xy = map_proj_forward(data.grid.x, data.grid.y, MAP_STRUCTURE=*self.map_structure)
      data.grid.x = reform(xy[0,*], size(data.grid.x, /DIMENSIONS))
      data.grid.y = reform(xy[1,*], size(data.grid.y, /DIMENSIONS))
      mgh_undefine, xy
      xy = map_proj_forward(x_rho, y_rho, MAP_STRUCTURE=*self.map_structure)
      x_rho = reform(xy[0,*], size(x_rho, /DIMENSIONS))
      y_rho = reform(xy[1,*], size(y_rho, /DIMENSIONS))
      mgh_undefine, xy
   endif

   ;; The variable will be plotted on a density plane; specify its
   ;; style and vertex positions, also the value to use at masked points.

   case strjoin(data.grid.dims.horizontal, ' ') of
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

   x_var = style eq 0 ? mgh_stagger(data.grid.x, DELTA=[1,1]) : data.grid.x
   y_var = style eq 0 ? mgh_stagger(data.grid.y, DELTA=[1,1]) : data.grid.y

   ;; The POST_MULTIPLIER and POST_TRANSFORMATION arguments are applied *after*
   ;; statistical processing

   if n_elements(post_multiplier) gt 0 then $
      data.values *= post_multiplier
   if n_elements(post_transformation) gt 0 then $
      data.values = call_function(post_transformation, data.values)

   ;; Apply the mask value to the data

   if mgh_struct_has_tag(data.grid, 'mask') then $
      data.values[where(data.grid.mask eq 0, /NULL)] = mask_value

   ;; Default graph aspect ratio, can be overridden via GRAPH_PROPERTIES

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(x_var)
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(y_var)

   aspect = mgh_aspect(x_range, y_range, LONLAT=self.lonlat && (~ use_map_structure))
   aspect = (aspect > 0.4) < 1.5

   ;; Create base graph

   ograph = obj_new('MGHgrGraph2D', ASPECT=aspect, XMARGIN=[0.30,0.45], $
      NAME='ROMS horizontal slice statistics', $
      _STRICT_EXTRA=graph_properties)

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Add mask around plot area

   ograph->NewMask

   ;; Add title

   if keyword_set(show_title) then begin
      case show_time of
         0: begin
            ttt = string(FORMAT='(%"%s %s")', self.variable, self.parameter)
         end
         1: begin
            ttt = string(FORMAT='(%"%s %s: %0.3f to %0.3f days")', self.variable, self.parameter, time_range)
         end
         2: begin
            ttt = string(FORMAT='(%"%s %s: %s to %s")', self.variable, self.parameter, mgh_dt_string(time_range+data.time_offset, FORMAT=dt_format))
         end
      endcase
      if strlen(title) gt 0 then $
         ttt = string(FORMAT='(%"%s %s")', title, ttt)
      ograph->NewTitle, ttt
   endif
   ;; Draw axes

   if self.lonlat && (~ use_map_structure) then begin
      xap = {tickformat: 'mgh_tf_longitude', tickfrmtdata: {format:'(F10.1)'}}
      yap = {tickformat: 'mgh_tf_latitude', tickfrmtdata: {format:'(F10.1)'} }
   endif else begin
      xap = {title: 'X (km)', tickformat: 'mgh_tf_linear', tickfrmtdata: {scale:1.E-3, format:'(F0.1)'}}
      yap = {title: 'Y (km)', tickformat: 'mgh_tf_linear', tickfrmtdata: {scale:1.E-3, format:'(F0.1)'}}
   endelse

   ograph->NewAxis, 0, $
      RANGE=x_range, /EXACT, EXTEND=0, $
      _STRICT_EXTRA=mgh_struct_merge(xap, xaxis_properties)
   ograph->NewAxis, 1, $
      RANGE=y_range, /EXACT, EXTEND=0, $
      _STRICT_EXTRA=mgh_struct_merge(yap, yaxis_properties)

   ;; Add walls. For each wall we extract the 2 rows/columns on each
   ;; side of the physical boundary from the x_rho & y_rho arrays into
   ;; variables xr & yr, call mgh_stagger to get xw & yw along the
   ;; wall, then trim the interior points off.

   walls = ohis->GetWalls()

   for i_wall=0,3 do begin

      if walls[i_wall] gt 0 then begin

         case i_wall of
            0: begin
               ;; Western wall
               xw = mgh_stagger(x_rho[0:1,*], DELTA=[1,1])
               yw = mgh_stagger(y_rho[0:1,*], DELTA=[1,1])
               xw = xw[0:1,*]
               yw = yw[0:1,*]
            end
            1: begin
               ;; Southern wall
               xw = mgh_stagger(x_rho[*,0:1], DELTA=[1,1])
               yw = mgh_stagger(y_rho[*,0:1], DELTA=[1,1])
               xw = xw[*,0:1]
               yw = yw[*,0:1]
            end
            2: begin
               ;; Eastern wall
               xw = mgh_stagger(x_rho[dim_rho[0]-2:dim_rho[0]-1,*], DELTA=[1,1])
               yw = mgh_stagger(y_rho[dim_rho[0]-2:dim_rho[0]-1,*], DELTA=[1,1])
               xw = xw[1:2,*]
               yw = yw[1:2,*]
            end
            3: begin
               ;; Northern wall
               xw = mgh_stagger(x_rho[*,dim_rho[1]-2:dim_rho[1]-1], DELTA=[1,1])
               yw = mgh_stagger(y_rho[*,dim_rho[1]-2:dim_rho[1]-1], DELTA=[1,1])
               xw = xw[*,1:2]
               yw = yw[*,1:2]
            end
         endcase

         zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
         ograph->NewAtom, 'IDLgrSurface', $
            DATAX=temporary(xw), DATAY=temporary(yw), DATAZ=temporary(zw), $
            STYLE=2, COLOR=[127,127,127]

      endif

   endfor

   ;; Draw land mask.

   if ohis->HasVar('mask_rho') then begin

      mask_missing = round(ohis->VarGet('mask_rho'))

      ograph->NewAtom, 'MGHgrColorPlane', RESULT=oland, $
         STYLE=0, DEFAULT_COLOR=mgh_color('grey'), $
         ZVALUE=-2*deltaz, MISSING_POINTS=mask_missing, $
         DATAX=mgh_stagger(x_rho, DELTA=[1,1]), $
         DATAY=mgh_stagger(y_rho, DELTA=[1,1]), $
         NAME='Land mask', $
         _STRICT_EXTRA=land_properties

      ;; Just playing really
      if oland->QueryProperty('STYLE') then $
         oland->SetPropertyAttribute, 'STYLE', SENSITIVE=0

   endif

   ;; Draw data

   ograph->NewPalette, 'Matlab Jet', RESULT=palette, $
      _STRICT_EXTRA=palette_properties

   ograph->NewColorBar, RESULT=obar, FONT=ograph->GetFont(POS=1), $
      DATA_RANGE=data_range, LOGARITHMIC=logarithmic, PALETTE=palette, $
      SHOW_CONTOUR=show_contour, CONTOUR_PROPERTIES=contour_properties, $
      NAME='Colour bar', _STRICT_EXTRA=colorbar_properties
   self.bar = obar

   ograph->NewAtom, 'MGHgrDensityPlane', RESULT=oplane, $
      HIDE=(~ show_plane), $
      DATAX=x_var, DATAY=y_var, DATA_VALUES=data.values, $
      STYLE=style, ZVALUE=-5*deltaz , NAME='Data plane', $
      COLORSCALE=obar, /STORE_DATA
   self.plane = oplane

   ograph->NewAtom, 'IDLgrContour', RESULT=ocont, $
      HIDE=(~ show_contour), $
      GEOMZ=deltaz, /PLANAR, DATA=data.values, GEOMX=data.grid.x, GEOMY=data.grid.y, $
      _STRICT_EXTRA=contour_properties
   self.contour = ocont

   ;; Load graph into window

   ma = ['Magnify','Translate','Context']
   ok = self->MGH_Window::Init(GRAPHICS_TREE=ograph, MOUSE_ACTION=ma, _STRICT_EXTRA=extra)
   if ~ ok then $
      message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   ;; Finalise plot & return

   self->Finalize, 'Mgh_Roms_Plot_Hstats'

   return, 1

end

; Mgh_Roms_Plot_Hstats::Cleanup
;
pro Mgh_Roms_Plot_Hstats::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Window::Cleanup

end

; Mgh_Roms_Plot_Hstats::GetProperty
;
pro Mgh_Roms_Plot_Hstats::GetProperty, $
     ALL=all, BAR=bar, BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
     HISTORY_FILE=history_file, LONLAT=lonlat, MAP_STRUCTURE=map_structure, $
     PALETTE=palette, PARAMETER=parameter, PLANE=plane, STYLE=style, VARIABLE=variable, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::GetProperty, ALL=all, _STRICT_EXTRA=extra

   bar = self.bar

   plane = self.plane

   history_file = self.history_file

   lonlat = self.lonlat

   if arg_present(all) || arg_present(map_structure) then $
        map_structure = ptr_valid(self.map_structure) ? *self.map_structure : -1

   parameter = self.parameter

   variable = self.variable

   self.plane->GetProperty, STYLE=style

   self.bar->GetProperty, $
        BYTE_RANGE=byte_range, DATA_RANGE=data_range, PALETTE=palette

   if arg_present(all) then $
        all = create_struct(all, 'bar', bar, 'byte_range', byte_range, $
                            'data_range', data_range, $
                            'history_file', history_file, $
                            'lonlat', lonlat, 'map_structure', map_structure, $
                            'palette', palette, 'parameter', parameter, 'style', style, $
                            'variable', variable)

end

; Mgh_Roms_Plot_Hstats::SetProperty
;
pro Mgh_Roms_Plot_Hstats::SetProperty, $
     BYTE_RANGE=byte_range, DATA_RANGE=data_range, STYLE=style, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(byte_range) gt 0 || n_elements(data_range) gt 0 then begin

      if n_elements(byte_range) eq 0 then $
           self->GetProperty, BYTE_RANGE=byte_range
      if n_elements(data_range) eq 0 then $
           self->GetProperty, DATA_RANGE=data_range

      if obj_valid(self.bar) then $
           self.bar->SetProperty, BYTE_RANGE=byte_range, DATA_RANGE=data_range

      self.plane->SetProperty, BYTE_RANGE=byte_range, DATA_RANGE=data_range

   endif

   self.plane->SetProperty, STYLE=style

   self->MGH_Window::SetProperty, _STRICT_EXTRA=extra

end

; Mgh_Roms_Plot_Hstats::About
;
pro Mgh_Roms_Plot_Hstats::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::About, lun

   self->GetProperty, BAR=bar, HISTORY_FILE=history_file, $
        PALETTE=palette, VARIABLE=variable

   if obj_valid(history_file) then begin
      printf, lun, FORMAT='(%"%s: my history file sequence is %s")', $
           mgh_obj_string(self), mgh_obj_string(history_file)
      history_file->GetProperty, FILE_NAME=file_name
      printf, lun, self, ': its files are:', file_name
   endif

   printf, lun, FORMAT='(%"%s: my variable name is %s")', $
        mgh_obj_string(self), variable

   if obj_valid(self.bar) then begin
      printf, lun, self, ': the colour bar is ' , $
              mgh_obj_string(self.bar, /SHOW_NAME)
   endif

   if obj_valid(self.plane) then begin
      printf, lun, self, ': the density plane object is ', $
              mgh_obj_string(self.plane, /SHOW_NAME)
   endif

   if obj_valid(palette) then begin
      printf, lun, self, ': the palette is ', $
              mgh_obj_string(palette, /SHOW_NAME)
   endif

   if obj_valid(self.contour) then begin
      printf, lun, self, ': the contour is ' , $
              mgh_obj_string(self.contour, /SHOW_NAME)
   endif

end

; Mgh_Roms_Plot_Hstats::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Plot_Hstats::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::BuildMenuBar

   obar = mgh_widget_self(self.menu_bar)

   obar->NewItem, PARENT='File.Export', 'NetCDF...'

   obar->NewItem, PARENT='Tools', SEPARATOR=[1,0,0,0,1], MENU=[1,0,0,0,0], $
        ['Data Range','Edit Palette...','Set Style...', $
         'View Colour Scale...','View Data Values...']

   obar->NewItem, PARENT='Tools.Data Range', ['Set...','Fit Data']

end

; Mgh_Roms_Plot_Hstats::EventMenuBar
;
function Mgh_Roms_Plot_Hstats::EventMenuBar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      'FILE.EXPORT.NETCDF': begin
         self.graphics_tree->GetProperty, NAME=name
         ext = '.nc'
         default_file = strlen(name) gt 0 ? mgh_str_vanilla(name)+ext : ''
         filename = dialog_pickfile(/WRITE, FILE=default_file, FILTER='*'+ext)
         if strlen(filename) gt 0 then begin
            widget_control, HOURGLASS=1
            if !mgh_prefs.sticky then begin
               dir = file_dirname(filename)
               if strlen(dir) gt 0 then begin
                  cd, CURRENT=old_dir
                  if dir ne old_dir then begin
                     message, /INFORM, string(dir, FORMAT='(%"Changing to directory %s)")')
                     cd, dir
                  endif
               endif
            endif
            self->ExportToNcFile, filename
         endif
         return, 0
      end

      'TOOLS.DATA RANGE.SET': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Range', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, IMMEDIATE=0, $
                  N_ELEMENTS=2, PROPERTY_NAME='DATA_RANGE'
         return, 0
      end

      'TOOLS.DATA RANGE.FIT DATA': begin
         self.plane->GetProperty, DATA_VALUES=data_values
         data_range = mgh_minmax(data_values, /NAN)
         if min(finite(data_range)) gt 0 then begin
            if data_range[0] eq data_range[1] then data_range += [-1,1]
            self->SetProperty, DATA_RANGE=data_range
            self->Update
         endif
         return, 0
      end

      'TOOLS.EDIT PALETTE': begin
         self->GetProperty, PALETTE=palette
         mgh_new, 'MGH_GUI_Palette_Editor', palette, CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, /IMMEDIATE
         return, 0
      end

      'TOOLS.SET STYLE': begin
         mgh_new, 'MGH_GUI_SetList', CAPTION='Style', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, IMMEDIATE=0, $
                  ITEM_STRING=['Block','Interpolated'], $
                  PROPERTY_NAME='STYLE'
         return, 0
      end

      'TOOLS.VIEW COLOUR SCALE': begin
         mgh_new, 'MGH_GUI_ColorScale', CLIENT=self, /FLOATING, $
                  GROUP_LEADER=self.base
         return, 0
      end

      'TOOLS.VIEW DATA VALUES': begin
         self.plane->GetProperty, DATA_VALUES=data_values
         data_dims = size(data_values, /DIMENSIONS)
         ;; Call REFORM so that XVAREDIT cannot modify values
         xvaredit, reform(data_values), GROUP=self.base, $
                   X_SCROLL_SIZE=(data_dims[0] < 8), $
                   Y_SCROLL_SIZE=(data_dims[1] < 30)
         return, 0
      end

      else: return, self->MGH_Window::EventMenuBar(event)

   endcase

end

; Mgh_Roms_Plot_Hstats::ExportData
;
pro Mgh_Roms_Plot_Hstats::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::ExportData, values, labels

   self.plane->GetProperty, $
         DATA_VALUES=data_values, DATAX=datax, DATAY=datay

   if n_elements(data_values) gt 0 then begin
      labels = [labels,'Data Values']
      values = [values,ptr_new(data_values, /NO_COPY)]
   endif

   if n_elements(datax)*n_elements(datay) gt 0 then begin
      if self.lonlat then begin
         if ptr_valid(self.map_structure) then begin
           dim = size(datax, /DIMENSIONS)
           lonlat = map_proj_inverse(datax, datay, $
                                     MAP_STRUCTURE=*self.map_structure)
           lon = reform(lonlat[0,*], dim)
           lat = reform(lonlat[1,*], dim)
           mgh_undefine, lonlat
           labels = [labels,'Vertex Lon','Vertex Lat', $
                     'Vertex X','Vertex Y','Map structure']
           values = [values,ptr_new(lon, /NO_COPY),ptr_new(lat, /NO_COPY), $
                     ptr_new(datax, /NO_COPY),ptr_new(datay, /NO_COPY), $
                     ptr_new(*self.map_structure)]
         endif else begin
           labels = [labels,'Vertex Lon','Vertex Lat']
           values = [values,ptr_new(datax, /NO_COPY),ptr_new(datay, /NO_COPY)]
         endelse
      endif else begin
         labels = [labels,'Vertex X','Vertex Y']
         values = [values,ptr_new(datax, /NO_COPY),ptr_new(datay, /NO_COPY)]
      endelse
   endif

end

; Mgh_Roms_Plot_Hstats::ExportToNcFile
;
pro Mgh_Roms_Plot_Hstats::ExportToNcFile, file

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->GetProperty, $
        LONLAT=lonlat, MAP_STRUCTURE=map_structure, $
        PARAMETER=parameter, VARIABLE=variable

   self.plane->GetProperty, $
        DATA_VALUES=data_values, DATAX=datax, DATAY=datay, STYLE=style

   if style eq 0 then begin
      datax = mgh_stagger(datax, DELTA=[-1,-1])
      datay = mgh_stagger(datay, DELTA=[-1,-1])
   endif

   dim = size(data_values, /DIMENSIONS)

   if lonlat && size(map_structure, /TYPE) eq 8 then begin
      ll = map_proj_inverse(datax, datay, MAP_STRUCTURE=map_structure)
      datax = (reform(ll[0,*], dim)+360) mod 360
      datay = reform(ll[1,*], dim)
   endif

   l_miss = where(~ finite(data_values), n_miss)
   if n_miss gt 0 then data_values[l_miss] = mgh_ncdf_fill()

   fmt = '(%"Writing data (dim %d x %d) to output file %s")'
   message, /INFORM, string(dim, file, FORMAT=fmt)

   onc = obj_new('MGHncFile', file, /CREATE, /CLOBBER)

   onc->AttAdd, /GLOBAL, 'title', 'ROMS hslice statistics data'

   onc->AttAdd, /GLOBAL, 'history', $
                'Generated by routine MGH_ROMS_Plot_Hstats::ExportToNcFile at '+ $
                mgh_dt_string(mgh_dt_now())

   onc->AttAdd, /GLOBAL, 'variable', variable
   onc->AttAdd, /GLOBAL, 'parameter', parameter

   onc->DimAdd, 'xi', dim[0]
   onc->DimAdd, 'eta', dim[1]

   x_name = lonlat ? 'lon' : 'x'
   y_name = lonlat ? 'lat' : 'y'
   v_name = mgh_str_vanilla(variable+'_'+parameter)

   onc->VarAdd, x_name, ['xi','eta'], /DOUBLE
   onc->VarAdd, y_name, ['xi','eta'], /DOUBLE
   onc->VarAdd, v_name, ['xi','eta'], /FLOAT

   onc->VarPut, x_name, datax
   onc->VarPut, y_name, datay
   onc->VarPut, v_name, data_values

   obj_destroy, onc


end

; Mgh_Roms_Plot_Hstats::PickReport
;
pro Mgh_Roms_Plot_Hstats::PickReport, pos, LUN=lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(lun) eq 0 then lun = -1

   if n_elements(pos) ne 2 then $
        message, 'Parameter POS must be a 2-element vector'

   self->GetProperty, GRAPHICS_TREE=ograph

   if ~ obj_valid(ograph) then begin
      printf, lun, FORMAT='(%"%s: no graphics tree")', mgh_obj_string(self)
      return
   endif

   printf, lun, FORMAT='(%"%s: graphics tree %s")', $
           mgh_obj_string(self), mgh_obj_string(ograph, /SHOW_NAME)

   atoms = self->Select(ograph, pos)
   valid = where(obj_valid(atoms), n_atoms)

   if n_atoms eq 0 then begin
      printf, lun, FORMAT='(%"%s: no atoms selected")', mgh_obj_string(self)
      return
   endif

   atoms = atoms[valid]

   for j=0,n_atoms-1 do begin
      atom = atoms[j]
      status = self->PickData(ograph, atom, pos, data)
      case (atom eq self.plane) of
         0: begin
            printf, lun, FORMAT='(%"%s: atom %s, success: %d, value: %f %f %f")', $
                    mgh_obj_string(self), mgh_obj_string(atom,/SHOW_NAME), $
                    status, double(data)
         end
         1: begin
            ;; If the selected atom is the density plane, report the
            ;; data value at the selected location.

            ;; Locate the selection point in the index space of the
            ;; density planes' pixel vertices.
            self.plane->GetProperty, $
                 DATAX=datax, DATAY=datay, DATA_VALUES=data_values, STYLE=style
            xy2d = size(datax, /N_DIMENSIONS) eq 2
            case xy2d of
               0: begin
                  ii = mgh_locate(datax, XOUT=data[0])
                  jj = mgh_locate(datay, XOUT=data[1])
               end
               1: begin
                  loc = mgh_locate2a(datax, datay, XOUT=data[[0]], YOUT=data[[1]], $
                                     MISSING=-1)
                  ii = reform(loc[0,*,*])
                  jj = reform(loc[1,*,*])
                  mgh_undefine, loc
               end
            endcase

            ;; If style is 0, allow for offset of data locations
            ;; (pixel centres) and use nearest-neighbour interpolation
            if style eq 0 then begin
               ii = round(ii-0.5)
               jj = round(jj-0.5)
            endif
            ;; Interpolate & report
            value = mgh_interpolate(data_values, ii, jj, GRID=(xy2d eq 0), $
                                    MISSING=!values.f_nan)
            printf, lun, FORMAT='(%"%s: atom %s, success: %d, value: %f %f %f")', $
                    mgh_obj_string(self), mgh_obj_string(atom,/SHOW_NAME), $
                    status, double(data)
            printf, lun, FORMAT='(%"%s: atom %s, index: %f %f, value: %f")', $
                    mgh_obj_string(self), mgh_obj_string(atom,/SHOW_NAME), $
                    ii,jj,value
         end
      endcase
   endfor

end

pro Mgh_Roms_Plot_Hstats__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, {Mgh_Roms_Plot_Hstats, inherits MGH_Window, $
                 history_file: obj_new(), $
                 variable: '', parameter: '', $
                 lonlat: !false, map_structure: ptr_new(), $
                 bar: obj_new(), plane: obj_new(), contour: obj_new()}

end
