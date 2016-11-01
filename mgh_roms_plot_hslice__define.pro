;+
; CLASS NAME:
;  MGH_ROMS_PLot_Hslice
;
; PURPOSE:
;   This class generates and displays a single graph showing a Hslice (as
;   implemented in the MGHromsHistory class) through a ROMS 2D or 3D output
;   field.
;
;   Cf MGH_ROMS_Movie_Hslice, which shows an animated sequence of Hslices,
;   and MGH_ROMS_Plot_Hstats, which shows a statistical quantity derived from
;   a sequence of Hslices.
;
; CATEGORY:
;   Ocean models.
;   Object graphics
;
; CALLING SEQUENCE:
;   mgh_new, 'mgh_roms_plot_hslice', history, variable
;
; POSITIONAL PARAMETERS:
;   history
;     A reference to a ROMS history sequence object or a string arrary
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
;   variable
;     The name of a 2-D or 3-D variable in the netCDF file.
;
; KEYWORD PARAMETERS:
;   DATA_MULTIPLIER (input, numeric scalar)
;     Number by which data values are multiplied before they are loaded into
;     the density surface. Default depends on the variable and is calculated by
;     MGH_ROMS_RESOLVE_DATA
;
;   DATA_RANGE (input, numeric 2-element vector)
;     Data range for the density surface. Default depends on the variable and
;     is calculated by MGH_ROMS_RESOLVE_DATA
;
;   DATA_TRANSFORMATION (input, numeric string)
;     Function applied to data values before they are loaded into the density
;     surface.
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
;     sigma values of a constant-sigma surface on which data are to be
;     plotted. This keyword should be specified only for variables having
;     a depth coordinate and it cannot be used together with DEPTH or LEVEL.
;
;   XI_RANGE
;   ETA_RANGE
;     Use these keywords to display a subset of the domain. They are
;     interpreted relative to the rho grid.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2016-01:
;     Written based on MGH_ROMS_Movie_Hslice.
;-
function mgh_roms_plot_hslice::Init, $
     history, variable, $
     BYTE_RANGE=byte_range, $
     DATA_MULTIPLIER=data_multiplier, $
     DATA_RANGE=data_range, LOGARITHMIC=logarithmic, $
     DATA_TRANSFORMATION=data_transformation, $
     DEPTH=depth, LAYER=layer, LEVEL=level, SIGMA=sigma, $
     MAP_STRUCTURE=map_structure, MASK_VALUE=mask_value, $
     RECORD=record, STYLE=style, $
     SHOW_COLORBAR=show_colorbar, SHOW_CONTOUR=show_contour, SHOW_PLANE=show_plane, $
     SHOW_TITLE=show_title, SHOW_TIME=show_time, TITLE=title, DT_FORMAT=dt_format, $
     USE_BATH=use_bath, USE_ZETA=use_zeta, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
     X_RANGE=x_range, Y_RANGE=y_range, $
     COLORBAR_PROPERTIES=colorbar_properties, $
     CONTOUR_PROPERTIES=contour_properties, $
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

   if n_elements(variable) eq 0 then variable = 'zeta'

   if n_elements(variable) gt 1 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'variable'

   if size(variable, /TNAME) ne 'STRING' then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'variable'

   self.variable = variable

   ;; Other defaults

   mgh_roms_resolve_data, self.variable, $
      DATA_RANGE=data_range, DATA_MULTIPLIER=data_multiplier

   if n_elements(use_zeta) eq 0 then use_zeta = 0B

   if n_elements(show_colorbar) eq 0 then show_colorbar = 1B

   if n_elements(show_contour) eq 0 then show_contour = 0B

   if n_elements(show_plane) eq 0 then show_plane = 1B

   if n_elements(map_structure) gt 0 then $
      self.map_structure = ptr_new(map_structure)

   use_map_structure = ptr_valid(self.map_structure)

   if n_elements(show_title) eq 0 then show_title = !true

   if keyword_set(show_title) then begin
      if n_elements(title) eq 0 then title = ''
      if n_elements(show_time) eq 0 then show_time = 1S
      if n_elements(dt_format) eq 0 then $
         dt_format = '(C(CYI4.4,"-",CMOI2.2,"-",CDI2.2))'
   endif

   ;; Set ETA_RANGE and XI_RANGE relative to the rho grid.
   ;; The default is to show all points.

   dim_rho = [ohis->DimInfo('xi_rho', /DIMSIZE),ohis->DimInfo('eta_rho', /DIMSIZE)]

   if n_elements(xi_range) eq 0 then xi_range = [0,dim_rho[0]-1]
   if n_elements(eta_range) eq 0 then eta_range = [0,dim_rho[1]-1]

   if xi_range[0] lt 0 then xi_range[0] += dim_rho[0]
   if xi_range[1] lt 0 then xi_range[1] += dim_rho[0]

   if eta_range[0] lt 0 then eta_range[0] += dim_rho[1]
   if eta_range[1] lt 0 then eta_range[1] += dim_rho[1]

   ;; Get grid data required for horizontal slice retrievals

   var_xi_range = xi_range
   var_eta_range = eta_range

   var_dims = ohis->VarDims(variable)

   vdh0 = var_dims.horizontal[0]
   if vdh0 eq 'xi_u' || vdh0 eq 'xi_psi' then var_xi_range += [0,-1]

   vdh1 = var_dims.horizontal[1]
   if vdh1 eq 'eta_v' || vdh1 eq 'eta_psi' then var_eta_range += [0,-1]

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

   ;; Establish record to be plotted (if applicable). Establish name of
   ;; time variable--the logic here is a bit shaky, as conventions
   ;; for defining time dimensions and variables in ROMS files are
   ;; not entirely consistent.

   has_time = strlen(grid.dims.time) gt 0

   if has_time then begin
      if n_elements(record) eq 0 then record = 0
      case !true of
         ohis->HasVar(grid.dims.time): $
            time_var = grid.dims.time
         ohis->HasVar('ocean_time'): $
            time_var = 'ocean_time'
         ohis->HasVar('scrum_time'): $
            time_var = 'scrum_time'
         else: $
            message, 'Time variable not found'
      endcase
   endif

   ;; Get time units from the history file

   if has_time then begin
      if ohis->HasAtt(time_var, 'units') then begin
         time_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
      endif else begin
         time_units = {scale: 1, offset: 0}
      endelse
   endif

   ;; Get data

   if has_time then begin
      data = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, RECORD=record, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
   endif else begin
      data = ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
   endelse
   if n_elements(data_multiplier) gt 0 then $
      data *= data_multiplier
   if n_elements(data_transformation) gt 0 then $
      data = call_function(data_transformation, data)

   ;; Default graph aspect ratio, can be overridden via GRAPH_PROPERTIES

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(x_var)
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(y_var)

   aspect = mgh_aspect(x_range, y_range, LONLAT=self.lonlat && (~ use_map_structure))
   aspect = (aspect > 0.4) < 1.5

   ;; Create base graph

   xmargin = show_colorbar ? [0.375,0.4] : [0.375,0.15]

   ograph = obj_new('MGHgrGraph2D', ASPECT=aspect, XMARGIN=xmargin, $
                    NAME='ROMS horizontal slice plot', $
                    _STRICT_EXTRA=graph_properties)

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Add mask around plot area

   ograph->NewMask

   ;; Add title

   if keyword_set(show_title) then begin
      if has_time && show_time gt 0 then begin
         case show_time of
            1: begin
               time = ohis->VarGet(time_var, OFFSET=record, COUNT=1, AUTOSCALE=0)*time_units.scale
               ttt = string(FORMAT='(%"%0.3f days")', time)
            end
            2: begin
               time = ohis->VarGet(time_var, OFFSET=record, COUNT=1, AUTOSCALE=0)*time_units.scale
               ttt = string(FORMAT='(%"%s")', mgh_dt_string(time+time_units.offset, FORMAT=dt_format))
            end
            3: begin
               time = ohis->VarGet(time_var, OFFSET=record, COUNT=1, AUTOSCALE=0)*time_units.scale
               ttt = string(FORMAT='(%"%s (%0.3f days)")', mgh_dt_string(time+time_units.offset, FORMAT=dt_format), time)
            end
            4: begin
               time = ohis->VarGet('time_range', OFFSET=[0,record], COUNT=[2,1], AUTOSCALE=0)*time_units.scale
               ttt = string(FORMAT='(%"%0.3f to %0.3f days")', time)
            end
            5: begin
               time = ohis->VarGet('time_range', OFFSET=[0,record], COUNT=[2,1], AUTOSCALE=0)*time_units.scale
               ttt = string(FORMAT='(%"%s to %s")', mgh_dt_string(time+time_units.offset, FORMAT=dt_format))
            end
         endcase
         if strlen(title) gt 0 then $
            ttt = string(FORMAT='(%"%s: %s")', title, ttt)
         ograph->NewTitle, ttt
      endif else begin
         ograph->NewTitle, title
      endelse
   endif

   ;; Draw axes

   if grid.lonlat && (~ use_map_structure) then begin
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
         NAME='Land mask', STYLE=0, DEFAULT_COLOR=mgh_color('grey'), $
         ZVALUE=-2*deltaz, MISSING_POINTS=mask_missing, $
         DATAX=mgh_stagger(x_rho, DELTA=[1,1]), $
         DATAY=mgh_stagger(y_rho, DELTA=[1,1]), $
         _STRICT_EXTRA=land_properties

      ;; Just playing really
      if oland->QueryProperty('STYLE') then $
           oland->SetPropertyAttribute, 'STYLE', SENSITIVE=0

   endif

   ;; Create the palette and add a colour bar

   ograph->NewPalette, 'Matlab Jet', RESULT=palette, $
      _STRICT_EXTRA=palette_properties

   ograph->NewColorBar, RESULT=obar, FONT=ograph->GetFont(POS=1), $
        DATA_RANGE=data_range, HIDE=(~ show_colorbar), $
        LOGARITHMIC=logarithmic, PALETTE=palette, $
        CONTOUR_PROPERTIES=contour_properties, $
        SHOW_CONTOUR=show_contour, $
        _STRICT_EXTRA=colorbar_properties
   self.bar = obar

   ograph->NewAtom, 'MGHgrDensityPlane', RESULT=oplane, $
      HIDE=(~ show_plane), $
      DATAX=x_var, DATAY=y_var, DATA_VALUES=data, $
      STYLE=style, ZVALUE=-5*deltaz , NAME='Data plane', $
      COLORSCALE=obar, /STORE_DATA
   self.plane = oplane

   ograph->NewAtom, 'IDLgrContour', RESULT=ocont, $
      HIDE=(~ show_contour), $
      GEOMZ=deltaz, /PLANAR, DATA=data, GEOMX=grid.x, GEOMY=grid.y, $
      _STRICT_EXTRA=contour_properties
   self.contour = ocont

   ;; Load graph into window

   ma = ['Magnify','Translate','Context']
   ok = self->MGH_Window::Init(GRAPHICS_TREE=ograph, MOUSE_ACTION=ma, _STRICT_EXTRA=extra)
   if ~ ok then $
      message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   ;; Finalise plot & return

   self->Finalize, 'Mgh_Roms_Plot_Hslice'

   return, 1

end

; mgh_roms_plot_hslice::Cleanup
;
pro mgh_roms_plot_hslice::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Window::Cleanup

end

; mgh_roms_plot_hslice::GetProperty
;
pro mgh_roms_plot_hslice::GetProperty, $
     ALL=all, BAR=bar, BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
     HISTORY_FILE=history_file, LONLAT=lonlat, MAP_STRUCTURE=map_structure, $
     PALETTE=palette, STYLE=style, VARIABLE=variable, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::GetProperty, ALL=all, _STRICT_EXTRA=extra

   bar = self.bar

   history_file = self.history_file

   lonlat = self.lonlat

   if arg_present(all) || arg_present(map_structure) then $
        map_structure = ptr_valid(self.map_structure) ? *self.map_structure : -1

   variable = self.variable

   self.plane->GetProperty, STYLE=style

   self.bar->GetProperty, $
        BYTE_RANGE=byte_range, DATA_RANGE=data_range, PALETTE=palette

   if arg_present(all) then $
        all = create_struct(all, $
                            'bar', bar, 'byte_range', byte_range, $
                            'data_range', data_range, $
                            'history_file', history_file, $
                            'lonlat', lonlat, 'map_structure', map_structure, $
                            'palette', palette, 'style', style, $
                            'variable', variable)

end

; mgh_roms_plot_hslice::SetProperty
;
pro mgh_roms_plot_hslice::SetProperty, $
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

; mgh_roms_plot_hslice::About
;
pro mgh_roms_plot_hslice::About, lun

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

; mgh_roms_plot_hslice::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro mgh_roms_plot_hslice::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::BuildMenuBar

   obar = mgh_widget_self(self.menu_bar)

   obar->NewItem, PARENT='Tools', SEPARATOR=[1,0,0,0,1], MENU=[1,0,0,0,0], $
        ['Data Range','Edit Palette...','Set Style...', $
         'View Colour Scale...','View Data Values...']

   obar->NewItem, PARENT='Tools.Data Range', ['Set...','Fit Data']

end


; mgh_roms_plot_hslice::EventMenuBar
;
function mgh_roms_plot_hslice::EventMenuBar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

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

; mgh_roms_plot_hslice::ExportData
;
pro mgh_roms_plot_hslice::ExportData, values, labels

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

; mgh_roms_plot_hslice::PickReport
;
pro mgh_roms_plot_hslice::PickReport, pos, LUN=lun

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

pro mgh_roms_plot_hslice__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, {mgh_roms_plot_hslice, inherits MGH_Window, $
                 history_file: obj_new(), $
                 variable: '', lonlat: !false, map_structure: ptr_new(), $
                 bar: obj_new(), plane: obj_new(), contour: obj_new()}

end
