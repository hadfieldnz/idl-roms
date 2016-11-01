;+
; CLASS NAME:
;   Mgh_Roms_Plot_Grid_Diff
;
; PURPOSE:
;   This class generates a graph of the horizontal grid from a ROMS file.
;
; OBJECT CREATION SEQUENCE
;   obj = obj_new('mgh_roms_plot_grid_diff', File)
;
; INPUTS:
;   Files
;     The name of a ROMS file (or list of files) containing grid data.
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2013-02:
;     Written.
function mgh_roms_plot_grid_diff::Init, history, $
     BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
     DATA_MULTIPLIER=data_multiplier, $
     MULTIPLIER=multiplier, $
     HISTORY_FILE=history_file, $
     MAP_STRUCTURE=map_structure, $
     VARIaBLE=variable, $
     X_RANGE=x_range, Y_RANgE=y_range, $
     SHOW_BAR=show_bar, SHOW_CONTOUR=show_contour, $
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

   ;; Process history_file argument

   if n_elements(history_file) eq 0 && n_elements(history) gt 0 then $
        history_file = history

   n_his = n_elements(history_file)

   if n_his eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'history_file'


   case size(history_file, /TNAME) of
      'STRING': begin
         self.history_file = list()
         for f=0,n_his-1 do $
              self.history_file[f] = obj_new('MGHromsHistory', history_file[f])
      end
      'OBJREF': begin
         self.history_file = list(history_file, /EXTRACT)
      end
      else: message, 'The history_file argument is of the wrong data type'
   endcase
   ohis = self.history_file

   ;; Select variable to be plotted

   self.variable = n_elements(variable) eq 1 ? variable : 'h'

   ;; Other defaults

   if n_elements(data_multiplier) eq 0 then data_multiplier = 1

   if n_elements(multiplier) eq 0 then multiplier = replicate(1.0/n_his, n_his)

   if n_elements(show_bar) eq 0 then show_bar = 1B

   if n_elements(show_contour) eq 0 then show_contour = 0B

   if n_elements(map_structure) gt 0 then $
        self.map_structure = ptr_new(map_structure)

   use_map_structure = ptr_valid(self.map_structure)

   ;; Are the grid locations in (x,y) or (lon,lat) coordinates?

   self.lonlat = ohis[0]->HasVar('lon_rho') && ohis[0]->HasVar('lat_rho')

   ;; Read x & y positions of rho points

   x0 = self.lonlat ? ohis[0]->VarGet('lon_rho') : ohis[0]->VarGet('x_rho')
   y0 = self.lonlat ? ohis[0]->VarGet('lat_rho') : ohis[0]->VarGet('y_rho')

   if use_map_structure then begin
      if ~ self.lonlat then $
           message, 'Cannot use map structure without lon, lat data'
      xy = map_proj_forward(x0, y0, MAP_STRUCTURE=*self.map_structure)
      x0 = reform(xy[0,*], size(x0, /DIMENSIONS))
      y0 = reform(xy[1,*], size(y0, /DIMENSIONS))
      mgh_undefine, xy
   endif

   ;; Set default graph aspect ratio

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(x0)
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(y0)

   aspect = mgh_aspect(x_range, y_range, LONLAT=self.lonlat && (~ use_map_structure))
   aspect = (aspect > 0.2) < 1.75

   ;; Default margin size depends on whether a colour bar is shown

   xmargin = keyword_set(show_bar) ? [0.375,0.40] : [0.375,0.15]

   ;; Create graph & add fonts

   mgh_new, 'MGHgrGraph2D', RESULT=ograph, $
            ASPECT=aspect, NAME='ROMS grid plot', XMaRGIN=xmargin, $
            _STRICT_EXTRA=graph_properties

   ograph->NewMask

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Draw axes

   if self.lonlat && (~ use_map_structure) then begin
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

   ograph->NewBackground

   ;; Add mask

   if ohis[0]->HasVar('mask_rho') then begin

      mask_rho = round(ohis[0]->VarGet('mask_rho', OFFSET=offset, COUNT=count))

      datax = mgh_roms_stagger(x0, TO='psi')
      datay = mgh_roms_stagger(y0, TO='psi')

      ograph->NewAtom, 'MGHgrColorPlane', NAME='Land mask', $
           DEFAULT_COLOR=!color.light_grey, DEPTH_OFFSET=1, $
           MISSING_POINTS=mask_rho, $
           DATAX=mgh_stagger(x0, DELTA=[1,1]), $
           DATAY=mgh_stagger(y0, DELTA=[1,1]), $
           _STRICT_EXTRA=land_properties

   endif

   ograph->NewPalette, 'Matlab Jet', RESULT=palette, $
        _STRICT_EXTRA=palette_properties

   ;; Retrieve and plot data

   data = data_multiplier*multiplier[0]*ohis[0]->VarGet(self.variable, /AUTOSCALE)

   for f=1,n_his-1 do begin
      xf = self.lonlat ? ohis[f]->VarGet('lon_rho') : ohis[f]->VarGet('x_rho')
      yf = self.lonlat ? ohis[f]->VarGet('lat_rho') : ohis[f]->VarGet('y_rho')
      if use_map_structure then begin
         xy = map_proj_forward(xf, yf, MAP_STRUCTURE=*self.map_structure)
         xf = reform(xy[0,*], size(x0, /DIMENSIONS))
         yf = reform(xy[1,*], size(y0, /DIMENSIONS))
         mgh_undefine, xy
      endif
      loc = mgh_locate2(xf, yf, XOUT=x0, YOUT=y0)
      dataf = data_multiplier*multiplier[f]*ohis[f]->VarGet(self.variable, /AUTOSCALE)
      data += mgh_interpolate(dataf, reform(loc[0,*,*]), reform(loc[1,*,*]))
   endfor

   if n_elements(data_range) eq 0 then begin
      data_range = mgh_minmax(data)
      if data_range[0] eq data_range[1] then $
           data_range += [-1,1]
   endif

   ograph->NewAtom, 'MGHgrDensityPlane', STYLE=0, DEPTH_OFFSET=1, $
        DATAX=mgh_stagger(x0, DELTA=[1,1]), $
        DATAY=mgh_stagger(y0, DELTA=[1,1]), $
        DATA_VALUES=data, $
        BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
        PALETTE=palette, NAME='Data', /STORE_DATA, ZVALUE=-20*deltaz, RESULT=osurf
   self.plane = osurf

   ograph->NewAtom, 'IDLgrContour', RESULT=ocont, $
        HIDE=(~ show_contour), COLOR=mgh_color('white'), $
        GEOMZ=-10*deltaz, /PLANAR, DATA=data, GEOMX=x0, GEOMY=y0, $
        _STRICT_EXTRA=contour_properties
   self.contour = ocont

   ;; Add a colour bar

   if keyword_set(show_bar) then begin
      ograph->NewColorBar, RESULT=obar, $
           FONT=ograph->GetFont(POS=1), COLORSCALE=osurf, $
           SHOW_CONTOUR=show_contour, CONTOUR_PROPERTIES=contour_properties, $
           _STRICT_EXTRA=colorbar_properties
      self.bar = obar
   endif

   ;; Initialise window

   a = ['Magnify','Pick','Context']

   ok = self->MGH_Window::Init(ograph, MOUSE_ACTION=a, _STRICT_EXTRA=extra)
   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   self->Finalize, 'mgh_roms_plot_grid_diff'

   return, 1

end

; mgh_roms_plot_grid_diff::GetProperty
;
pro mgh_roms_plot_grid_diff::GetProperty, $
     ALL=all, HISTORY_FILE=history_file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR

   self->MGH_Window::GetProperty, ALL=all, _STRICT_EXTRA=extra

   history_file = self.history_file

   if arg_present(all) then $
        all = create_struct(all, 'history_file', history_file)

end

; mgh_roms_plot_grid_diff::About
;
pro mgh_roms_plot_grid_diff::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::About, lun

   self->GetProperty, HISTORY_FILE=history_file

   if obj_valid(history_file) then begin
      printf, lun, FORMAT='(%"%s: my history file sequence is %s")', $
           mgh_obj_string(self), mgh_obj_string(history_file)
      history_file->GetProperty, FILE_NAME=file_name
      printf, lun, self, ': its files are:', file_name
   endif

end

; mgh_roms_plot_grid_diff::Cleanup
;
pro mgh_roms_plot_grid_diff::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Window::Cleanup

end

; mgh_roms_plot_grid_diff::SetProperty
;
pro mgh_roms_plot_grid_diff::SetProperty, $
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

; mgh_roms_plot_grid_diff::GetProperty
;
pro mgh_roms_plot_grid_diff::GetProperty, $
     BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
     HISTORY_FILE=history_file, LONLAT=lonlat, PALETTE=palette, $
     STATION_FILE=station_file, VARIaBLE=variable, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   history_file = self.history_file

   lonlat = self.lonlat

   variable = self.variable

   self.plane->GetProperty, BYTE_RANGE=byte_range, $
        DATA_RANGE=data_range, PALETTE=palette

   self->MGH_Window::GetProperty, _STRICT_EXTRA=extra

end


; mgh_roms_plot_grid_diff::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro mgh_roms_plot_grid_diff::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::BuildMenuBar

   obar = mgh_widget_self(self.menu_bar)

   obar->NewItem, PARENT='Tools', SEPARATOR=[1,0,0,1], MENU=[1,0,0,0], $
        ['Data Range','Edit Palette...', $
         'View Colour Scale...','View Data Values...']

   obar->NewItem, PARENT='Tools.Data Range', ['Set...','Fit Data']

end

; mgh_roms_plot_grid_diff::EventMenuBar
;
function mgh_roms_plot_grid_diff::EventMenuBar, event

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

;      'TOOLS.SET DATA RANGE': begin
;         mgh_new, 'MGH_GUI_SetArray', CAPTION='Range', CLIENT=self, $
;                  /FLOATING, GROUP_LEADER=self.base, IMMEDIATE=0, $
;                  N_ELEMENTS=2, PROPERTY_NAME='DATA_RANGE'
;         return, 0
;      end

      'TOOLS.EDIT PALETTE': begin
         self->GetProperty, PALETTE=palette
         mgh_new, 'MGH_GUI_Palette_Editor', palette, CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, /IMMEDIATE
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
         xvaredit, reform(data_values), GROUP=self.base, $
                   X_SCROLL_SIZE=(data_dims[0] < 12), $
                   Y_SCROLL_SIZE=(data_dims[1] < 30)
         return, 0
      end

      else: return, self->MGH_Window::EventMenubar(event)

   endcase

end

; mgh_roms_plot_grid_diff::ExportData
;
pro mgh_roms_plot_grid_diff::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::ExportData, values, labels

   labels = [labels, 'History-File Object', 'Station-File Object']
   values = [values, ptr_new(self.history_file), ptr_new(self.station_file)]

   self.plane->GetProperty, DATA_VALUES=data_values
   if n_elements(data_values) gt 0 then begin
      labels = [labels, 'Data Values']
      values = [values, ptr_new(data_values)]
   endif

end

; mgh_roms_plot_grid_diff::PickReport
;
pro mgh_roms_plot_grid_diff::PickReport, pos, LUN=lun

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
               ii = round(temporary(ii)-0.5)
               jj = round(temporary(jj)-0.5)
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

pro mgh_roms_plot_grid_diff__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, {mgh_roms_plot_grid_diff, inherits MGH_Window, $
                 history_file: obj_new(), $
                 lonlat: 0B, map_structure: ptr_new(), $
                 variable: '', bar: obj_new(), $
                 plane: obj_new(), contour: obj_new()}

end
