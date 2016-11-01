;+
; CLASS NAME:
;   Mgh_Roms_Plot_Grid
;
; PURPOSE:
;   This class generates a graph of the horizontal grid from a ROMS file.
;
; OBJECT CREATION SEQUENCE
;   obj = obj_new('Mgh_Roms_Plot_Grid', File)
;
; INPUTS:
;   Files
;     The name of a ROMS file (or list of files) containing 3D grid data.
;
; KEYWORD PARAMETERS:
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2000-01:
;     Written.
;   Mark Hadfield, 2009-11:
;     Added LAND_PROPERTIES keyword.
;   Mark Hadfield, 2008-11:
;     Added support for map projections via the MAP_STRUCTURE
;     keyword.
;   Mark Hadfield, 2011-08:
;     Removed the history_destroy and station_destroy tags from
;     the class structure.
;   Mark Hadfield, 2014-07:
;     Longitude & latitude labels on axes now have two decimal places rather
;     than one.
;-
function MGH_ROMS_plot_grid::Init, history, station, $
     BATHYMETRY_VERSION=bathymetry_version, $
     BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
     DATA_MULTIPLIER=data_multiplier, $
     HISTORY_FILE=history_file, $
     MAP_STRUCTURE=map_structure, $
     STATION_FILE=station_file, VARIaBLE=variable, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
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

   if n_elements(history_file) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'history_file'

   case size(history_file, /TNAME) of
      'STRING': begin
         self.history_file = obj_new('MGHromsHistory', history_file)
      end
      'OBJREF': begin
         self.history_file = history_file
      end
      else: message, 'The history_file argument is of the wrong data type'
   endcase
   ohis = self.history_file

   ;; Process station_file argument

   if n_elements(station_file) eq 0 && n_elements(station) gt 0 then $
        station_file = station

   case size(station_file, /TNAME) of
      'STRING': begin
         self.station_file = obj_new('MGHromsStation', station_file)
      end
      'OBJREF': begin
         self.station_file = station_file
      end
      'UNDEFINED':              ; Do nothing
      else: message, 'The station_file argument is of the wrong data type'
   endcase
   osta = self.station_file

   ;; Select variable to be plotted

   self.variable = n_elements(variable) eq 1 ? variable : 'h'

   mgh_roms_resolve_data, self.variable, $
        DATA_MULTIPLIER=data_multiplier

   ;; Other defaults

   if n_elements(show_bar) eq 0 then show_bar = 1B

   if n_elements(show_contour) eq 0 then show_contour = 0B

   if n_elements(map_structure) gt 0 then $
        self.map_structure = ptr_new(map_structure)

   use_map_structure = ptr_valid(self.map_structure)

   ;; Set ETA_RANGE and XI_RANGE relative to the rho grid.
   ;; The default is to show all points (cf. all interior points
   ;; for the MGH_ROMS_Movie_Hslice class.

   dim_rho = [ohis->DimInfo('xi_rho', /DIMSIZE), $
              ohis->DimInfo('eta_rho', /DIMSIZE)]

   if n_elements(xi_range) eq 0 then xi_range = [0,dim_rho[0]-1]
   if n_elements(eta_range) eq 0 then eta_range = [0,dim_rho[1]-1]

   if xi_range[0] lt 0 then xi_range[0] += dim_rho[0]
   if xi_range[1] lt 0 then xi_range[1] += dim_rho[0]

   if eta_range[0] lt 0 then eta_range[0] += dim_rho[1]
   if eta_range[1] lt 0 then eta_range[1] += dim_rho[1]

   ;; Get dimensions for the variable to be plotted.
   ;; Check out consistency of keywords & provide defaults.

   var_dims = ohis->VarDimNames(self.variable)

   nvdims = n_elements(var_dims)

   if nvdims lt 2 || nvdims gt 3 then $
        message, 'The bathymetry variable must have 3 or 4 dimensions'

   if var_dims[0] ne 'xi_rho' then $
        message, 'The 1st dimension of the bathymetry variable must be xi_rho'
   if var_dims[1] ne 'eta_rho' then $
        message, 'The 2nd dimension of the bathymetry variable must be eta_rho'
   if nvdims eq 3 then begin
      if var_dims[2] ne 'bath' then begin
         message, 'The final dimension of the bathymetry variable must be bath'
      endif
      if n_elements(bathymetry_version) eq 0 then begin
         bathymetry_version = ohis->DimInfo('bath', /DIMSIZE) - 1
      endif
   endif

   ;; Extract a list of CPP options

   options = ohis->GetCppOptions()

   ;; Are the grid locations in (x,y) or (lon,lat) coordinates?

   self.lonlat = ohis->HasVar('lon_rho') && ohis->HasVar('lat_rho')

   ;; Read x & y positions of rho points

   offset = [xi_range[0],eta_range[0]]
   count = [xi_range[1]-xi_range[0]+1,eta_range[1]-eta_range[0]+1]

   x_rho = self.lonlat ? $
        ohis->VarGet('lon_rho', OFFSET=offset, COUNT=count) : $
        ohis->VarGet('x_rho', OFFSET=offset, COUNT=count)
   y_rho = self.lonlat ? $
        ohis->VarGet('lat_rho', OFFSET=offset, COUNT=count) : $
        ohis->VarGet('y_rho', OFFSET=offset, COUNT=count)

   if use_map_structure then begin
      if ~ self.lonlat then $
           message, 'Cannot use map structure without lon, lat data'
      xy = map_proj_forward(x_rho, y_rho, MAP_STRUCTURE=*self.map_structure)
      x_rho = reform(xy[0,*], size(x_rho, /DIMENSIONS))
      y_rho = reform(xy[1,*], size(y_rho, /DIMENSIONS))
      mgh_undefine, xy
   endif

   ;; Set default graph aspect ratio

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(x_rho)
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(y_rho)

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
      xap = {tickformat: 'mgh_tf_longitude', tickfrmtdata: {format:'(F10.2)'}}
      yap = {tickformat: 'mgh_tf_latitude', tickfrmtdata: {format:'(F10.2)'} }
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

   ;; Add walls (but only if that boundary is included in the plot).
   ;; For each wall we extract the 2 rows/columns on each
   ;; side of the physical boundary from the x_rho & y_rho arrays into
   ;; variables xr & yr, call mgh_stagger to get xw & yw along the
   ;; wall, then trim the interior points off.

   walls = ohis->GetWalls()

   if walls[0] && xi_range[0] eq 0 then begin
      ;; Western wall
      xw = mgh_stagger(x_rho[0:1,*], DELTA=[1,1])
      yw = mgh_stagger(y_rho[0:1,*], DELTA=[1,1])
      xw = xw[0:1,*]
      yw = yw[0:1,*]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', $
           DATAX=temporary(xw), DATAY=temporary(yw), DATAZ=temporary(zw), $
           STYLE=2, COLOR=[127,127,127]
   endif

   if walls[1] && eta_range[0] eq 0 then begin
      ;; Southern wall
      xw = mgh_stagger(x_rho[*,0:1], DELTA=[1,1])
      yw = mgh_stagger(y_rho[*,0:1], DELTA=[1,1])
      xw = xw[*,0:1]
      yw = yw[*,0:1]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', $
           DATAX=temporary(xw), DATAY=temporary(yw), DATAZ=temporary(zw), $
           STYLE=2, COLOR=[127,127,127]
   endif

   if walls[2] && xi_range[1] eq (dim_rho[0]-1) then begin
      ;; Eastern wall
      xw = mgh_stagger(x_rho[-2:-1,*], DELTA=[1,1])
      yw = mgh_stagger(y_rho[-2:-1,*], DELTA=[1,1])
      xw = xw[1:2,*]
      yw = yw[1:2,*]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', $
           DATAX=temporary(xw), DATAY=temporary(yw), DATAZ=temporary(zw), $
           STYLE=2, COLOR=[127,127,127]
   endif

   if walls[3] && eta_range[1] eq (dim_rho[1]-1) then begin
      ;; Northern wall
      xw = mgh_stagger(x_rho[*,-2:-1], DELTA=[1,1])
      yw = mgh_stagger(y_rho[*,-2:-1], DELTA=[1,1])
      xw = xw[*,1:2]
      yw = yw[*,1:2]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', $
           DATAX=temporary(xw), DATAY=temporary(yw), DATAZ=temporary(zw), $
           STYLE=2, COLOR=[127,127,127]
   endif

   ;; Add mask

   if ohis->HasVar('mask_rho') then begin

      mask_rho = round(ohis->VarGet('mask_rho', OFFSET=offset, COUNT=count))

      ograph->NewAtom, 'MGHgrColorPlane', NAME='Land mask', $
           DEFAULT_COLOR=!color.light_grey, DEPTH_OFFSET=1, $
           MISSING_POINTS=mask_rho, $
           DATAX=mgh_stagger(x_rho, DELTA=[1,1]), $
           DATAY=mgh_stagger(y_rho, DELTA=[1,1]), $
           _STRICT_EXTRA=land_properties

   endif

   ograph->NewPalette, 'Matlab Jet', RESULT=palette, $
        _STRICT_EXTRA=palette_properties

   ;; Retrieve and plot variable values

   case nvdims of
      2: begin
         vardata = ohis->VarGet(self.variable, OFFSET=offset, COUNT=count, /AUTOSCALE)
      end
      3: begin
         vardata = ohis->VarGet(self.variable, /AUTOSCALE, $
                                OFFSET=[offset,bathymetry_version], $
                                COUNT=[count,1])
      end
   endcase

   vardata *= data_multiplier

   if n_elements(data_range) eq 0 then begin
      data_range = mgh_minmax(vardata)
      if data_range[0] eq data_range[1] then $
           data_range += [-1,1]
   endif

   ograph->NewAtom, 'MGHgrDensityPlane', STYLE=0, DEPTH_OFFSET=1, $
        DATAX=mgh_stagger(x_rho, DELTA=[1,1]), $
        DATAY=mgh_stagger(y_rho, DELTA=[1,1]), $
        DATA_VALUES=vardata, $
        BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
        PALETTE=palette, NAME='Data', /STORE_DATA, ZVALUE=-20*deltaz, RESULT=osurf
   self.plane = osurf

   ograph->NewAtom, 'IDLgrContour', RESULT=ocont, $
        HIDE=(~ show_contour), COLOR=mgh_color('white'), $
        GEOMZ=-10*deltaz, /PLANAR, DATA=vardata, GEOMX=x_rho, GEOMY=y_rho, $
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

   ;; Plot stations

   if obj_valid(osta) then begin
      nsta = osta->DimInfo('station', /DIMSIZE)
      if self.lonlat then begin
         xs = osta->VarGet('lon_rho')
         ys = osta->VarGet('lat_rho')
         if use_map_structure then begin
            xy = map_proj_forward(xs, ys, MAP_STRUCTURE=*self.map_structure)
            if size(xy, /N_DIMENSIONS) eq 1 then begin
               xs = xy[[0]]
               ys = xy[[1]]
            endif else begin
               xs = reform(xy[0,*], size(xs, /DIMENSIONS))
               ys = reform(xy[1,*], size(ys, /DIMENSIONS))
            endelse
            mgh_undefine, xy
         endif
      endif else begin
         xs = osta->VarGet('x_rho')
         ys = osta->VarGet('y_rho')
      endelse
      zs = replicate(-5*deltaz, nsta)
      ograph->NewText, strtrim(sindgen(nsta), 2), $
           ALIGN=0.5, VERTICAL_ALIGN=0.5, COLOR=mgh_color('black'), $
           LOCATIONS=transpose([[xs],[ys],[zs]]), NAME='Stations'
   endif

   ;; Initialise window

   a = ['Magnify','Pick','Context']

   ok = self->MGH_Window::Init(ograph, MOUSE_ACTION=a, _STRICT_EXTRA=extra)
   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   self->Finalize, 'Mgh_Roms_Plot_Grid'

   return, 1

end

; Mgh_Roms_Plot_Grid::GetProperty
;
pro Mgh_Roms_Plot_Grid::GetProperty, $
     ALL=all, HISTORY_FILE=history_file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR

   self->MGH_Window::GetProperty, ALL=all, _STRICT_EXTRA=extra

   history_file = self.history_file

   if arg_present(all) then $
        all = create_struct(all, 'history_file', history_file)

end

; Mgh_Roms_Plot_Grid::About
;
pro Mgh_Roms_Plot_Grid::About, lun

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

; Mgh_Roms_Plot_Grid::Cleanup
;
pro Mgh_Roms_Plot_Grid::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Window::Cleanup

end

; Mgh_Roms_Plot_Grid::SetProperty
;
pro Mgh_Roms_Plot_Grid::SetProperty, $
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

; Mgh_Roms_Plot_Grid::GetProperty
;
pro Mgh_Roms_Plot_Grid::GetProperty, $
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

   station_file = self.station_file

   variable = self.variable

   self.plane->GetProperty, BYTE_RANGE=byte_range, $
        DATA_RANGE=data_range, PALETTE=palette

   self->MGH_Window::GetProperty, _STRICT_EXTRA=extra

end


; Mgh_Roms_Plot_Grid::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Plot_Grid::BuildMenuBar

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

; Mgh_Roms_Plot_Grid::EventMenuBar
;
function Mgh_Roms_Plot_Grid::EventMenuBar, event

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

; Mgh_Roms_Plot_Grid::ExportData
;
pro Mgh_Roms_Plot_Grid::ExportData, values, labels

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

; Mgh_Roms_Plot_Grid::PickReport
;
pro Mgh_Roms_Plot_Grid::PickReport, pos, LUN=lun

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

pro Mgh_Roms_Plot_Grid__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, {Mgh_Roms_Plot_Grid, inherits MGH_Window, $
                 history_file: obj_new(), station_file: obj_new(), $
                 lonlat: !false, map_structure: ptr_new(), $
                 variable: '', bar: obj_new(), $
                 plane: obj_new(), contour: obj_new()}

end
