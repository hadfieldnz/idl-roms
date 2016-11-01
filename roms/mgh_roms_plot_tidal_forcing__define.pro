;+
; CLASS NAME:
;   Mgh_Roms_Plot_Tidal_Forcing
;
; PURPOSE:
;   This class generates a graph of tidal forcing data from a ROMS forcing
;   file.
;
; OBJECT CREATION SEQUENCE
;   obj = obj_new('Mgh_Roms_Plot_Tidal_Forcing', File)
;
; INPUTS:
;   Files
;     The name of a ROMS file (or list of files) containing tidal forcing data
;
; KEYWORD PARAMETERS:
;
;   MAP_STRUCTURE
;     If a map structure is supplied via this keyword, and if there
;     are longitude and latitude data in the file, then the longitudes
;     and latitudes are converted with MAP_PROJ_FORWARD before display.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, ????:
;     Written.
;   Mark Hadfield, 2010-09:
;     Constituent now specified by name. This requires the constituent names
;     to have been written to the forcing file: routine
;     MGH_ROMS_FRC_LOAD_EEZ_TIDE has been modified to do this.
;-

function Mgh_Roms_Plot_Tidal_Forcing::Init, ffile, varname, $
     CONSTITUENT=constituent, $
     BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
     DATA_MULTIPLIER=data_multiplier, $
     CONTOUR_PROPERTIES=contour_properties, $
     FORCING_FILE=forcing_file, $
     MAP_STRUCTURE=map_structure, $
     SHOW_CONTOUR=show_contour, $
     GRAPH_PROPERTIES=graph_properties, $
     PALETTE_PROPERTIES=palette_properties, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process forcing_file argument

   if n_elements(forcing_file) eq 0 && n_elements(ffile) gt 0 then forcing_file = ffile

   if n_elements(forcing_file) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'forcing_file'

   case size(forcing_file, /TNAME) of
      'STRING': begin
         self.forcing_file = obj_new('MGHromsHistory', forcing_file)
         self.forcing_destroy = 1B
      end
      'OBJREF': begin
         self.forcing_file = forcing_file
         self.forcing_destroy = 0B
      end
      else: message, 'The forcing_file argument is of the wrong data type'
   endcase
   ofrc = self.forcing_file

   ;; Select variable to be plotted

   self.varname = n_elements(varname) eq 1 ? varname : 'tide_Eamp'

   mgh_roms_resolve_data, self.varname, $
        DATA_MULTIPLIER=data_multiplier

   ;; Other defaults

   if n_elements(constituent) eq 0 then $
        constituent = 'M2'

   if n_elements(show_contour) eq 0 then show_contour = 0B

   use_map_structure = n_elements(map_structure) gt 0

   ;; Get dimensions for the variable to be plotted.
   ;; Check out consistency of keywords & provide defaults.

   var_dims = ofrc->VarDimNames(self.varname)

   nvdims = n_elements(var_dims)

   if nvdims ne 3 then $
        message, 'The tidal variable must have 3 dimensions'

   if var_dims[0] ne 'xi_rho' then $
        message, 'The 1st dimension of the tidal variable must be xi_rho'
   if var_dims[1] ne 'eta_rho' then $
        message, 'The 2nd dimension of the tidal variable must be eta_rho'
   if var_dims[2] ne 'tide_period' then begin
        message, 'The final dimension of the tidal variable must be tide_period'
   endif

   ;; Locate the constituent

   tide_name = ofrc->VarGet('tide_name')

   l_con = where(strmatch(tide_name, constituent), n_match)
   if n_match eq 0 then $
        message, 'Constituent '+constituent+' not found.'
   if n_match gt 1 then $
        message, 'I was not expecting that!'

   ;; Read file dimensions

   dim_rho = ofrc->DimRho()

   ;; Are the grid locations in (x,y) or (lon,lat) coordinates?

   self.lonlat = ofrc->HasVar('lon_rho') && ofrc->HasVar('lat_rho')

   ;; Read x & y positions of rho points

   case self.lonlat of
      0: begin
         x_rho = ofrc->VarGet('x_rho')
         y_rho = ofrc->VarGet('y_rho')
      end
      1: begin
         x_rho = ofrc->VarGet('lon_rho')
         y_rho = ofrc->VarGet('lat_rho')
      end
   endcase

   ;; Convert all position data to map projection if appropriate

   if keyword_set(use_map_structure) then begin
      if ~ self.lonlat then $
           message, 'Cannot use map structure without lon, lat data'
      xy = map_proj_forward(x_rho, y_rho, MAP_STRUCTURE=map_structure)
      x_rho = reform(xy[0,*], size(x_rho, /DIMENSIONS))
      y_rho = reform(xy[1,*], size(y_rho, /DIMENSIONS))
      mgh_undefine, xy
   endif

   ;; Default graph aspect ratio, can be overridden via GRAPH_PROPERTIES

   x_range = mgh_minmax(x_rho)
   y_range = mgh_minmax(y_rho)

   aspect = double(mgh_diff(y_range))/ $
            double(mgh_diff(x_range))
   if self.lonlat && (~ use_map_structure) then $
        aspect /= cos(!dtor*mgh_avg(y_range))
   aspect = (aspect > 0.2) < 1.75

   ;; Create graph & add fonts

   mgh_new, 'MGHgrGraph2D', RESULT=ograph, $
            ASPECT=aspect, NAME='ROMS tidal forcing plot', XMARGIN=[0.30,0.45], $
            _STRICT_EXTRA=graph_properties

   ograph->NewMask

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Draw axes

   if self.lonlat && (~ use_map_structure) then begin
      ograph->NewAxis, 0, $
           RANGE=mgh_minmax(x_rho), /EXACT, /EXTEND, $
           TICKFORMAT='mgh_tf_longitude', TICKFRMTDATA={format:'(F10.1)'}, $
           _STRICT_EXTRA=xaxis_properties
      ograph->NewAxis, 1, $
           RANGE=mgh_minmax(y_rho), /EXACT, /EXTEND, $
           TICKFORMAT='mgh_tf_latitude', TICKFRMTDATA={format:'(F10.1)'}, $
           _STRICT_EXTRA=yaxis_properties
   endif else begin
      ograph->NewAxis, 0, $
           RANGE=mgh_minmax(x_rho), /EXACT, /EXTEND, TITLE='X (km)', $
           TICKFORMAT='mgh_tf_linear', TICKFRMTDATA={scale:1.E-3, format:'(F10.1)'}, $
           _STRICT_EXTRA=xaxis_properties
      ograph->NewAxis, 1, $
           RANGE=mgh_minmax(y_rho), /EXACT, /EXTEND, TITLE='Y (km)', $
           TICKFORMAT='mgh_tf_linear', TICKFRMTDATA={scale:1.E-3, format:'(F10.1)'}, $
           _STRICT_EXTRA=yaxis_properties
   endelse


   ;; Add mask

   if ofrc->HasVar('mask_rho') then begin

      mask_rho = round(ofrc->VarGet('mask_rho'))

      ograph->NewAtom, 'MGHgrColorPlane', DEFAULT_COLOR=mgh_color('grey'), $
           MISSING_POINTS=mask_rho, DEPTH_OFFSET=1, $
           DATAX=mgh_stagger(x_rho, DELTA=[1,1]), $
           DATAY=mgh_stagger(y_rho, DELTA=[1,1]), $
           NAME='Land mask'

   endif

   ograph->NewPalette, 'Matlab Jet', RESULT=palette, $
    _STRICT_EXTRA=palette_properties

   ;; Retrieve and plot tidal data

   vardata = ofrc->VarGet(self.varname, /AUTOSCALE, $
                          COUNT=[0,0,1], $
                          OFFSET=[0,0,l_con])

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
        PALETTE=palette, NAME='Data', /STORE_DATA, ZVALUE=-5*deltaz, RESULT=osurf
   self.plane = osurf

   ograph->NewAtom, 'IDLgrContour', RESULT=ocont, $
        HIDE=(~ show_contour), COLOR=mgh_color('white'), $
        GEOMZ=-2*deltaz, /PLANAR, DATA=vardata, GEOMX=x_rho, GEOMY=y_rho, $
        _STRICT_EXTRA=contour_properties
   self.contour = ocont

   ;; Add a colour bar

   ograph->NewColorBar, RESULT=obar, $
        FONT=ograph->GetFont(POS=1), COLORSCALE=osurf, $
        SHOW_CONTOUR=show_contour, CONTOUR_PROPERTIES=contour_properties
   self.bar = obar

   ;; Initialise window

   a = ['Magnify','Pick','Context']

   ok = self->MGH_Window::Init(ograph, MOUSE_ACTION=a, _STRICT_EXTRA=extra)
   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   self->Finalize, 'Mgh_Roms_Plot_Tidal_Forcing'

   return, 1

end

; Mgh_Roms_Plot_Tidal_Forcing::GetProperty
;
pro Mgh_Roms_Plot_Tidal_Forcing::GetProperty, $
     ALL=all, FORCING_FILE=forcing_file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR

   self->MGH_Window::GetProperty, ALL=all, _STRICT_EXTRA=extra

   forcing_file = self.forcing_file

   if arg_present(all) then $
        all = create_struct(all, 'forcing_file', forcing_file)

end

; Mgh_Roms_Plot_Tidal_Forcing::About
;
pro Mgh_Roms_Plot_Tidal_Forcing::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::About, lun

   self->GetProperty, FORCING_FILE=forcing_file

   if obj_valid(forcing_file) then begin
      printf, lun, FORMAT='(%"%s: my forcing file sequence is %s")', $
           mgh_obj_string(self), mgh_obj_string(forcing_file)
      forcing_file->GetProperty, FILE_NAME=file_name
      printf, lun, self, ': its files are:', file_name
   endif

end

; Mgh_Roms_Plot_Tidal_Forcing::Cleanup
;
pro Mgh_Roms_Plot_Tidal_Forcing::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if self.forcing_destroy then obj_destroy, self.forcing_file

   self->MGH_Window::Cleanup

end

; Mgh_Roms_Plot_Tidal_Forcing::SetProperty
;
pro Mgh_Roms_Plot_Tidal_Forcing::SetProperty, $
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

; Mgh_Roms_Plot_Tidal_Forcing::GetProperty
;
pro Mgh_Roms_Plot_Tidal_Forcing::GetProperty, $
     BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
     FORCING_FILE=forcing_file, LONLAT=lonlat, PALETTE=palette, $
     VARNAME=varname, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   forcing_file = self.forcing_file

   lonlat = self.lonlat

   varname = self.varname

   self.plane->GetProperty, BYTE_RANGE=byte_range, $
        DATA_RANGE=data_range, PALETTE=palette

   self->MGH_Window::GetProperty, _STRICT_EXTRA=extra

end


; Mgh_Roms_Plot_Tidal_Forcing::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Plot_Tidal_Forcing::BuildMenuBar

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

; Mgh_Roms_Plot_Tidal_Forcing::EventMenuBar
;
function Mgh_Roms_Plot_Tidal_Forcing::EventMenuBar, event

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

; Mgh_Roms_Plot_Tidal_Forcing::ExportData
;
pro Mgh_Roms_Plot_Tidal_Forcing::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::ExportData, values, labels

   labels = [labels, 'Forcing-File Object']
   values = [values, ptr_new(self.forcing_file)]

   self.plane->GetProperty, DATA_VALUES=data_values
   if n_elements(data_values) gt 0 then begin
      labels = [labels, 'Data Values']
      values = [values, ptr_new(data_values)]
   endif

end

; Mgh_Roms_Plot_Tidal_Forcing::PickReport
;
pro Mgh_Roms_Plot_Tidal_Forcing::PickReport, pos, LUN=lun

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

pro Mgh_Roms_Plot_Tidal_Forcing__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, {Mgh_Roms_Plot_Tidal_Forcing, inherits MGH_Window, $
                 forcing_file: obj_new(), forcing_destroy: 0B, $
                 lonlat: 0B, varname: '', bar: obj_new(), $
                 plane: obj_new(), contour: obj_new()}

end
