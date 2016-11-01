;+
; CLASS NAME:
;   Mgh_Roms_Plot_Geke
;
; PURPOSE:
;   This class generates and displays a graph showing
;   geostrophic eddy kinetic energy (and other statistics) derived
;   from a series of Hslices of ROMS zeta data.
;
; CATEGORY:
;   Ocean models.
;   Object graphics
;
; CALLING SEQUENCE:
;   mgh_new, 'mgh_roms_plot_geke', history
;
; POSITIONAL PARAMETERS:
;   history
;     A reference to a ROMS history sequence object or a string arrary
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
; KEYWORD PARAMETERS:
;   DATA_RANGE (input, 2-element numeric)
;     Data range for the density surface.
;
;   MAP_STRUCTURE
;     If a map structure is supplied via this keyword, and if there
;     are longitude and latitude data in the file, then the longitudes
;     and latitudes are converted with MAP_PROJ_FORWARD before display.
;
;   RECORD_RANGE (input, 2-element numeric)
;     Range of record numbers (0-based indices in the time dimension) over
;     which statistics are to be calculated. The record range may also be
;     specified indirectly via the TIME_RANGE keyword.
;
;   TIME_RANGE (input, 2-element numeric)
;     Time range (days) over which statistics are to be calculated. The
;     record range may also be specified directly via the RECORD_RANgE
;     keyword.
;
;   XI_RANGE
;   ETA_RANGE
;     Use these keywords to display a subset of the domain. They are
;     interpreted relative to the rho grid.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2009-04:
;     Written.
;   Mark Hadfield, 2009-10:
;     Removed calls to widget_event(/NOWAIT).
;   Mark Hadfield, 2012-08:
;     - Removed HISTORY_DESTROY functionality.
;     - Added support for other parameters and implemented one of them:
;       mean eddy speed (MES)
;     - Removed code for a one-pass calculation
;   Mark Hadfield, 2015-11:
;     - Added TMPFILE property, allowing the caller to specify a file for
;       temporary storage of statistics.
;   Mark Hadfield, 2016-01:
;     - Added X_RANGE & Y_RANGE keywords.
;-
function mgh_roms_plot_geke::Init, history, $
     BYTE_RANGE=byte_range, $
     COLORBAR_PROPERTIES=colorbar_properties, $
     CONTOUR_PROPERTIES=contour_properties, $
     DATA_MULTIPLIER=data_multiplier, $
     DATA_RANGE=data_range, $
     GRAPH_PROPERTIES=graph_properties, $
     LAND_PROPERTIES=land_properties, $
     MAP_STRUCTURE=map_structure, $
     PALETTE_PROPERTIES=palette_properties, $
     PARAMETER=parameter, $
     RECORD_RANGE=record_range, $
     SHOW_CONTOUR=show_contour, $
     SHOW_PLANE=show_plane, $
     STYLE=style, TITLE=title, TIME_RANGE=time_range, $
     RECALC=recalc, TMPFILE=tmpfile, USE_ZETA=use_zeta, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
     X_RANGE=x_range, Y_RANGE=y_range, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process history argument.

   case size(history, /TYPE) of
      7: begin  ;;; String
         ohis = obj_new('MGHromsHistory', history)
         self.history_file = ohis
      end
      11: begin   ;;; Object reference
         ohis = history
         self.history_file = history
      end
      0: $
         message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'history'
      else: $
         message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'history'
   endcase

   ;; Other defaults


   if n_elements(parameter) eq 0 then parameter = 'GEKE'

   if n_elements(data_range) eq 0 then begin
      case strupcase(parameter) of
         'GEKE': data_range = [-3.5,-1.5]
         'GMES': data_range = [0,0.5]
      endcase
   endif

   if n_elements(style) eq 0 then style = 0

   if n_elements(show_contour) eq 0 then show_contour = 0B

   if n_elements(show_plane) eq 0 then show_plane = 1B

   if n_elements(map_structure) gt 0 then $
      self.map_structure = ptr_new(map_structure)

   use_map_structure = ptr_valid(self.map_structure)

   use_tmpfile = n_elements(tmpfile) gt 0

   if use_tmpfile then begin
      if ~ file_test(tmpfile) then recalc = 1B
   endif else begin
      recalc = 1B
   endelse

   ;; Get x & y positions for all RHO points (used for plotting walls
   ;; and land mask)

   self.lonlat = ohis->HasVar('lon_rho') && ohis->HasVar('lat_rho')

   x_rho = self.lonlat ? ohis->VarGet('lon_rho') : ohis->VarGet('x_rho')
   y_rho = self.lonlat ? ohis->VarGet('lat_rho') : ohis->VarGet('y_rho')

   mask_rho = ohis->VarGet('mask_rho')

   dim_rho = size(x_rho, /DIMENSIONS)

   ;; The default is to calculate geostrophic velocity from all zeta points.

   if n_elements(xi_range) eq 0 then xi_range = [0,dim_rho[0]-1]
   if n_elements(eta_range) eq 0 then eta_range = [0,dim_rho[1]-1]

   if xi_range[0] lt 0 then xi_range[0] += dim_rho[0]
   if xi_range[1] lt 0 then xi_range[1] += dim_rho[0]

   if eta_range[0] lt 0 then eta_range[0] += dim_rho[1]
   if eta_range[1] lt 0 then eta_range[1] += dim_rho[1]

   ;; Abbreviations for horizontal range

   xra0 = xi_range[0]
   xra1 = xi_range[1]
   xran = xra1-xra0+1
   era0 = eta_range[0]
   era1 = eta_range[1]
   eran = era1-era0+1

   ;; Get grid data for calculating geostrophic velocity

   pm = ohis->VarGet('pm', OFFSET=[xra0,era0], COUNT=[xran,eran])
   pn = ohis->VarGet('pn', OFFSET=[xra0,era0], COUNT=[xran,eran])

   f = ohis->VarGet('f', OFFSET=[xra0,era0], COUNT=[xran,eran])

   mask = mask_rho[xra0:xra1,era0:era1]

   ;; For plotting purposes, calculate locations for GEKE/GMES data

   x_val = mgh_stagger(x_rho[xra0:xra1,era0:era1], DELTA=[-1,-1])
   y_val = mgh_stagger(y_rho[xra0:xra1,era0:era1], DELTA=[-1,-1])

   ;; Convert all position data to map projection if appropriate

   if keyword_set(use_map_structure) then begin
      if ~ self.lonlat then $
         message, 'Cannot use map structure without lon, lat data'
      xy = map_proj_forward(x_val, y_val, MAP_STRUCTURE=*self.map_structure)
      x_val = reform(xy[0,*], size(x_val, /DIMENSIONS))
      y_val = reform(xy[1,*], size(y_val, /DIMENSIONS))
      mgh_undefine, xy
      xy = map_proj_forward(x_rho, y_rho, MAP_STRUCTURE=*self.map_structure)
      x_rho = reform(xy[0,*], size(x_rho, /DIMENSIONS))
      y_rho = reform(xy[1,*], size(y_rho, /DIMENSIONS))
      mgh_undefine, xy
   endif

   ;; Establish records to be processed

   case 1B of
      ohis->HasVar('ocean_time'): $
         time_var = 'ocean_time'
      ohis->HasVar('scrum_time'): $
         time_var = 'scrum_time'
      else: $
         message, 'Time variable not found'
   endcase

   time = ohis->VarGet(time_var, /AUTOSCALE)
   if n_elements(time_range) gt 0 then begin
      record_range = mgh_subset(time, time_range)
   endif
   if n_elements(record_range) eq 0 then begin
      n_time = ohis->DimInfo('ocean_time', /DIMSIZE)
      record_range = [0,n_time-1]
   endif
   if n_elements(time_range) eq 0 then $
      time_range = time[record_range]
   msg = ['Getting zeta data between records', $
      strtrim(record_range,2), $
      'times', mgh_format_float(time[record_range])]
   message, /INFORM, strjoin(temporary(msg), ' ')

   rra0 = record_range[0]
   rra1 = record_range[1]
   rran = rra1-rra0+1

   ;; Default graph aspect ratio, can be overridden via
   ;; GRAPH_PROPERTIES

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(x_val)
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(y_val)

   aspect = mgh_aspect(x_range, y_range, LONLAT=self.lonlat && (~ use_map_structure))
   aspect = (aspect > 0.4) < 1.5

   ;; Create base graph

   ograph = obj_new('MGHgrGraph2D', ASPECT=aspect, XMARGIN=[0.30,0.45], $
      NAME='ROMS horizontal slice geostrophic EKE', $
      _STRICT_EXTRA=graph_properties)

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Add mask around plot area

   ograph->NewMask

   ;; Add title

   if n_elements(title) eq 0 then begin
      fmt = '(%"ROMS %s %s-%s d")'
      title = string(FORMaT=fmt, parameter, mgh_format_float(time_range))
   endif

   ograph->NewTitle, title

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

   ;; Add walls. For each wall we extract the 2 rows/columns on each
   ;; side of the physical boundary from the x_rho & y_rho arrays into
   ;; variables xr & yr, call mgh_stagger to get xw & yw along the
   ;; wall, then trim the interior points off.

   walls = ohis->GetWalls()

   if walls[0] then begin
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

   if walls[1] then begin
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

   if walls[2] then begin
      ;; Eastern wall
      xw = mgh_stagger(x_rho[dim_rho[0]-2:dim_rho[0]-1,*], DELTA=[1,1])
      yw = mgh_stagger(y_rho[dim_rho[0]-2:dim_rho[0]-1,*], DELTA=[1,1])
      xw = xw[1:2,*]
      yw = yw[1:2,*]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', $
         DATAX=temporary(xw), DATAY=temporary(yw), DATAZ=temporary(zw), $
         STYLE=2, COLOR=[127,127,127]
   endif

   if walls[3] then begin
      ;; Northern wall
      xw = mgh_stagger(x_rho[*,dim_rho[1]-2:dim_rho[1]-1], DELTA=[1,1])
      yw = mgh_stagger(y_rho[*,dim_rho[1]-2:dim_rho[1]-1], DELTA=[1,1])
      xw = xw[*,1:2]
      yw = yw[*,1:2]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', $
         DATAX=temporary(xw), DATAY=temporary(yw), DATAZ=temporary(zw), $
         STYLE=2, COLOR=[127,127,127]
   endif

   ;; Draw land mask.

   if ohis->HasVar('mask_rho') then begin

      mask_missing = round(ohis->VarGet('mask_rho', OFFSET=[1,1], COUNT=dim_rho-2))

      ograph->NewAtom, 'MGHgrColorPlane', RESULT=oland, $
         STYLE=0, DEFAULT_COLOR=mgh_color('grey'), ZVALUE=-2*deltaz, $
         MISSING_POINTS=mask_missing, $
         DATAX=mgh_stagger(x_rho, DELTA=[-1,-1]), $
         DATAY=mgh_stagger(y_rho, DELTA=[-1,-1]), $
         NAME='Land mask'

      ;; Just playing really
      if oland->QueryProperty('STYLE') then $
         oland->SetPropertyAttribute, 'STYLE', SENSITIVE=0

   endif

   if keyword_set(recalc) then begin

      message, /INFORM, 'Calculating GEKE and GMES data'

      ;; Two calculations used to be available: one-pass and
      ;; two-pass. The latter is slower but potentially more accurate. It
      ;; ius not clear whether the accuracy of the one-pass method
      ;; was adequate in practice; however with the introduction of
      ;; calculations of MES, I have dropped the one-pass code anyway.

      fac = 1/double(rran)

      ;; The first pass accumulates mean geostrophic velocities.

      mean_u = fltarr(xran-1, eran-1)
      mean_v = fltarr(xran-1, eran-1)

      for r=rra0,rra1 do begin

         zeta = ohis->VarGet('zeta', OFFSET=[xra0,era0,r], COUNT=[xran,eran,1])

         ugs = mgh_roms_gsvel(zeta, pn, f, DIRECTION=0, MASK=mask)
         vgs = mgh_roms_gsvel(zeta, pm, f, DIRECTION=1, MASK=mask)

         ugs = mgh_stagger(temporary(ugs), DELTA=[-1,0])
         vgs = mgh_stagger(temporary(vgs), DELTA=[0,-1])

         mgh_undefine, zeta

         mean_u += fac*temporary(ugs)
         mean_v += fac*temporary(vgs)

      endfor

      ;; The second pass accumulates GEKE and GMES

      geke = dblarr(xran-1, eran-1)
      gmes = dblarr(xran-1, eran-1)

      for r=rra0,rra1 do begin

         zeta = ohis->VarGet('zeta', OFFSET=[xra0,era0,r], COUNT=[xran,eran,1])

         ugs = mgh_roms_gsvel(zeta, pn, f, DIRECTION=0, MASK=mask)
         vgs = mgh_roms_gsvel(zeta, pm, f, DIRECTION=1, MASK=mask)

         ugs = mgh_stagger(temporary(ugs), DELTA=[-1,0])
         vgs = mgh_stagger(temporary(vgs), DELTA=[0,-1])

         mgh_undefine, zeta

         geke += 0.5D0*fac*((ugs-mean_u)^2+(vgs-mean_v)^2)
         gmes += fac*sqrt((ugs-mean_u)^2+(vgs-mean_v)^2)

         mgh_undefine, ugs, vgs

      endfor

      mgh_undefine, mean_u, mean_v

      if use_tmpfile then begin
         message, /INFORM, 'Saving GEKE & GMES data to '+tmpfile
         save, FILE=tmpfile, geke, gmes
      endif

   endif else begin

      message, /INFORM, 'Restoring GEKE & GMES data from '+tmpfile
      restore, FILE=tmpfile

   endelse


   ;; Draw data

   case strupcase(parameter) of
      'GEKE': begin
         color_table = 'Matlab Jet'
         bar_title = 'Log!D10!N EKE (m!U2!Ns!U!Z(2212)2!N)'
         values = alog10(geke)
      end
      'GMES': begin
         color_table = 'MGH Speed'
         bar_title = 'Speed (m!U!Ns!U!Z(2212)1!N)'
         values = gmes
      end
   endcase

   ograph->NewPalette, color_table, RESULT=palette, $
      _STRICT_EXTRA=palette_properties

   ograph->NewColorBar, RESULT=obar, $
      DATA_RANGE=data_range, PALETTE=palette, $
      FONT=ograph->GetFont(POS=1), TITLE=bar_title, $
      SHOW_CONTOUR=show_contour, CONTOUR_PROPERTIES=contour_properties, $
      NAME='Colour bar', _STRICT_EXTRA=colorbar_properties
   self.bar = obar

   ograph->NewAtom, 'MGHgrDensityPlane', RESULT=oplane, $
      HIDE=(~ show_plane), STYLE=style, $
      DATAX=style eq 0 ? mgh_stagger(x_val, DELTA=[1,1]) : x_val, $
      DATAY=style eq 0 ? mgh_stagger(y_val, DELTA=[1,1]) : y_val, $
      DATA_VALUES=values, $
      ZVALUE=-5*deltaz , NAME='Data plane', $
      COLORSCALE=obar, /STORE_DATA
   self.plane = oplane

   ograph->NewAtom, 'IDLgrContour', RESULT=ocont, $
      HIDE=(~ show_contour), $
      GEOMZ=deltaz, /PLANAR, DATA=values, GEOMX=x_val, GEOMY=y_val, $
      _STRICT_EXTRA=contour_properties
   self.contour = ocont

   ;; Load graph into window

   ma = ['Magnify','Translate','Context']

   ok = self->MGH_Window::Init(GRAPHICS_TREE=ograph, MOUSE_ACTION=ma, _STRICT_EXTRA=extra)
   if ~ ok then $
      message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   ;; Finalise plot & return

   self->Finalize, 'Mgh_Roms_Plot_Geke'

   return, 1

end

; Mgh_Roms_Plot_Geke::Cleanup
;
pro Mgh_Roms_Plot_Geke::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Window::Cleanup

end

; Mgh_Roms_Plot_Geke::GetProperty
;
pro Mgh_Roms_Plot_Geke::GetProperty, $
     ALL=all, BAR=bar, BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
     HISTORY_FILE=history_file, LONLAT=lonlat, MAP_STRUCTURE=map_structure, $
     PALETTE=palette, PARAMETER=parameter, STYLE=style, _REF_EXTRA=extra

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

   self.plane->GetProperty, STYLE=style

   self.bar->GetProperty, $
        BYTE_RANGE=byte_range, DATA_RANGE=data_range, PALETTE=palette

   if arg_present(all) then $
        all = create_struct(all, 'bar', bar, 'byte_range', byte_range, $
                            'data_range', data_range, $
                            'history_file', history_file, $
                            'lonlat', lonlat, 'map_structure', map_structure, $
                            'palette', palette, 'parameter', parameter, 'style', style)

end

; Mgh_Roms_Plot_Geke::SetProperty
;
pro Mgh_Roms_Plot_Geke::SetProperty, $
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

; Mgh_Roms_Plot_Geke::About
;
pro Mgh_Roms_Plot_Geke::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::About, lun

   self->GetProperty, BAR=bar, HISTORY_FILE=history_file, $
        PALETTE=palette

   if obj_valid(history_file) then begin
      printf, lun, FORMAT='(%"%s: my history file sequence is %s")', $
           mgh_obj_string(self), mgh_obj_string(history_file)
      history_file->GetProperty, FILE_NAME=file_name
      printf, lun, self, ': its files are:', file_name
   endif

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

; Mgh_Roms_Plot_Geke::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Plot_Geke::BuildMenuBar

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

; Mgh_Roms_Plot_Geke::EventMenuBar
;
function Mgh_Roms_Plot_Geke::EventMenuBar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      'FILE.EXPORT.NETCDF': begin
         self.graphics_tree->GetProperty, NAME=name
         ext = '.nc'
         case strlen(name) of
            0: default_file = ''
            else: default_file = mgh_str_vanilla(name)+ext
         endcase
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

; Mgh_Roms_Plot_Geke::ExportData
;
pro Mgh_Roms_Plot_Geke::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::ExportData, values, labels

   self.plane->GetProperty, DATA_VALUES=data_values
   if n_elements(data_values) gt 0 then begin
      labels = [labels,'Data Values']
      values = [values,ptr_new(data_values, /NO_COPY)]
   endif

   self.plane->GetProperty, DATAX=datax, DATAY=datay
   if n_elements(datax)*n_elements(datay) gt 0 then begin
      labels = [labels,'Vertex X','Vertex Y']
      values = [values,ptr_new(datax, /NO_COPY),ptr_new(datay, /NO_COPY)]
   endif

   if self.lonlat && ptr_valid(self.map_structure) then begin
      labels = [labels,'Map structure']
      values = [values,ptr_new(*self.map_structure)]
   endif

end

; Mgh_Roms_Plot_Geke::ExportToNcFile
;
pro Mgh_Roms_Plot_Geke::ExportToNcFile, file

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->GetProperty, $
        LONLAT=lonlat, MAP_STRUCTURE=map_structure, $
        PARAMETER=parameter

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

   onc = obj_new('MGHncFile', file, /CREATE, /CLOBBER)

   onc->AttAdd, /GLOBAL, 'title', 'ROMS geostrophic EKE data'

   onc->AttAdd, /GLOBAL, 'history', $
                'Generated by routine MGH_ROMS_Plot_Geke::ExportToNcFile at '+ $
                mgh_dt_string(mgh_dt_now())

   onc->AttAdd, /GLOBAL, 'parameter', parameter

   onc->DimAdd, 'xi', dim[0]
   onc->DimAdd, 'eta', dim[1]

   x_name = lonlat ? 'lon' : 'x'
   y_name = lonlat ? 'lat' : 'y'
   v_name = mgh_str_vanilla(parameter)

   onc->VarAdd, x_name, ['xi','eta'], /DOUBLE
   onc->VarAdd, y_name, ['xi','eta'], /DOUBLE
   onc->VarAdd, v_name, ['xi','eta'], /FLOAT

   onc->VarPut, x_name, datax
   onc->VarPut, y_name, datay
   onc->VarPut, v_name, data_values

   obj_destroy, onc


end

; Mgh_Roms_Plot_Geke::PickReport
;
pro Mgh_Roms_Plot_Geke::PickReport, pos, LUN=lun

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

pro Mgh_Roms_Plot_Geke__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, {Mgh_Roms_Plot_Geke, inherits MGH_Window, $
                 history_file: obj_new(), $
                 lonlat: 0B, map_structure: ptr_new(), $
                 bar: obj_new(), plane: obj_new(), contour: obj_new()}

end
