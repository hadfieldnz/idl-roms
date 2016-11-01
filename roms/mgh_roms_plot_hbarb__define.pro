;+
; CLASS NAME:
;   Mgh_Roms_Plot_Hbarb
;
; PURPOSE:
;   This procedure plots a graph showing a ROMS 2D or 3D velocity field on an x-y surface
;
; CALLING SEQUENCE:
;   mgh_new, 'mgh_roms_Plot_Hbarb', history
;
; POSITIONAL PARAMETERS:
;   history (input)
;     A reference to a ROMS history sequence object or a string arrary
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
; KEYWORD PARAMETERS:
;   BARB_SCALE
;     Scale for the velocity barbs
;
;   DEPTH
;     Set this keyword to a scalar numeric value to specify the
;     depth of a z surface on which data are to be plotted. This
;     keyword should be specified only for variables having a depth
;     coordinate and it cannot be used together with LEVEL or SIGMA.
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
;   SHOW_TIME (input, integer)
;     Show the time or date in the title of each frame.
;
;   VARIABLE (input, optional)
;     A 2-element string array specifying the names of a pair of u & v
;     type variables. Default is ['u','v'] if available, otherwise
;     ['ubar','vbar'].
;
;   XI_RANGE
;   ETA_RANGE
;     Use these keywords to display a subset of the domain. They are
;     interpreted relative to the rho grid.
;
;###########################################################################
; Copyright (c) 2000-2012 NIWA:
;   http://www.niwa.co.nz/
; Licensed under the MIT open source license:
;   http://www.opensource.org/licenses/mit-license.php
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2016-01:
;     Written, based on MGH_ROMS_Movie_Hbarb and MGH_ROMS_Plot_Hslice.
;   Mark Hadfield, 2016-02:
;     Fixed bug: confusion between the symbols created for the velocity plot
;     and the key.
;-
function Mgh_Roms_Plot_Hbarb::Init, $
     history, $
     BARB_SCALE=barb_scale, $
     DEPTH=depth, LEVEL=level, SIGMA=sigma, $
     KEY_BARB_MAGNITUDE=key_barb_magnitude, $
     KEY_UNITS=key_units, $
     MAP_STRUCTURE=map_structure, $
     RECORD=record, $
     SHOW_BATHYMETRY=show_bathymetry, $
     SHOW_GSHHS=show_gshhs, $
     SHOW_TITLE=show_title, SHOW_TIME=show_time, $
     TITLE=title, DT_FORMaT=dt_format, $
     HEAD_SIZE=head_size, SHOW_HEAD=show_head, SHOW_SYMBOL=show_symbol, $
     VARIaBLE=variable, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
     XI_STRIDE=xi_stride, ETA_STRIDE=eta_stride, $
     X_RANGE=x_range, Y_RANGE=y_range, $
     BATH_PROPERTIES=bath_properties, $
     GRAPH_PROPERTIES=graph_properties, $
     SYMBOL_PROPERTIES=symbol_properties, $
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

   if n_elements(variable) eq 0 then begin
      if ohis->HasVar('u') && ohis->HasVar('v') then begin
         variable = ['u','v']
      endif else begin
         variable = ['ubar','vbar']
      endelse
   endif

   if n_elements(variable) ne 2 then $
        message, block='mgh_mblk_motley', name='mgh_m_wrgnumelem', 'variable'
   if size(variable, /TNAME) ne 'STRING' then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'variable'

   self.variable = variable

   ;; Other defaults

   self.barb_scale = n_elements(barb_scale) gt 0 ? barb_scale : 0.1

   if n_elements(key_barb_magnitude) eq 0 then key_barb_magnitude = 0.2
   if n_elements(key_units) eq 0 then key_units = 'm/s'

   if n_elements(xi_stride) eq 0 then xi_stride = 1
   if n_elements(eta_stride) eq 0 then eta_stride = 1

   if n_elements(map_structure) gt 0 then $
        self.map_structure = ptr_new(map_structure)

   use_map_structure = ptr_valid(self.map_structure)

   if n_elements(show_title) eq 0 then show_title = 1B

   if keyword_set(show_title) then begin
      if n_elements(title) eq 0 then title = ''
      if n_elements(show_time) eq 0 then show_time = 1B
      if n_elements(dt_format) eq 0 then $
         dt_format = '(C(CYI4.4,"-",CMOI2.2,"-",CDI2.2))'
   endif

   if n_elements(show_bathymetry) eq 0 then show_bathymetry = 1B

   ;; Get x & y positions for all RHO points (used for plotting walls
   ;; and land mask)

   self.lonlat = ohis->HasVar('lon_rho') && ohis->HasVar('lat_rho')

   x_rho = self.lonlat ? ohis->VarGet('lon_rho') : ohis->VarGet('x_rho')
   y_rho = self.lonlat ? ohis->VarGet('lat_rho') : ohis->VarGet('y_rho')

   dim_rho = size(x_rho, /DIMENSIONS)

   ;; We are getting velocity data at U and V grid points and
   ;; interpolating them to RHO points so we have to do some of
   ;; dimension calculations normally handled by the HsliceData and
   ;; HsliceGrid methods ourselves.

   ;;The default is to show all interior RHO points.

   if n_elements(xi_range) eq 0 then xi_range = [1,dim_rho[0]-2]
   if n_elements(eta_range) eq 0 then eta_range = [1,dim_rho[1]-2]

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

   ;; Get grid data required for horizontal slice retrievals

   grid_u = ohis->HsliceGrid(self.variable[0], ETA_RANGE=eta_range, $
                             LONLAT=self.lonlat, XI_RANGE=xi_range-[1,0])
   grid_v = ohis->HsliceGrid(self.variable[1], ETA_RANGE=eta_range-[1,0], $
                             LONLAT=self.lonlat, XI_RANGE=xi_range)

   if grid_u.dims.time ne grid_v.dims.time then $
        message, 'Inconsistent time dimensions'

   ;; Calculate x & y locations for slice

   x_var = mgh_stagger(grid_u.x, DELTA=[-1,0])
   y_var = mgh_stagger(grid_u.y, DELTA=[-1,0])

   ;; Convert all position data to map projection if appropriate

   if use_map_structure then begin
      if ~ self.lonlat then $
           message, 'Cannot use map structure without lon, lat data'
      xy = map_proj_forward(x_var, y_var, MAP_STRUCTURE=*self.map_structure)
      x_var = reform(xy[0,*], size(x_var, /DIMENSIONS))
      y_var = reform(xy[1,*], size(y_var, /DIMENSIONS))
      mgh_undefine, xy
      xy = map_proj_forward(x_rho, y_rho, MAP_STRUCTURE=*self.map_structure)
      x_rho = reform(xy[0,*], size(x_rho, /DIMENSIONS))
      y_rho = reform(xy[1,*], size(y_rho, /DIMENSIONS))
      mgh_undefine, xy
   endif

   ;; Calculate angle between xi axis & the graph's x/lon axis from position data.
   ;; Ignore values stored in the history file as they may be incorrect when a
   ;; map structure is used.

   mgh_roms_xy_to_metric, x_var, y_var, $
      LONLAT=(self.lonlat && (~ use_map_structure)), ANGLE=angle

   ;; Establish record to be plotted (if applicable). Establish name of
   ;; time variable--the logic here is a bit shaky, as conventions
   ;; for defining time dimensions and variables in ROMS files are
   ;; not entirely consistent.

   has_time = strlen(grid_u.dims.time) gt 0

   if has_time then begin
      if n_elements(record) eq 0 then record = 0
      case !true of
         ohis->HasVar(grid_u.dims.time): $
            time_var = grid_u.dims.time
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
         time_units = {scale: 1}
      endelse
   endif

   ;; Get velocity data

   if has_time then begin
      uu = ohis->HsliceData(self.variable[0], GRID=grid_u, MASK_VALUE=0, RECORD=record, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
      vv = ohis->HsliceData(self.variable[1], GRID=grid_v, MASK_VALUE=0, RECORD=record, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
   endif else begin
      uu = ohis->HsliceData(self.variable[0], GRID=grid_u, MASK_VALUE=0, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
      vv = ohis->HsliceData(self.variable[1], GRID=grid_v, MASK_VALUE=0, DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
   endelse

   uu[where(~ finite(uu), /NULL)] = 0
   vv[where(~ finite(vv), /NULL)] = 0

   ;; Interpolate to rho grid & rotate to the graph's coordinate system

   uv = complex(mgh_stagger(temporary(uu), DELTA=[-1,0]), $
                mgh_stagger(temporary(vv), DELTA=[0,-1]))*exp(!const.i*angle)

   ;; Subsample the velocity

   if xi_stride*eta_stride gt 1 then begin
      uv = uv[(xi_stride-1)/2:*:xi_stride,(eta_stride-1)/2:*:eta_stride]
   endif

   ;; Default graph aspect ratio, can be overridden via GRAPH_PROPERTIES

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(x_var)
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(y_var)

   aspect = mgh_aspect(x_range, y_range, LONLAT=self.lonlat && (~ use_map_structure))
   aspect = (aspect > 0.4) < 1.5

   ;; Create base graph

   ograph = obj_new('MGHgrGraph2D', ASPECT=aspect, $
                    NAME='ROMS barb plot', $
                    _STRICT_EXTRA=graph_properties)

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize, $
        PLOT_RECT=prect, YMARGIN=ymargin

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

   ;; Draw axes. (Not sure about those limits! Should calculate them
   ;; from x_rho & y_rho.)

   if self.lonlat && (~ use_map_structure) then begin
      xap = {tickformat: 'mgh_tf_longitude', tickfrmtdata: {format:'(F10.1)'}}
      yap = {tickformat: 'mgh_tf_latitude', tickfrmtdata: {format:'(F10.1)'} }
   endif else begin
      xap = {title: 'X (km)', tickformat: 'mgh_tf_linear', tickfrmtdata: {scale:1.E-3, format:'(F10.1)'}}
      yap = {title: 'Y (km)', tickformat: 'mgh_tf_linear', tickfrmtdata: {scale:1.E-3, format:'(F10.1)'}}
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

   ;; Draw bathymetry contours.

   if keyword_set(show_bathymetry) && ohis->HasVar('h') then begin

      hr = ohis->VarGet('h', OFFSET=[1,1], COUNT=dim_rho-2)

      ograph->NewAtom, 'IDLgrContour', RESULT=obath, $
           DATA=hr, GEOMZ=-8*deltaz, /PLANAR, $
           GEOMX=x_rho[1:dim_rho[0]-2,1:dim_rho[1]-2], $
           GEOMY=y_rho[1:dim_rho[0]-2,1:dim_rho[1]-2], $
           C_COLOR=mgh_color(['blue','blue']), $
           NAME='Bathymetry', _STRICT_EXTRA=bath_properties

   endif

   if keyword_set(show_gshhs) then begin

      if ~ self.lonlat then $
           message, 'Cannot draw GSHHS coastline without lon, lat data'

      ocoast = mgh_gshhs_get_region(BOUNDARIES={lon: mgh_minmax(x_rho), $
                                                lat: mgh_minmax(y_rho)}, $
                                    FILL=0, RESOLUTION=2, CLIP=0, $
                                    _STRICT_EXTRA=gshhs_properties)

      while ocoast->Count() gt 0 do begin

         opoly = ocoast->Get()
         ocoast->Remove, opoly

         if use_map_structure then begin
            opoly->GetProperty, DATA=data
            data[0:1,*] = map_proj_forward(data[0:1,*], MAP_STRUCTURE=map_structure)
            opoly->SetProperty, DATA=data
         endif

         opoly->SetProperty, COLOR=!color.green, THICK=2

         ograph->AddAtom, opoly

       endwhile

       obj_destroy, ocoast

   endif

   ;; Add a velocity-key barb

   kb = [prect[0]+0.5*prect[2]-0.1,prect[1]-ymargin[0]+0.1,10*deltaz]

   if keyword_set(show_symbol) then begin
      ograph->NewSymbol, CLASS='MGHgrSymbol', 0, RESULT=osym, $
         /FILL, NORM_SIZE=0.003, COLOR=mgh_color('blue'), $
         XAXIS=0, YAXIS=0
   endif

   ograph->NewAtom, 'MGHgrBarb', XAXIS=0, YAXIS=0, $
      DATAX=kb[0], DATAY=kb[1], DATAZ=kb[2], $
      DATAU=key_barb_magnitude, DATAV=0, $
      SCALE=self.barb_scale, /NORM_SCALE, COLOR=mgh_color('red'), $
      SHOW_HEAD=show_head, HEAD_SIzE=head_size, $
      SYMBOL=temporary(osym), NAME='Key barb', RESULT=okey
   self.key_barb = okey
   ograph->NewText, XAXIS=0, YAXIS=0, COLOR=mgh_color('red'), /ENABLE_FORMATTING, $
      STRINGS=mgh_format_float(key_barb_magnitude)+' '+key_units, $
      LOCATIONS=kb+[-0.02,0,0], ALIGNMENT=1, VERTICAL_ALIGNMENT=0.5

   ;; Create a symbol for the base of each barb.

   if keyword_set(show_symbol) then begin
      ograph->NewSymbol, CLASS='MGHgrSymbol', 0, RESULT=osym, $
           /FILL, NORM_SIZE=0.003, COLOR=mgh_color('blue'), $
           _STRICT_EXTRA=symbol_properties
   endif

   ;; Show velocity in a barb object

   ograph->NewAtom, 'MGHgrBarb', $
        DATAX=x_var[(xi_stride-1)/2:*:xi_stride,(eta_stride-1)/2:*:eta_stride], $
        DATAY=y_var[(xi_stride-1)/2:*:xi_stride,(eta_stride-1)/2:*:eta_stride], $
        DATAU=real_part(uv), DATAV=imaginary(uv), $
        DATAZ=-1*deltaz, $
        SCALE=self.barb_scale, /NORM_SCALE, COLOR=mgh_color('red'), $
        SYMBOL=temporary(osym), HEAD_SIZE=head_size, SHOW_HEAD=show_head, NAME='Barbs', RESULT=obarb
   self.barb = obarb

   ;; Load graph into window

   ma = ['Magnify','Translate','Context']
   ok = self->MGH_Window::Init(GRAPHICS_TREE=ograph, MOUSE_ACTION=ma, _STRICT_EXTRA=extra)
   if ~ ok then $
      message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   ;; Finalise plot & return

   self->Finalize, 'Mgh_Roms_Plot_Hbarb'

   return, 1


end

; Mgh_Roms_Plot_Hbarb::Cleanup
;
pro Mgh_Roms_Plot_Hbarb::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Window::Cleanup

end

; Mgh_Roms_Plot_Hbarb::GetProperty
;
pro Mgh_Roms_Plot_Hbarb::GetProperty, $
     BARB_SCALE=barb_scale, HISTORY_FILE=history_file, $
     LONLAT=lonlat, VARIaBLE=variable, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::GetProperty, _STRICT_EXTRA=extra

   barb_scale = self.barb_scale

   history_file = self.history_file

   lonlat = self.lonlat

   variable = self.variable

end

; Mgh_Roms_Plot_Hbarb::SetProperty
;
pro Mgh_Roms_Plot_Hbarb::SetProperty, $
     BARB_SCALE=barb_scale, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::SetProperty, _STRICT_EXTRA=extra

   if n_elements(barb_scale) gt 0 then begin
      self.barb_scale = barb_scale
      self.barb->SetProperty, SCALE=self.barb_scale
      self.key_barb->SetProperty, SCALE=self.barb_scale
   endif

end

; Mgh_Roms_Plot_Hbarb::About
;
;   Print information about the window and its contents
;
pro Mgh_Roms_Plot_Hbarb::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::About, lun

   if obj_valid(self.history_file) then begin
      printf, lun, self, ': the history file sequence is ', self.history_file
      self.history_file->GetProperty, FILE_NAME=file_name
      printf, lun, self, ': the files are:', file_name
   endif

   printf, lun, self, ': vector variables are '+strjoin(self.variable)

   printf, lun, self, ': the barb plot scale is '+strjoin(self.barb_scale)

end

; Mgh_Roms_Plot_Hbarb::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Plot_Hbarb::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

    self->MGH_Window::BuildMenuBar

    obar = mgh_widget_self(self.menu_bar)

    obar->NewItem, PARENT='Tools', SEPARATOR=[1,1,0], $
                   ['Set Barb Scale...','View U Data...','View V Data...']

end

; Mgh_Roms_Plot_Hbarb::EventMenubar
;
function Mgh_Roms_Plot_Hbarb::EventMenubar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      'TOOLS.SET BARB SCALE': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Scale', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, IMMEDIATE=0, $
                  N_ELEMENTS=1, PROPERTY_NAME='BARB_SCALE'
         return, 0
      end

      'TOOLS.VIEW U DATA': begin
         self.barb->GetProperty, DATAU=datau
         data_dims = size(datau, /DIMENSIONS)
         xvaredit, datau, GROUP=self.base, $
                   X_SCROLL_SIZE=(data_dims[0] < 12), $
                   Y_SCROLL_SIZE=(data_dims[1] < 30)
         return, 0
      end

      'TOOLS.VIEW V DATA': begin
         self.barb->GetProperty, DATAV=datav
         data_dims = size(datav, /DIMENSIONS)
         xvaredit, datav, GROUP=self.base, $
                   X_SCROLL_SIZE=(data_dims[0] < 12), $
                   Y_SCROLL_SIZE=(data_dims[1] < 30)
         return, 0
      end

      else: return, self->MGH_Window::EventMenubar(event)

   endcase

end

; Mgh_Roms_Plot_Hbarb::ExportData
;
pro Mgh_Roms_Plot_Hbarb::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::ExportData, values, labels

   self.barb->GetProperty, $
      DATAU=datau, DATAV=datav, DATAX=datax, DATAY=datay

   labels = [labels,'U, V Data']
   values = [values,ptr_new(complex(datau, datav))]

   if n_elements(datax)*n_elements(datay) gt 0 then begin
      if self.lonlat then begin
         if ptr_valid(self.map_structure) then begin
            dim = size(datax, /DIMENSIONS)
            lonlat = map_proj_inverse(datax, datay, MAP_STRUCTURE=*self.map_structure)
            lon = reform(lonlat[0,*], dim)
            lat = reform(lonlat[1,*], dim)
            mgh_undefine, lonlat
            labels = [labels,'Lon','Lat','X','Y','Map structure']
            values = [values,ptr_new(lon, /NO_COPY),ptr_new(lat, /NO_COPY), $
               ptr_new(datax, /NO_COPY),ptr_new(datay, /NO_COPY), $
               ptr_new(*self.map_structure)]
         endif else begin
            labels = [labels,'Lon','Lat']
            values = [values,ptr_new(datax, /NO_COPY),ptr_new(datay, /NO_COPY)]
         endelse
      endif else begin
         labels = [labels,'Vertex X','Vertex Y']
         values = [values,ptr_new(datax, /NO_COPY),ptr_new(datay, /NO_COPY)]
      endelse
   endif


end

pro Mgh_Roms_Plot_Hbarb__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
      {Mgh_Roms_Plot_Hbarb, inherits MGH_Window, $
       barb_scale: 0.0, variable: strarr(2), lonlat: !false, $
       history_file: obj_new(), map_structure: ptr_new(), $
       barb: obj_new(), key_barb: obj_new()}

end
