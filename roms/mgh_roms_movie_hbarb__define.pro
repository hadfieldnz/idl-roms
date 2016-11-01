;+
; CLASS NAME:
;   Mgh_Roms_Movie_Hbarb
;
; PURPOSE:
;   This procedure generates and displays a series of animated graphs showing
;   a ROMS 2D or 3D velocity field on an x-y surface
;
; CALLING SEQUENCE:
;   mgh_new, 'mgh_roms_movie_hbarb', history
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
; Copyright (c) 2000-2016 NIWA:
;   http://www.niwa.co.nz/
; Licensed under the MIT open source license:
;   http://www.opensource.org/licenses/mit-license.php
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2000-01:
;     Written.
;   Mark Hadfield, 2001-06:
;     Updated for recent changes in MGH_GUI_Base event-handling conventions.
;   Mark Hadfield, 2001-08:
;     Updated for IDL 5.5.
;   Mark Hadfield, 2009-10:
;     Fixed a bug: dividing by RECORD_AVERAGE twice.
;   Mark Hadfield, 2011-03:
;     Updated SHOW_TIME functionality.
;   Mark Hadfield, 2011-07:
;     Removed history_destroy field: unnecessary with automatic garbage
;     collection.
;   Mark Hadfield, 2011-10:
;     Added SHOW_GSHHS & SHOW_LOGO functionality.
;   Mark Hadfield, 2012-10:
;     The variable argument is now a keyword.
;   Mark Hadfield, 2014-04:
;     Further enhancements to SHOW_TIME functionality.
;   Mark Hadfield, 2015-03:
;     Changed the vertical location of the barb object from -5*deltaz
;     to -1*deltaz so that it stays above the mask.
;   Mark Hadfield, 2016-03:
;     - Added SHOW_NZLAM_WIND functionality.
;-
function Mgh_Roms_Movie_Hbarb::Init, $
     history, $
     BARB_SCALE=barb_scale, $
     DEPTH=depth, LEVEL=level, SIGMA=sigma, $
     ETA_RANGE=eta_range, ETA_STRIDE=eta_stride, $
     KEY_BARB_MAGNITUDE=key_barb_magnitude, $
     KEY_UNITS=key_units, $
     MAP_STRUCTURE=map_structure, $
     RECORD_AVERAGE=record_average, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     RECORDS=records, $
     SHOW_LOGO=show_logo, $
     SHOW_BATHYMETRY=show_bathymetry, $
     SHOW_GSHHS=show_gshhs, $
     SHOW_TITLE=show_title, SHOW_TIME=show_time, $
     SHOW_NZLAM_WIND=show_nzlam_wind, $
     NZLAM_WIND_LOCATION=nzlam_wind_location, $
     NZLAM_WIND_N_SMOOTH=nzlam_wind_n_smooth, $
     TITLE=title, DT_FORMaT=dt_format, $
     HEAD_SIZE=head_size, SHOW_HEAD=show_head, SHOW_SYMBOL=show_symbol, $
     VARIaBLE=variable, $
     XI_RANGE=xi_range, XI_STRIDE=xi_stride, $
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
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_WRGNUMELEM', 'variable'

   if size(variable, /TNAME) ne 'STRING' then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'variable'

   self.variable = variable

   ;; Other defaults

   self.barb_scale = n_elements(barb_scale) gt 0 ? barb_scale : 0.1

   if n_elements(key_barb_magnitude) eq 0 then key_barb_magnitude = 0.2
   if n_elements(key_units) eq 0 then key_units = 'm/s'

   if n_elements(xi_stride) eq 0 then xi_stride = 1

   if n_elements(eta_stride) eq 0 then eta_stride = 1

   if n_elements(map_structure) gt 0 then self.map_structure = ptr_new(map_structure)

   use_map_structure = ptr_valid(self.map_structure)

   if n_elements(show_title) eq 0 then show_title = !true

   if keyword_set(show_title) then begin
      if n_elements(title) eq 0 then title = ''
      if n_elements(show_time) eq 0 then show_time = 1B
      if n_elements(dt_format) eq 0 then $
         dt_format = '(C(CYI4.4,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))'
   endif

   if n_elements(show_bathymetry) eq 0 then show_bathymetry = !true

   ;; Get x & y positions for all RHO points (used for plotting walls
   ;; and land mask)

   self.lonlat = ohis->HasVar('lon_rho') && ohis->HasVar('lat_rho')

   x_rho = self.lonlat ? ohis->VarGet('lon_rho') : ohis->VarGet('x_rho')
   y_rho = self.lonlat ? ohis->VarGet('lat_rho') : ohis->VarGet('y_rho')

   dim_rho = size(x_rho, /DIMENSIONS)

   ;; We are getting velocity data at U and V grid points and
   ;; interpolating them to RHO points so we have to do some of
   ;; dimension calculations normally handled by the HsliceData and
   ;; HsliceGrid methods ourselves. The default is to show all
   ;; interior RHO points.

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

   if grid_u.dims.time ne grid_v.dims.time then message, 'Inconsistent time dimensions'

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

   ;; Establish name of time variable & records to be plotted (if applicable)

   if n_elements(record_average) eq 0 then record_average = 1

   has_time = strlen(grid_u.dims.time) gt 0

   if has_time then begin
      n_time = ohis->DimInfo(grid_u.dims.time, /DIMSIZE)
      mgh_resolve_indices, n_time, record_range, record_stride, records
      n_records = n_elements(records)
      time_var = ohis->TimeVarName(grid_u.dims.time)
      if isa(time_var, /NULL) then message, 'Time variable not found'
   endif else begin
      n_records = 1
      ;; Check other record keywords for validity here?
      if record_average gt 1 then $
           message, 'Cannot average over records for time-independent variable'
   endelse

   n_frames = long(n_records)/long(record_average)

   ;; Get time units from the history file. Also retrieve the time time series.

   if has_time then begin
      if ohis->HasAtt(time_var, 'units') then begin
         time_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
      endif else begin
         time_units = {scale: 1, offset: 0}
      endelse
   endif

   time = ohis->VarGet(time_var, AUTOSCALE=0)*time_units.scale
   time = time[records]

   ;; Calculate angle between xi axis & the graph's x/lon axis from position data.
   ;; Ignore values stored in the history file as they may be incorrect when a
   ;; map structure is used.

   mgh_roms_xy_to_metric, x_var, y_var, $
        LONLAT=(self.lonlat && (~ use_map_structure)), ANGLE=angle

   ;; Get a wind time series for the animated NZLAM wind arrow

   if keyword_set(show_nzlam_wind) then begin
      if ~ self.lonlat then $
         message, 'Cannot draw a wind barb for a ROMS file without longitude-latitude data'
      if time_units.offset lt mgh_dt_julday('1900') then $
         message, 'Cannot extract NZLAM data for non-calendar times'
      ;; Determine location, which is specified in graph coordinates
      if n_elements(nzlam_wind_location) eq 0 then begin
         dim = size(x_var, /DIMENSIONS)
         nzlam_wind_location = [x_var[dim[0]/2,dim[1]/2],y_var[dim[0]/2,dim[1]/2]]
      endif
      ;; Determine the longitude & latitude at which NZLAM data are to be extracted, and the
      ;; angle between the map projection and the graph's coordinate system.
      if use_map_structure then begin
         ll = map_proj_inverse(nzlam_wind_location, MAP_STRUCTURE=*self.map_structure)
         ang = mgh_map_proj_angle(nzlam_wind_location[0], nzlam_wind_location[1], MAP_STRUCTURE=*self.map_structure)
      endif else begin
         ll = nzlam_wind_location
         ang = 0
      endelse
      ;; Extract NZLAM data. Spacing is 3-hourly, but with larger gaps, so interpolate
      ;; to a 3-hourly time vector with uniform spacing.
      nzlam_time_range = mgh_minmax(time+time_units.offset)+[-2,2]
      nzlam_time_3hr = timegen(START=nzlam_time_range[0], FINAL=nzlam_time_range[1], STEP_SIZE=3, UNITS='hours')
      nzlam = mgh_nzlam_sls_series(ll[0], ll[1], TIME_RANgE=nzlam_time_range)
      nzlam_wind_uv = interpol(complex(nzlam.u, nzlam.v), nzlam.time, nzlam_time_3hr)*exp(-!const.i*ang)
      mgh_undefine, nzlam
      ;; Optionally smooth the NZLAM (u,v) data
      if n_elements(nzlam_wind_n_smooth) gt 0 then $
         nzlam_wind_uv = smooth(nzlam_wind_uv, nzlam_wind_n_smooth, /EDGE_TRUNCATE)
      ;; Interpolate NZLAM wind data to model times
      nzlam_wind_uv = interpol(nzlam_wind_uv, nzlam_time_3hr, time+time_units.offset)
   endif

   ;; Default graph aspect ratio, can be overridden via GRAPH_PROPERTIES

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(mgh_stagger(x_var, DELTA=[1,1]))
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(mgh_stagger(y_var, DELTA=[1,1]))

   aspect = mgh_aspect(x_range, y_range, LONLAT=self.lonlat && (~ use_map_structure))
   aspect = (aspect > 0.4) < 1.5

   ;; Create base graph

   mgh_new, 'MGHgrGraph2D', ASPECT=aspect, $
            NAME='ROMS velocity barb animation', $
            _STRICT_EXTRA=graph_properties, RESULT=ograph

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize, $
        PLOT_RECT=prect, YMARGIN=ymargin

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Add mask around plot area

   ograph->NewMask

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
         ograph->NewAtom, 'IDLgrSurface', DATAX=temporary(xw), DATAY=temporary(yw), DATAZ=temporary(zw), $
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

      ocoast = mgh_gshhs_get_region(BOUNDARIES={lon: mgh_minmax(x_rho),lat: mgh_minmax(y_rho)}, $
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

   ;; Create a symbol for the base of each barb. This tends to slow
   ;; the animation down, so make the symbol reasonably simple.

   if keyword_set(show_symbol) then begin
      ograph->NewSymbol, CLASS='MGHgrSymbol', 0, RESULT=osym, $
           /FILL, NORM_SIZE=0.003, COLOR=mgh_color('blue'), $
           _STRICT_EXTRA=symbol_properties
   endif

   ;; Add a text object and a barb plot object to be animated

   if keyword_set(show_title) then $
      otitle = ograph->NewTitle(title)

   ograph->NewAtom, 'MGHgrBarb', $
        DATAX=x_var[(xi_stride-1)/2:*:xi_stride,(eta_stride-1)/2:*:eta_stride], $
        DATAY=y_var[(xi_stride-1)/2:*:xi_stride,(eta_stride-1)/2:*:eta_stride], $
        DATAZ=-1*deltaz, $
        SCALE=self.barb_scale, /NORM_SCALE, COLOR=mgh_color('red'), $
        SYMBOL=osym, HEAD_SIZE=head_size, SHOW_HEAD=show_head, NAME='Barbs', RESULT=obarb
   self.barb = obarb

   ;; ...NZLAM wind barb.

   if keyword_set(show_nzlam_wind) then begin
      ograph->NewSymbol, 0, RESULT=osym1, $
         XAXIS=0, YAXIS=0, FILL=1, NORM_SIZE=0.015, COLOR=!color.red
      ;; I'll sort out how to control where the barb appears on the page later.
      ograph->NewAtom, 'MGHgrBarb', RESULT=owind, $
         XAXIS=0, YAXIS=0, SCALE=0.01, COLOR=!color.red, /SHOW_HEAD, SYMBOL=osym1, $
         DATAX=prect[0]+0.2*prect[2], DATAY=prect[1]+0.8*prect[3], DATAU=0, DATAV=0
   endif

   ;; Create an animator window to display and manage the movie.

   ma = ['Magnify','Translate','Context']
   ok = self->MGH_Datamator::Init(CHANGEABLE=!false, GRAPHICS_TREE=ograph, MOUSE_ACTION=ma, _STRICT_EXTRA=extra)
   if ~ ok then message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Datamator'

   ;; Add a velocity-key barb

   kb = [prect[0]+0.5*prect[2]-0.1,prect[1]-ymargin[0]+0.1,10*deltaz]

   ograph->NewSymbol, CLASS='MGHgrSymbol', 0, RESULT=osym, $
           /FILL, NORM_SIZE=0.003, COLOR=mgh_color('blue'), $
           XAXIS=0, YAXIS=0

   ograph->NewAtom, 'MGHgrBarb', XAXIS=0, YAXIS=0, $
        DATAX=kb[0], DATAY=kb[1], DATAZ=kb[2], $
        DATAU=key_barb_magnitude, DATAV=0, $
        SCALE=self.barb_scale, /NORM_SCALE, COLOR=mgh_color('red'), $
        SHOW_HEAD=show_head, HEAD_SIzE=head_size, $
        SYMBOL=osym, NAME='Key barb', RESULT=okey
   self.key_barb = okey
   ograph->NewText, XAXIS=0, YAXIS=0, COLOR=mgh_color('red'), /ENABLE_FORMATTING, $
        STRINGS=mgh_format_float(key_barb_magnitude)+' '+key_units, $
        LOCATIONS=kb+[-0.02,0,0], ALIGNMENT=1, VERTICAL_ALIGNMENT=0.5

   ;; Step through records in the netCDF file, generating new frames &
   ;; plotting data

   oframe = objarr(3)

   ra = record_average

   for f=0,n_frames-1 do begin

      if self->Finished() then break

      if has_time then begin
         rec0 = ra*f
         frame_time = 0
         frame_wind_uv = 0
         frame_u = 0  &  frame_v = 0
         frame_wind_uv = 0
         for r=rec0,rec0+ra-1 do begin
            frame_time += time[r]/double(ra)
            if keyword_set(show_nzlam_wind) then frame_wind_uv += nzlam_wind_uv[r]/double(ra)
            uu = ohis->HsliceData(self.variable[0], GRID=grid_u, MASK_VALUE=0, $
                                  RECORD=records[r], DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
            uu[where(~ finite(uu), /NULL)] = 0
            frame_u += temporary(uu)/float(ra)
            vv = ohis->HsliceData(self.variable[1], GRID=grid_v, MASK_VALUE=0, $
                                  RECORD=records[r], DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
            vv[where(~ finite(vv), /NULL)] = 0
            frame_v += temporary(vv)/float(ra)
         endfor
      endif else begin
         frame_u = ohis->HsliceData(self.variable[0], GRID=grid_u, MASK_VALUE=0, $
                                    DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
         frame_u[where(~ finite(frame_u), /NULL)] = 0
         frame_v = ohis->HsliceData(self.variable[1], GRID=grid_v, MASK_VALUE=0, $
                                    DEPTHS=depth, LEVELS=level, SIGMAS=sigma)
         frame_v[where(~ finite(frame_v), /NULL)] = 0
      endelse

      ;; Interpolate to rho grid & rotate to the graph's coordinate system

      uv = complex(mgh_stagger(temporary(frame_u), DELTA=[-1,0]), $
                   mgh_stagger(temporary(frame_v), DELTA=[0,-1]))*exp(!const.i*angle)
      if xi_stride*eta_stride gt 1 then begin
         uv = uv[(xi_stride-1)/2:*:xi_stride,(eta_stride-1)/2:*:eta_stride]
      endif
      oframe[0] = obj_new('MGH_Command', OBJECT=obarb, $
                          'SetProperty', DATAU=real_part(uv), DATAV=imaginary(uv))
      mgh_undefine, uv

      ;; Update title with time or date, if applicable

      if keyword_set(show_title) && has_time && (show_time gt 0) then begin
         case show_time of
            1: ttt = string(FORMAT='(%"%0.3f days")', frame_time)
            2: ttt = string(FORMAT='(%"%s")', mgh_dt_string(frame_time+time_units.offset, FORMAT=dt_format))
            3: ttt = string(FORMAT='(%"%s (%0.3f days)")', mgh_dt_string(frame_time+time_units.offset, FORMAT=dt_format), frame_time)
         endcase
         if strlen(title) gt 0 then $
            ttt = string(FORMAT='(%"%s: %s")', title, ttt)
         oframe[1] = obj_new('MGH_Command', OBJECT=otitle, 'SetProperty', STRINGS=temporary(ttt))
      endif

      if keyword_set(show_nzlam_wind) then begin
         oframe[2] = obj_new('MGH_Command', OBJECT=owind, 'SetProperty', DATAU=real_part(frame_wind_uv), DATAV=imaginary(frame_wind_uv))
      endif

      ;; Add frame to animator & display

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; Mgh_Roms_Movie_Hbarb::Cleanup
;
pro Mgh_Roms_Movie_Hbarb::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Datamator::Cleanup

end

; Mgh_Roms_Movie_Hbarb::GetProperty
;
pro Mgh_Roms_Movie_Hbarb::GetProperty, $
     BARB_SCALE=barb_scale, HISTORY_FILE=history_file, $
     LONLAT=lonlat, VARIaBLE=variable, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::GetProperty, _STRICT_EXTRA=extra

   barb_scale = self.barb_scale

   history_file = self.history_file

   lonlat = self.lonlat

   variable = self.variable

end

; Mgh_Roms_Movie_Hbarb::SetProperty
;
pro Mgh_Roms_Movie_Hbarb::SetProperty, $
     BARB_SCALE=barb_scale, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=extra

   if n_elements(barb_scale) gt 0 then begin
      self.barb_scale = barb_scale
      self.barb->SetProperty, SCALE=self.barb_scale
      self.key_barb->SetProperty, SCALE=self.barb_scale
   endif

end

; Mgh_Roms_Movie_Hbarb::About
;
;   Print information about the window and its contents
;
pro Mgh_Roms_Movie_Hbarb::About, lun

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

   printf, lun, self, ': vector variables are '+strjoin(self.variable)

   printf, lun, self, ': the barb plot scale is '+strjoin(self.barb_scale)

end

; Mgh_Roms_Movie_Hbarb::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Movie_Hbarb::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

    self->MGH_Datamator::BuildMenuBar

    obar = mgh_widget_self(self.menu_bar)

    obar->NewItem, PARENT='File.Export Animation', 'NetCDF...'

    obar->NewItem, PARENT='Tools', SEPARATOR=[1,1,0], $
                   ['Set Barb Scale...','View U Data...','View V Data...']

end

; Mgh_Roms_Movie_Hbarb::EventMenubar
;
function Mgh_Roms_Movie_Hbarb::EventMenubar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      'FILE.EXPORT ANIMATION.NETCDF': begin
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
            self->ExportAnimationDataToNcFile, filename
         endif
         return, 0
      end

      'TOOLS.SET BARB SCALE': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Scale', CLIENT=self, $
            /FLOATING, GROUP_LEADER=self.base, IMMEDIATE=0, $
            N_ELEMENTS=1, PROPERTY_NAME='BARB_SCALE'
         return, 0
      end

      'TOOLS.VIEW U DATA': begin
         self->GetProperty, POSITION=position
         oframe = self.animation->GetFrame(POSITION=position)
         oframe[1]->GetProperty, KEYWORDS=keywords
         data_dims = size(keywords.datau, /DIMENSIONS)
         xvaredit, keywords.datau, GROUP=self.base, $
            X_SCROLL_SIZE=(data_dims[0] < 12), $
            Y_SCROLL_SIZE=(data_dims[1] < 30)
         return, 0
      end

      'TOOLS.VIEW V DATA': begin
         self->GetProperty, POSITION=position
         oframe = self.animation->GetFrame(POSITION=position)
         oframe[1]->GetProperty, KEYWORDS=keywords
         data_dims = size(keywords.datav, /DIMENSIONS)
         xvaredit, keywords.datav, GROUP=self.base, $
            X_SCROLL_SIZE=(data_dims[0] < 12), $
            Y_SCROLL_SIZE=(data_dims[1] < 30)
         return, 0
      end

      else: return, self->MGH_Datamator::EventMenubar(event)

   endcase

end

; Mgh_Roms_Movie_Hbarb::ExportData
;
pro Mgh_Roms_Movie_Hbarb::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::ExportData, values, labels

   self->GetProperty, $
        ANIMATION=animation, HISTORY_FILE=history_file, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[1]->GetProperty, KEYWORDS=keywords

   labels = [labels,'History Object','U,V Data']
   values = [values,ptr_new(history_file), $
             ptr_new(complex(keywords.datau,keywords.datav))]

end

; Mgh_Roms_Movie_Hbarb::ExportAnimationDataToNcFile
;
pro Mgh_Roms_Movie_Hbarb::ExportAnimationDataToNcFile, file

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Query the appropriate objects for info

   self->GetProperty, $
      LONLAT=lonlat, MAP_STRUCTURE=map_structure, VARIABLE=variable

   self.animation->GetProperty, $
      GRAPHICS_TREE=ograph, N_FRAMES=n_frames

   self.barb->GetProperty, $
      DATAX=datax, DATAY=datay

   self.animator->GetPlayBack, $
      RANGE=play_range, USE_RANGE=play_use_range

   ;; Sort out grid data

   dim = size(datax, /DIMENSIONS)

   if lonlat && size(map_structure, /TYPE) eq 8 then begin
      ll = map_proj_inverse(datax, datay, MAP_STRUCTURE=map_structure)
      datax = (reform(ll[0,*], dim)+360) mod 360
      datay = reform(ll[1,*], dim)
   endif

   ;; Sort out frames to be read

   if play_use_range then begin
      if n_elements(range) eq 0 then range = play_range[0:1]
      if n_elements(stride) eq 0 then stride = play_range[2]
   endif else begin
      if n_elements(range) eq 0 then range = [0,n_frames-1]
      if n_elements(stride) eq 0 then stride = 1
   endelse

   ;; Set up netCDF file

   onc = obj_new('MGHncFile', file, /CREATE, /CLOBBER)

   onc->AttAdd, /GLOBAL, 'title', 'ROMS hslice vector data'

   onc->AttAdd, /GLOBAL, 'history', $
      'Generated by routine Mgh_Roms_Movie_Hbarb::ExportAnimationDataToNcFile at '+ $
      mgh_dt_string(mgh_dt_now())

   onc->AttAdd, /GLOBAL, 'variables', strjoin(variable, ' ')

   onc->DimAdd, 'xi', dim[0]
   onc->DimAdd, 'eta', dim[1]
   onc->DimAdd, 'time'

   x_name = lonlat ? 'lon' : 'x'
   y_name = lonlat ? 'lat' : 'y'
   v_name = mgh_str_vanilla(variable)

   onc->VarAdd, x_name, ['xi','eta'], /DOUBLE
   onc->VarAdd, y_name, ['xi','eta'], /DOUBLE
   onc->VarAdd, v_name[0], ['xi','eta','time'], /FLOAT
   onc->VarAdd, v_name[1], ['xi','eta','time'], /FLOAT

   onc->VarPut, x_name, datax
   onc->VarPut, y_name, datay

   ;; Work through frames, retrieving and write data

   n_rec = 1 + (range[1]-range[0])/stride

   fmt ='(%"Writing %d frames of %d x %d vector data to netCDF file %s")'
   message, /INFORM, string(n_rec, dim, file, FORMAT=fmt)

   for r=0,n_rec-1 do begin

      pos = range[0]+stride*r

      oframe = self.animation->GetFrame(POSITION=pos)
      oframe[0]->GetProperty, KEYWORDS=keywords

      onc->VarPut, v_name[0], keywords.datau, OFFSET=[0,0,r]
      onc->VarPut, v_name[1], keywords.datav, OFFSET=[0,0,r]

   endfor

   obj_destroy, onc

   fmt ='(%"Finished saving netCDF file %s")'
   message, /INFORM, string(file, FORMAT=fmt)

end

pro Mgh_Roms_Movie_Hbarb__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {Mgh_Roms_Movie_Hbarb, inherits MGH_Datamator, $
         barb_scale: 0.0, variable: strarr(2), lonlat: !false, map_structure: ptr_new(), $
         history_file: obj_new(), barb: obj_new(), key_barb: obj_new()}

end
