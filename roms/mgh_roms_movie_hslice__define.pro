;+
; CLASS NAME:
;  MGH_ROMS_Movie_Hslice
;
; PURPOSE:
;   This class generates and displays an animated sequence of graphs
;   showing an Hslice through a ROMS 2D or 3D output field on an x-y surface.
;
; CALLING SEQUENCE:
;   mgh_new, 'mgh_roms_movie_hslice', history, variable
;
; POSITIONAL PARAMETERS:
;   history (input)
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
;     sigma values of a copnstant-sigma surface on which data are to be
;     plotted. This keyword should be specified only for variables having
;     a depth coordinate and it cannot be used together with DEPTH or LEVEL.
;
;   SHOW_TIME (input, integer)
;     Show the time or date in the title of each frame.
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
;   Mark Hadfield, Apr 1999:
;     Written as ROMS_HSLICE_MOVIE.
;   Mark Hadfield, Aug 1999:
;     Miscellaneous improvements including sequences of files &
;     LEVEL keyword.
;   Mark Hadfield, Jan 2000:
;     Converted from a procedure (ROMS_HSLICE_MOVIE) to an
;     object. The front-end procedure (now MGH_ROMS_MOVIE_HSLICE)
;     can be used to create the object.
;   Mark Hadfield, Aug 2000:
;     Class name is now Mgh_Roms_Movie_Hslice and it can be created
;     easily using MGH_NEW. The history file sequence now *must* be
;     in the form of an object.
;   Mark Hadfield, Sep 2000:
;     Added (reinstated?) code to accept a list or file names or an
;     MGHromsHistory object as the first positional argument. In
;     either case an MGHromsHistory object is created (if necessary)
;     and stored with the present object.
;   Mark Hadfield, 2000-11:
;     Implemented ETA_RANGE & XI_RANGE keywords to allow selection
;     of subsets. Increased flexibility of the code handling the
;     variable's dimensions, so that variables without a time
;     dimension can be plotted.
;   Mark Hadfield, 2000-12:
;     Rewritten to take advantage of MGHromsHistory's new HsliceGrid
;     & HsliceData methods.
;   Mark Hadfield, 2004-04:
;     Added support for overplotted contours via SHOW_CONTOUR &
;     CONTOUR_PROPERTIES keywords.
;   Mark Hadfield, 2008-11:
;     Added support for map projections via the MAP_STRUCTURE
;     keyword.
;   Mark Hadfield, 2009-03:
;     Added support for sigma slices via the SIGMA keyword.
;   Mark Hadfield, 2010-11:
;     Wall plotting code cleaned up.
;   Mark Hadfield, 2010-12:
;     Default STYLE is now 0 for all variables.
;   Mark Hadfield, 2011-03:
;     Updated SHOW_TIME functionality.
;   Mark Hadfield, 2011-05:
;     - Added support for sediment bed-layer variables (LAYER keyword)
;       and dynamic bathymetry (USE_BATH keyword).
;     - The defaults for USE_BATH and USE_ZETA are now provided by
;       MGHromsHistory::HsliceData.
;   Mark Hadfield, 2012-02:
;     - Fixed minor bug in code that exports slice data.
;   Mark Hadfield, 2012-10:
;     - Added LOGARITHMIC property.
;   Mark Hadfield, 2013-02:
;     - Added SHOW_COLORBAR property.
;   Mark Hadfield, 2013-09:
;     - Removed extraneous references to "parameter" in the GetProperty method.
;   Mark Hadfield, 2014-03:
;     Further enhancements to SHOW_TIME functionality.
;   Mark Hadfield, 2015-05:
;     Even more enhancements to SHOW_TIME functionality.
;   Mark Hadfield, 2015-06:
;     The format of the time string (with SHOW_TIME equal to 1 or 2) can now be
;     specified.
;   Mark Hadfield, 2015-12:
;     - The default XI_RANGE and ETA_RANgE are now set to retrieve all data.
;     - All cells of the land mask are now plotted.
;   Mark Hadfield, 2016-01:
;     - The default colour table is now "Matlab Jet"
;   Mark Hadfield, 2016-03:
;     - Added SHOW_NZLAM_WIND functionality.
;   Mark Hadfield, 2016-05:
;     - Added "File.Export Animation.NetCDF..." functionality.
;   Mark Hadfield, 2016-07:
;     - Added a keyword to control the position (in normalised coordinates) of the
;       NZLAM wind barb. The keyword is called NZLAM_BARB_ORIGIN and the default is
;       the position )converted to normalised coordinates) at which the data are extracted.
;-
function mgh_roms_movie_hslice::Init, $
     history, variable, $
     BYTE_RANGE=byte_range, $
     DATA_MULTIPLIER=data_multiplier, $
     DATA_RANGE=data_range, LOGARITHMIC=logarithmic, $
     DATA_TRANSFORMATION=data_transformation, $
     DEPTH=depth, LAYER=layer, LEVEL=level, SIGMA=sigma, $
     MAP_STRUCTURE=map_structure, $
     MASK_VALUE=mask_value, $
     RECORD_AVERAGE=record_average, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     RECORDS=records, $
     STYLE=style, $
     SHOW_COLORBAR=show_colorbar, SHOW_CONTOUR=show_contour, $
     SHOW_TITLE=show_title, SHOW_TIME=show_time, $
     SHOW_NZLAM_WIND=show_nzlam_wind, $
     NZLAM_BARB_ORIGIN=nzlam_barb_origin, $
     NZLAM_WIND_LOCATION=nzlam_wind_location, $
     NZLAM_WIND_N_SMOOTH=nzlam_wind_n_smooth, $
     TITLE=title, DT_FORMaT=dt_format, $
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

   if n_elements(show_colorbar) eq 0 then show_colorbar = !true

   if n_elements(show_contour) eq 0 then show_contour = !false

   if n_elements(map_structure) gt 0 then self.map_structure = ptr_new(map_structure)

   use_map_structure = ptr_valid(self.map_structure)

   if n_elements(show_title) eq 0 then show_title = !true

   if keyword_set(show_title) then begin
      if n_elements(title) eq 0 then title = ''
      if n_elements(show_time) eq 0 then show_time = 1
      if n_elements(dt_format) eq 0 then $
         dt_format = '(C(CYI4.4,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))'
   endif

   ;; Set ETA_RANGE and XI_RANGE relative to the rho grid.
   ;; The default is to show all interior points.

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
      if ~ self.lonlat then $
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

   ;; Get time units from the history file. Also retrieve the time series.

   if has_time then begin
      if ohis->HasAtt(time_var, 'units') then begin
         time_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
      endif else begin
         time_units = {scale: 1, offset: 0}
      endelse
   endif

   time = ohis->VarGet(time_var, AUTOSCALE=0)*time_units.scale
   time = time[records]

   ;; Get a wind time series for the animated NZLAM wind arrow

   if keyword_set(show_nzlam_wind) then begin
      if ~ self.lonlat then $
         message, 'Cannot draw a wind barb for a ROMS file without longitude-latitude data'
      if time_units.offset lt mgh_dt_julday('1900') then $
         message, 'Cannot extract NZLAM data for non-calendar times'
      ;; Determine location, which is specified in graph coordinates. The default is the
      ;; position of the grid point in the centre of the domain.
      if n_elements(nzlam_wind_location) eq 0 then begin
         dim = size(grid.x, /DIMENSIONS)
         nzlam_wind_location = [grid.x[dim[0]/2,dim[1]/2],grid.y[dim[0]/2,dim[1]/2]]
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

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(x_var)
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(y_var)

   aspect = mgh_aspect(x_range, y_range, LONLAT=self.lonlat && (~ use_map_structure))
   aspect = (aspect > 0.4) < 1.5

   ;; Create base graph

   xmargin = show_colorbar ? [0.375,0.4] : [0.375,0.15]

   ograph = obj_new('MGHgrGraph2D', ASPECT=aspect, XMARGIN=xmargin, $
                    NAME='ROMS horizontal slice animation', $
                    _STRICT_EXTRA=graph_properties)

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize, PLOT_RECT=prect

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Add mask around plot area

   ograph->NewMask

   ;; Draw axes

   if self.lonlat && (~ use_map_structure) then begin
      xap = {tickformat: 'mgh_tf_longitude', tickfrmtdata: {format:'(F10.1)'}}
      yap = {tickformat: 'mgh_tf_latitude', tickfrmtdata: {format:'(F10.1)'} }
   endif else begin
      xap = {title: 'X (km)', tickformat: 'mgh_tf_linear', tickfrmtdata: {scale:1.E-3, format:'(F10.1)'}}
      yap = {title: 'Y (km)', tickformat: 'mgh_tf_linear', tickfrmtdata: {scale:1.E-3, format:'(F10.1)'}}
   endelse

   ograph->NewAxis, 0, $
        RANGE=x_range, /EXACT, $
        _STRICT_EXTRA=mgh_struct_merge(xap, xaxis_properties)
   ograph->NewAxis, 1, $
        RANGE=y_range, /EXACT, $
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

      mask_missing = round(ohis->VarGet('mask_rho'))

      ograph->NewAtom, 'MGHgrColorPlane', RESULT=oland, $
         NAME='Land mask', STYLE=0, DEFAULT_COLOR=mgh_color('grey'), $
         ZVALUE=-2*deltaz, MISSING_POINTS=mask_missing, $
         DATAX=mgh_stagger(x_rho, DELTA=[1,1]), $
         DATAY=mgh_stagger(y_rho, DELTA=[1,1]), $
         _STRICT_eXTRA=land_properties

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

   ;; Add various graphics objects that will (or might) be animated...

   ;; ...title

   if keyword_set(show_title) then otitle = ograph->NewTitle(title)

   ;; ...density plane showing data values.

   ograph->NewAtom, 'MGHgrDensityPlane', RESULT=oplane, STYLE=style, ZVALUE=-5*deltaz , $
        DATAX=x_var, DATAY=y_var, DATA_VALUES=mgh_reproduce(0,grid.x), $
        COLORSCALE=self.bar, /STORE_DATA
   self.plane = oplane

   ;; ...contour showing data values.

   if show_contour then begin
      ograph->NewAtom, 'IDLgrContour', RESULT=ocont, /PLANAR, GEOMZ=-2*deltaz , $
           GEOMX=grid.x, GEOMY=grid.y, DATA=mgh_reproduce(0, grid.x), $
           _STRICT_EXTRA=contour_properties
      self.contour = ocont
   endif

   ;; ...NZLAM wind barb.

   if keyword_set(show_nzlam_wind) then begin
      ograph->NewSymbol, 0, RESULT=osym1, $
         XAXIS=0, YAXIS=0, FILL=1, NORM_SIZE=0.015, COLOR=!color.red
      if n_elements(nzlam_barb_origin) eq 0 then $
         nzlam_barb_origin = ograph->NormPosition(nzlam_wind_location)
      ograph->NewAtom, 'MGHgrBarb', RESULT=owind, $
         XAXIS=0, YAXIS=0, SCALE=0.01, COLOR=!color.red, /SHOW_HEAD, SYMBOL=osym1, $
         DATAX=nzlam_barb_origin[0], DATAY=nzlam_barb_origin[1], DATAU=0, DATAV=0
   endif

   ;; Create an animator window to display and manage the movie.

   ma = ['Magnify','Translate','Context']
   ok = self->MGH_Datamator::Init(CHANGEABLE=!false, GRAPHICS_TREE=ograph, MOUSE_ACTION=ma, _STRICT_EXTRA=extra)
   if ~ ok then message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Datamator'

   ;; Step through the netCDF file, generating new frames & plotting data

   oframe = objarr(4)

   ra = record_average
   dm = data_multiplier

   for f=0,n_frames-1 do begin

      if self->Finished() then break

      rec0 = ra*f

      frame_time = 0
      frame_slice = 0

      frame_wind_uv = 0

      for r=rec0,rec0+ra-1 do begin

         if has_time then begin
            frame_time += time[r]/double(ra)
            if keyword_set(show_nzlam_wind) then frame_wind_uv += nzlam_wind_uv[r]/double(ra)
            frame_slice += ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                                      RECORD=records[r], $
                                      DEPTHS=depth, LAYERS=layer, LEVELS=level, $
                                      SIGMAS=sigma, USE_BATH=use_bath, USE_ZETA=use_zeta)
         endif else begin
            frame_slice += ohis->HsliceData(variable, GRID=grid, MASK_VALUE=mask_value, $
                                      DEPTHS=depth, LAYERS=layer, LEVELS=level, $
                                      SIGMAS=sigma, USE_BATH=use_bath, USE_ZETA=use_zeta)
         endelse

      endfor

      frame_slice = dm*frame_slice/float(ra)

      if n_elements(data_transformation) gt 0 then $
           frame_slice = call_function(data_transformation, frame_slice)

      ;; Note that if the LOGARITHMIC keyword is set, the logarithmic transformation
      ;; is handled within the methods of the MGHgrDensityPlane class.

      oframe[0] = obj_new('MGH_Command', OBJECT=self.plane, 'SetProperty', DATA_VALUES=frame_slice)
      if show_contour then $
           oframe[1] = obj_new('MGH_Command', OBJECT=self.contour, 'SetProperty', DATA=frame_slice)
     mgh_undefine, frame_slice

      ;; Update title with time or date, if applicable

      if keyword_set(show_title) && has_time && show_time gt 0 then begin
         case show_time of
            1: ttt = string(FORMAT='(%"%0.3f days")', frame_time)
            2: ttt = string(FORMAT='(%"%s")', mgh_dt_string(frame_time+time_units.offset, FORMAT=dt_format))
            3: ttt = string(FORMAT='(%"%s (%0.3f days)")', mgh_dt_string(frame_time+time_units.offset, FORMAT=dt_format), frame_time)
         endcase
         if strlen(title) gt 0 then $
            ttt = string(FORMAT='(%"%s: %s")', title, ttt)
         oframe[2] = obj_new('MGH_Command', OBJECT=otitle, 'SetProperty', STRINGS=temporary(ttt))
      endif

      if keyword_set(show_nzlam_wind) then begin
         oframe[3] = obj_new('MGH_Command', OBJECT=owind, 'SetProperty', DATAU=real_part(frame_wind_uv), DATAV=imaginary(frame_wind_uv))
      endif

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; Mgh_Roms_Movie_Hslice::Cleanup
;
pro Mgh_Roms_Movie_Hslice::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Datamator::Cleanup

end

; Mgh_Roms_Movie_Hslice::GetProperty
;
pro Mgh_Roms_Movie_Hslice::GetProperty, $
     ALL=all, BAR=bar, BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
     HISTORY_FILE=history_file, LONLAT=lonlat, MAP_STRUCTURE=map_structure, $
     PALETTE=palette, STYLE=style, VARIABLE=variable, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::GetProperty, ALL=all, _STRICT_EXTRA=extra

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

; Mgh_Roms_Movie_Hslice::SetProperty
;
pro Mgh_Roms_Movie_Hslice::SetProperty, $
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

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=extra

end

; Mgh_Roms_Movie_Hslice::About
;
pro Mgh_Roms_Movie_Hslice::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::About, lun

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

; Mgh_Roms_Movie_Hslice::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Movie_Hslice::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

    self->MGH_Datamator::BuildMenuBar

    obar = mgh_widget_self(self.menu_bar)

    obar->NewItem, PARENT='File.Export Animation', 'NetCDF...'

   obar->NewItem, PARENT='Tools', SEPARATOR=[1,0,0,0,1], MENU=[1,0,0,0,0], $
        ['Data Range','Edit Palette...','Set Style...', $
         'View Colour Scale...','View Data Values...']

   obar->NewItem, PARENT='Tools.Data Range', ['Set...','Fit this Frame']

end


; Mgh_Roms_Movie_Hslice::EventMenuBar
;
function Mgh_Roms_Movie_Hslice::EventMenuBar, event

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

      'TOOLS.DATA RANGE.SET': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Range', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, IMMEDIATE=0, $
                  N_ELEMENTS=2, PROPERTY_NAME='DATA_RANGE'
         return, 0
      end

      'TOOLS.DATA RANGE.FIT THIS FRAME': begin
         self->GetProperty, POSITION=position
         oframe = self.animation->GetFrame(POSITION=position)
         oframe[0]->GetProperty, KEYWORDS=keywords
         data_range = mgh_minmax(keywords.data_values, /NAN)
         if data_range[0] eq data_range[1] then data_range += [-1,1]
         self->SetProperty, DATA_RANGE=data_range
         self->Update
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
         self->GetProperty, POSITION=position
         oframe = self.animation->GetFrame(POSITION=position)
         oframe[0]->GetProperty, KEYWORDS=keywords
         data_dims = size(keywords.data_values, /DIMENSIONS)
         ;; Call REFORM so that XVAREDIT cannot modify values
         xvaredit, reform(keywords.data_values), GROUP=self.base, $
                   X_SCROLL_SIZE=(data_dims[0] < 8), $
                   Y_SCROLL_SIZE=(data_dims[1] < 30)
         return, 0
      end

      else: return, self->MGH_Datamator::EventMenuBar(event)

   endcase

end

; Mgh_Roms_Movie_Hslice::ExportData
;
pro Mgh_Roms_Movie_Hslice::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::ExportData, values, labels

   self->GetProperty, ANIMATION=animation, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[0]->GetProperty, KEYWORDS=keywords

   labels = [labels, 'Slice Data']
   values = [values, ptr_new(keywords.data_values)]

end

; Mgh_Roms_Movie_Hslice::PickReport
;
pro Mgh_Roms_Movie_Hslice::PickReport, pos, LUN=lun

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

; Mgh_Roms_Movie_Hslice::ExportAnimationDataToNcFile
;
pro Mgh_Roms_Movie_Hslice::ExportAnimationDataToNcFile, file

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Query the appropriate objects for info

   self->GetProperty, $
      LONLAT=lonlat, MAP_STRUCTURE=map_structure, VARIABLE=variable

   self.animation->GetProperty, $
      GRAPHICS_TREE=ograph, N_FRAMES=n_frames

   self.plane->GetProperty, $
      DATAX=datax, DATAY=datay, STYLE=style

   self.animator->GetPlayBack, $
      RANGE=play_range, USE_RANGE=play_use_range

   ;; Sort out grid data

   if style eq 0 then begin
      datax = mgh_stagger(temporary(datax), DELTA=[-1,-1])
      datay = mgh_stagger(temporary(datay), DELTA=[-1,-1])
   endif

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

   onc->AttAdd, /GLOBAL, 'title', 'ROMS hslice scalar data'

   onc->AttAdd, /GLOBAL, 'history', $
      'Generated by routine Mgh_Roms_Movie_Hslice::ExportAnimationDataToNcFile at '+ $
      mgh_dt_string(mgh_dt_now())

   onc->AttAdd, /GLOBAL, 'variable', variable

   onc->DimAdd, 'xi', dim[0]
   onc->DimAdd, 'eta', dim[1]
   onc->DimAdd, 'time'

   x_name = lonlat ? 'lon' : 'x'
   y_name = lonlat ? 'lat' : 'y'
   v_name = mgh_str_vanilla(variable)

   onc->VarAdd, x_name, ['xi','eta'], /DOUBLE
   onc->VarAdd, y_name, ['xi','eta'], /DOUBLE
   onc->VarAdd, v_name, ['xi','eta','time'], /FLOAT

   onc->VarPut, x_name, datax
   onc->VarPut, y_name, datay

   ;; Work through frames, retrieving and write data

   n_rec = 1 + (range[1]-range[0])/stride

   fmt ='(%"Writing %d frames of %d x %d data to netCDF file %s")'
   message, /INFORM, string(n_rec, dim, file, FORMAT=fmt)

   for r=0,n_rec-1 do begin

      pos = range[0]+stride*r

      oframe = self.animation->GetFrame(POSITION=pos)
      oframe[0]->GetProperty, KEYWORDS=keywords

      onc->VarPut, v_name, keywords.data_values, OFFSET=[0,0,r]

   endfor

   obj_destroy, onc

   fmt ='(%"Finished saving netCDF file %s")'
   message, /INFORM, string(file, FORMAT=fmt)

end

pro Mgh_Roms_Movie_Hslice__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
      {mgh_roms_movie_hslice, inherits MGH_Datamator, $
       variable: '', lonlat: !false, map_structure: ptr_new(), $
       history_file: obj_new(), $
       bar: obj_new(), plane: obj_new(), contour: obj_new()}

end
