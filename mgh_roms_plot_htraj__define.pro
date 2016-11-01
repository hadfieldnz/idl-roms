;+
; CLASS NAME:
;   MGH_ROMS_Plot_Htraj
;
; PURPOSE:
;   This procedure displays ROMS horizontal float trajectories in (x,y) or
;   (lon,lat) space
;
; CALLING SEQUENCE:
;   mgh_new, 'MGH_ROMS_Plot_Htraj', ffile, hfile
;
; SUPERCLASS:
;   MGH_Window
;
; INIT ARGUMENTS:
;  ffile
;    Synonym for FLOAT_FILE property.
;
;  hfile
;    Synonym for HISTORY_FILE property.
;
;  ifile
;    Synonym for FLTIN_FILE property.
;
; PROPERTIES:
;  ETA_RANGE (Init)
;    A 2-element integer vector specifying the grid points to be
;    plotted in the eta direction. The range is specified relative to
;    the rho grid.  Negative values are taken to be offsets from the
;    "northern" edge of the grid. The default is [0,-1] which is
;    equivalent to [0,dim(eta_rho)-1].
;
;  FLOAT_FILE (Init, Get)
;    A reference to a ROMS float-file object.
;
;  FLOAT_RANGE (Init)
;  FLOAT_STRIDE (Init)
;    Range and stride specifying a subset of floats to be
;    plotted. These are ignored if FLOATS is specified.
;
;  FLOATS (Init)
;    An integer vector specifying the floats to be plotted
;
;  FLTIN_FILE (Init, Get)
;    The name of a ROMS floats-input file.
;
;  GRAPH_PROPERTIES (Init)
;    A structure containing keywords to be passed to the graph.
;
;  HISTORY_FILE (Init, Get)
;    A reference to a ROMS history-file object.
;
;   MAP_STRUCTURE
;     If a map structure is supplied via this keyword, and if there
;     are longitude and latitude data in the file, then the positions
;     are converted with MAP_PROJ_FORWARD before display.
;
;  RECORD_RANGE (Init)
;  RECORD_STRIDE (Init)
;    Range and stride specifying a subset of records to be
;    plotted.
;
;  SYMBOL_PROPERTIES (Init)
;    A structure containing keywords to be passed to the plotting
;    symbol.
;
;  TIME_RANGE (Init)
;    This keyword can be used to clip the trajectory to a specified
;    time interval. Times are measured in days from the release time
;    of each float (approximated as the first time at which the float
;    has an unbounded position).
;
;  XAXIS_PROPERTIES (Init)
;  YAXIS_PROPERTIES (Init)
;    Structures containing keywords to be passed to the X & Y axis.
;
;  XI_RANGE (Init)
;    A 2-element integer vector specifying the grid points to be plotted
;    in the xi direction. The range is specified relative to the rho grid.
;    Negative values are taken to be offsets from the "eastern" edge of
;    the grid. The default is [0,-1] which is equivalent to [0,dim(xi_rho)-1].
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2002-02:
;     Written.
;   Mark Hadfield, 2009-09:
;     Added support for map projections via the MAP_STRUCTURE
;     keyword.
;   Mark Hadfield, 2010-11:
;     - The procedure now gets all trajectories then selects from them.
;     - Removed code relating to destroying float-file and history-file
;       objects: no longer necessary with IDL 8 (and never that
;       important),
;     - Corrected bug: failure to account for non-default XI_RANGE &
;       ETA_RANgE when calculating release positions.
;   Mark Hadfield, 2016-10:
;     - Fixed bugs in SHOW_ORIGIN and TIME_RANGE functionality.
;-
function MGH_ROMS_Plot_Htraj::Init, ffile, hfile, ifile, $
     BATH_PROPERTIES=bath_properties, $
     FLOAT_FILE=float_file, $
     FLOAT_RANGE=float_range, $
     FLOAT_STRIDE=float_stride, $
     FLOATS=floats, $
     FLTIN_FILE=fltin_file, $
     MAP_STRUCTURE=map_structure, $
     SHOW_ORIGIN=show_origin, $
     SHOW_SYMBOL=show_symbol, $
     TIME_RANGE=time_range, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
     X_RANGE=x_range, Y_RANGE=y_range, $
     HISTORY_FILE=history_file, $
     GRAPH_PROPERTIES=graph_properties, $
     LAND_PROPERTIES=land_properties, $
     PLOT_PROPERTIES=plot_properties, $
     SYMBOL_PROPERTIES=symbol_properties, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process float-file argument

   if n_elements(float_file) eq 0 && n_elements(ffile) gt 0 then $
        float_file = ffile

   if n_elements(float_file) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'float_file'

   case !true of
      isa(float_file, 'STRING'): begin
         self.float_file = obj_new('MGHromsFloat', float_file)
      end
      isa(float_file, 'OBJREF'): begin
         self.float_file = float_file
      end
   endcase
   oflt = self.float_file

   ;; Process history-file argument

   if n_elements(history_file) eq 0 && n_elements(hfile) gt 0 then $
        history_file = hfile

   if n_elements(history_file) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'history_file'

   case !true of
      isa(history_file, 'STRING'): begin
         self.history_file = obj_new('MGHromsHistory', history_file)
      end
      isa(history_file, 'OBJREF'): begin
         self.history_file = history_file
      end
   endcase
   ohis = self.history_file

   ;; Process floats-input file argument

   if n_elements(fltin_file) eq 0 && n_elements(ifile) gt 0 then $
        fltin_file = ifile

   ;; Process miscellaneous arguments

   if n_elements(time_range) gt 0 && n_elements(time_range) ne 2 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_WRGNUMELEM', 'time_range'

   ;; Get number of floats and resolve float-selection arguments

   n_float = oflt->DimInfo('drifter', /DIMSIZE)

   mgh_resolve_indices, n_float, float_range, float_stride, floats

   ;; Get grid size and resolve grid-subset arguments

   dim_rho = ohis->DimRho()

   if n_elements(xi_range) eq 0 then xi_range = [0,-1]
   if xi_range[0] lt 0 then xi_range[0] += dim_rho[0]
   if xi_range[1] lt 0 then xi_range[1] += dim_rho[0]

   if n_elements(eta_range) eq 0 then eta_range = [0,-1]
   if eta_range[0] lt 0 then eta_range[0] += dim_rho[1]
   if eta_range[1] lt 0 then eta_range[1] += dim_rho[1]

   xra0 = xi_range[0]
   xra1 = xi_range[1]
   xran = xra1-xra0+1
   era0 = eta_range[0]
   era1 = eta_range[1]
   eran = era1-era0+1

   ;; Are longitude & latitude data available?

   self.lonlat = ohis->HasVar('lon_rho') and ohis->HasVar('lat_rho')

   ;; Are we using a map projection?

   if n_elements(map_structure) gt 0 then begin
      if ~ self.lonlat then $
           message, 'Cannot use map structure without lon, lat data'
      self.map_structure = ptr_new(map_structure)
   endif

   use_map_structure = ptr_valid(self.map_structure)

   ;; Get x & y positions @ RHO points. This is used for walls and
   ;; land mask, and also possibly for float positions. Convert to map
   ;; projection if appropriate

   o = [xra0,era0]  &  c = [xran,eran]

   x_rho = ohis->VarGet(self.lonlat ? 'lon_rho' : 'x_rho', OFFSET=o, COUNT=c)
   y_rho = ohis->VarGet(self.lonlat ? 'lat_rho' : 'y_rho', OFFSET=o, COUNT=c)

   if keyword_set(use_map_structure) then begin
      xy = map_proj_forward(x_rho, y_rho, MAP_STRUCTURE=*self.map_structure)
      x_rho = reform(xy[0,*], size(x_rho, /DIMENSIONS))
      y_rho = reform(xy[1,*], size(y_rho, /DIMENSIONS))
      mgh_undefine, xy
   endif

   ;; If we are going to clip the trajectories in time, establish the
   ;; appropriate variable name and retrieve the data

   if n_elements(time_range) gt 0 then begin
      case !true of
         oflt->HasVar('Xgrid') && oflt->HasAtt('Xgrid','time'): $
              time_var = oflt->AttGet('Xgrid','time')
         oflt->HasVar('Ygrid') && oflt->HasAtt('Ygrid','time'): $
              time_var = oflt->AttGet('Ygrid','time')
         else: $
              time_var = 'ocean_time'
      endcase
      time = oflt->VarGet(time_var, /AUTOSCALE)
   endif

   ;; Now get float x,y data

   if self.lonlat then begin
      xvar = 'lon'
      yvar = 'lat'
   endif else begin
      xvar = 'x'
      yvar = 'y'
   endelse

   ;; Get float position data (currently all records)

   xgrid = oflt->VarGet('Xgrid', AUTOSCALE=0)
   xgrid = xgrid[floats,*]

   ygrid = oflt->VarGet('Ygrid', AUTOSCALE=0)
   ygrid = ygrid[floats,*]

   ;; Float grid-relative positions should be in the range
   ;; [0.5,dim_rho-1.5]; values outside this range imply the float
   ;; has not been released yet or has become unbounded. (The xgrid
   ;; & ygrid values in this case are usually 1.E35, but for an MPI
   ;; run they may be zero.)

   l_bound = where(xgrid ge 0.5 and xgrid le dim_rho[0]-1.5 and $
                   ygrid ge 0.5 and ygrid le dim_rho[1]-1.5, n_bound, $
                   COMPLEMENT=l_unbound, NCOMPLEMENT=n_unbound)

   if n_bound gt 0 then begin
      xgrid[l_bound] -= xra0
      ygrid[l_bound] -= era0
   end

   if n_unbound gt 0 then begin
      xgrid[l_unbound] = !values.f_nan
      ygrid[l_unbound] = !values.f_nan
   endif

   x = mgh_interpolate(x_rho, xgrid, ygrid, MISSING=!values.f_nan)
   y = mgh_interpolate(y_rho, xgrid, ygrid, MISSING=!values.f_nan)

   mgh_undefine, xgrid, ygrid

   ;; If information about the float-release time or location is
   ;; required, get them from the floats-input file

   if keyword_set(show_origin) || n_elements(time_range) gt 0 then begin

      flt = mgh_roms_read_fltin(fltin_file)

      if total(flt.src.n, /INTEGER) ne flt.n_float then $
           message, 'Mismatch in number of floats within input file'

;      if flt.n_float ne n_float then $
;           message, 'Mismatch in number of floats between input and output files'

      ;; Set up of release-location arrays. For now, ignore the case
      ;; where the number of floats is zero

      x0 = dblarr(n_float)
      y0 = dblarr(n_float)
      t0 = dblarr(n_float)

      ;; For each source (i.e. each entry in the float-release data)
      ;; establish release locations

      n0 = 0
      for s=0,flt.n_src-1 do begin
         n1 = n0 + flt.src[s].n - 1
         fx0 = flt.src[s].fx0
         if flt.src[s].fdx ne 0 then $
              fx0 += flt.src[s].fdx*dindgen(flt.src[s].n)
         fy0 = flt.src[s].fy0
         if flt.src[s].fdy ne 0 then $
              fy0 += flt.src[s].fdy*dindgen(flt.src[s].n)
         ft0 = flt.src[s].ft0
         if flt.src[s].fdt ne 0 then $
              ft0 += flt.src[s].fdt*dindgen(flt.src[s].n)
         case flt.src[s].c of
            0: begin
               x0[n0:n1] = mgh_interpolate(x_rho, fx0-xra0, fy0-era0, MISSING=!values.f_nan)
               y0[n0:n1] = mgh_interpolate(y_rho, fx0-xra0, fy0-era0, MISSING=!values.f_nan)
               mgh_undefine, fx0, fy0
            end
            1: begin
               x0[n0:n1] = temporary(fx0)
               y0[n0:n1] = temporary(fy0)
            end
         endcase
         t0[n0:n1] = temporary(ft0)
         n0 = n1 + 1
      endfor

      x0 = x0[floats]
      y0 = y0[floats]

      t0 = t0[floats]

      if keyword_set(use_map_structure) then begin
         xy = map_proj_forward(x0, y0, MAP_STRUCTURE=*self.map_structure)
         x0 = reform(xy[0,*])
         y0 = reform(xy[1,*])
         mgh_undefine, xy
      endif

      ;; Times in the floats-input file differ from times in the
      ;; floats-output file by the simulation parameter DSTART.

      t0 += oflt->VarGet('dstart')

   endif

   ;; Default graph aspect ratio, can be overridden via
   ;; GRAPH_PROPERTIES

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(x_rho)
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(y_rho)

   aspect = mgh_aspect(x_range, y_range, LONLAT=self.lonlat && (~ use_map_structure))
   aspect = (aspect > 0.4) < 1.5

   ;; Create base graph

   ograph = obj_new('MGHgrGraph2D', ASPECT=aspect, NAME='ROMS float trajectory', $
                    _STRICT_EXTRA=graph_properties)

   ograph->GetProperty, DELTAZ=deltaz

   ograph->NewFont

   ograph->NewMask

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
   ;; wall, then trim the interior points off xr and yr.

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
           STYLE=2, COLOR=[127,127,127], NAME='Wall west'
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
           STYLE=2, COLOR=[127,127,127], NAME='Wall south'
   endif

   if walls[2] then begin
      ;; Eastern wall
      xw = mgh_stagger(x_rho[-2:-1,*], DELTA=[1,1])
      yw = mgh_stagger(y_rho[-2:-1,*], DELTA=[1,1])
      xw = xw[1:2,*]
      yw = yw[1:2,*]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', $
           DATAX=temporary(xw), DATAY=temporary(yw), DATAZ=temporary(zw), $
           STYLE=2, COLOR=[127,127,127], NAME='Wall east'
   endif

   if walls[3] then begin
      ;; Northern wall
      xw = mgh_stagger(x_rho[*,-2:-1], DELTA=[1,1])
      yw = mgh_stagger(y_rho[*,-2:-1], DELTA=[1,1])
      xw = xw[*,1:2]
      yw = yw[*,1:2]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', $
           DATAX=temporary(xw), DATAY=temporary(yw), DATAZ=temporary(zw), $
           STYLE=2, COLOR=[127,127,127], NAME='Wall north'
   endif

   ;; Draw land mask.

   if ohis->HasVar('mask_rho') then begin

      mask_rho = ohis->VarGet('mask_rho', OFFSET=[xra0,era0], COUNT=[xran,eran])

      ograph->NewAtom, 'MGHgrColorPlane', RESULT=oland, $
           STYLE=0, DEFAULT_COLOR=mgh_color('grey'), ZVALUE=-2*deltaz, $
           MISSING_POINTS=round(mask_rho), $
           DATAX=mgh_stagger(x_rho, DELTA=[1,1]), $
           DATAY=mgh_stagger(y_rho, DELTA=[1,1]), $
           NAME='Land mask', _STRICT_EXTRA=land_properties

      ;; Just playing really
      if oland->QueryProperty('STYLE') then $
           oland->SetPropertyAttribute, 'STYLE', SENSITIVE=0

   endif

   ;; Draw bathymetry contours.

   if ohis->HasVar('h') then begin

      hr = ohis->VarGet('h', OFFSET=[xra0,era0], COUNT=[xran,eran])

      ograph->NewAtom, 'IDLgrContour', RESULT=obath, $
           DATA=hr, GEOMX=x_rho, GEOMY=y_rho, GEOMZ=-4*deltaz, /PLANAR, $
           C_COLOR=mgh_color(['blue','blue']), $
           NAME='Bathymetry', _STRICT_EXTRA=bath_properties

   endif

   ;; If SHOW_SYMBOL keyword is set, create a symbol to show float data points

   if keyword_set(show_symbol) then begin
      osym = ograph->NewSymbol(0, N_VERTICES=3, NORM_SIZE=0.001, $
                               COLOR=mgh_color('red'), _STRICT_EXTRA=symbol_properties)
   endif

   ;; Step through floats, extracting and plotting trajectories

   for f=0,n_elements(floats)-1 do begin

      if n_elements(time_range) gt 0 then begin
         rr = mgh_subset(time-t0[f], time_range, EMPTY=empty)
         if ~ empty then begin
            ograph->NewAtom, 'IDLgrPlot', x[f,rr[0]:rr[1]], y[f,rr[0]:rr[1]], $
                 NAME=string(FORMAT='(%"Float %d")', floats[f]), $
                 COLOR=mgh_color('red'), SYMBOL=osym, $
                 _STRICT_EXTRA=plot_properties
         endif
      endif else begin
         ograph->NewAtom, 'IDLgrPlot', x[f,*], y[f,*], $
              NAME=string(FORMAT='(%"Float %d")', floats[f]), $
              COLOR=mgh_color('red'), SYMBOL=osym, $
              _STRICT_EXTRA=plot_properties
      endelse

   endfor

   ;; Optionally plot a symbol at the origin (release-point) of each
   ;; trajectory.  Currently all release points are shown, however it
   ;; would be possible (and sometimes very much more efficient) to
   ;; remove duplicates. This post suggests how this could be done:
   ;;   http://tinyurl.com/pajpg5

   if keyword_set(show_origin) then begin
      ograph->NewAtom, 'IDLgrPlot', x0, y0, ZVALUE=20*deltaz, LINESTYLE=6, $
           SYMBOL=ograph->NewSymbol(0, /FILL, NORM_SIZE=0.015, COLOR=mgh_color('dark green'))
   endif

   ma = ['Magnify','Translate','Context']

   ok = self->MGH_Window::Init(GRAPHICS_TREE=ograph, $
                               MOUSE_ACTION=ma,  _STRICT_EXTRA=extra)

   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   ;; Finalise plot & return

   self->Finalize, 'MGH_ROMS_Plot_Htraj'

   return, 1

end

; MGH_ROMS_Plot_Htraj::Cleanup
;
pro MGH_ROMS_Plot_Htraj::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Window::Cleanup

end

; MGH_ROMS_Plot_Htraj::GetProperty
;
pro MGH_ROMS_Plot_Htraj::GetProperty, $
     ALL=all, FLOAT_FILE=float_file, HISTORY_FILE=history_file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::GetProperty, ALL=all, _STRICT_EXTRA=extra

   float_file = self.float_file
   history_file = self.history_file

   if arg_present(all) then $
        all = create_struct(all, $
                            'float_file', float_file, 'history_file', history_file)

end

; MGH_ROMS_Plot_Htraj::SetProperty
;
pro MGH_ROMS_Plot_Htraj::SetProperty, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::SetProperty, _STRICT_EXTRA=extra

end

; MGH_ROMS_Plot_Htraj::About
;
;   Print information about the window and its contents
;
pro MGH_ROMS_Plot_Htraj::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::About, lun

   if obj_valid(self.float_file) then begin
      printf, lun, self, ': the float file sequence is '+ $
              mgh_obj_string(self.float_file)+'. Its files are:'
      self.float_file->GetProperty, FILE_NAME=file_name
      print, file_name
   endif

   if obj_valid(self.history_file) then begin
      printf, lun, self, ': the history file sequence is '+ $
              mgh_obj_string(self.history_file)+'. Its files are:'
      self.history_file->GetProperty, FILE_NAME=file_name
      print, file_name
   endif

end

pro MGH_ROMS_Plot_Htraj::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::BuildMenuBar

end


function MGH_ROMS_Plot_Htraj::Event, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case widget_info(event.id, /UNAME) of

      else: return, self->MGH_Window::Event(event)

   endcase

end

pro MGH_ROMS_Plot_Htraj::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::ExportData, values, labels

   self->GetProperty, $
        FLOAT_FILE=float_file, HISTORY_FILE=history_file

   labels = [labels, 'Float Object', 'History Object']
   values = [values, ptr_new(float_file), ptr_new(history_file)]

end

pro MGH_ROMS_Plot_Htraj__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {MGH_ROMS_Plot_Htraj, inherits MGH_Window, $
         float_file: obj_new(), history_file: obj_new(), $
         lonlat: !false, map_structure: ptr_new()}

end
