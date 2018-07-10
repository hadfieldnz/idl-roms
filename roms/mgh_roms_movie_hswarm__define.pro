;+
; CLASS NAME:
;   MGH_ROMS_Movie_Hswarm
;
; PURPOSE:
;   This procedure generates and displays an animated graph showing
;   ROMS float locations on a horizontal surface (H-slice).
;
; CALLING SEQUENCE:
;   mgh_new, 'MGH_ROMS_Movie_Hswarm', ffile, hfile
;
; SUPERCLASS:
;   MGH_Datamator
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
;  DEPTH_SELeCT (Init)
;    A 2-element numeric vector specifying a range of float depths
;    (+ve upward). Floats initially outside this range are excluded.
;
;  ETA_RANGE (Init)
;    A 2-element integer vector specifying the grid points to be plotted
;    in the eta direction. The range is specified relative to the rho grid.
;    Negative values are taken to be offsets from the "northern" edge of
;    the grid. The default is [0,-1] which is equivalent to [0,dim(eta_rho)-1].
;
;  FLTIN_FILE (Init, Get)
;    The name of a ROMS floats-input file.
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
;    plotted. These are ignored if RECORDS is specified.
;
;  RECORDS (Init)
;    An integer vector specifying the records to be plotted
;
;  SYMBOL_PROPERTIES (Init)
;    A structure containing keywords to be passed to the plotting
;    symbol.
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
;   Mark Hadfield, 2001-03:
;     Written.
;   Mark Hadfield, 2001-07:
;     Now a subclass of MGH_Datamator. Updated for IDL 5.5.
;   Mark Hadfield, 2009-07:
;     Added support for map projections via the MAP_STRUCTURE
;     keyword.
;   Mark Hadfield, 2011-08:
;     Removed FLOAT_DESTROY and HISTORY_DESTROY keywords.
;   Mark Hadfield, 2011-10:
;     - Copied SHOW_TIME & TITLE functionality from the
;       mgh_roms_movie_hslice class, the former allowing optional
;       date-time formatting.
;     - Added SHOW_GSHHS & SHOW_LOGO functionality.
;   Mark Hadfield, 2014-04:
;     Further enhancements to SHOW_TIME functionality.
;-
function mgh_roms_movie_hswarm::Init, ffile, hfile, $
     DEPTH_SELECT=depth_select, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
     X_RANGE=x_range, Y_RANGE=y_range, $
     FLOAT_FILE=float_file, $
     FLOAT_RANGE=float_range, $
     FLOAT_STRIDE=float_stride, $
     FLOATS=floats, $
     FLTIN_FILE=fltin_file, $
     HISTORY_FILE=history_file, $
     MAP_STRUCTURE=map_structure, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     RECORDS=records, SHOW_LOGO=show_logo, $
     SHOW_TIME=show_time, SHOW_GSHHS=show_gshhs, TITLE=title, $
     BATH_PROPERTIES=bath_properties, $
     GRAPH_PROPERTIES=graph_properties, $
     LAND_PROPERTIES=land_properties, $
     SYMBOL_PROPERTIES=symbol_properties, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process float file argument

   if n_elements(float_file) eq 0 && n_elements(ffile) gt 0 then $
        float_file = ffile

   if n_elements(float_file) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'float_file'

   case size(float_file, /TNAME) of
      'STRING': begin
         self.float_file = obj_new('MGHromsFloat', float_file)
      end
      'OBJREF': begin
         self.float_file = float_file
      end
   endcase
   oflt = self.float_file

   ;; Process history file argument

   if n_elements(history_file) eq 0 && n_elements(hfile) gt 0 then $
        history_file = hfile

   if n_elements(history_file) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'history_file'

   case size(history_file, /TNAME) of
      'STRING': begin
         self.history_file = obj_new('MGHromsHistory', history_file)
      end
      'OBJREF': begin
         self.history_file = history_file
      end
   endcase
   ohis = self.history_file

   ;; Process floats-input file argument

   if n_elements(fltin_file) eq 0 && n_elements(ifile) gt 0 then $
        fltin_file = ifile

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

   self.lonlat = ohis->HasVar('lon_rho') && ohis->HasVar('lat_rho')

   ;; Are we using a map projection?

   if n_elements(map_structure) gt 0 then begin
      if ~ self.lonlat then $
           message, 'Cannot use map structure without lon, lat data'
      self.map_structure = ptr_new(map_structure)
   endif

   use_map_structure = ptr_valid(self.map_structure)

   if n_elements(title) eq 0 then title = ''
   if n_elements(show_time) eq 0 then show_time = 1B

   ;; Get x & y positions @ RHO points. This is used for float
   ;; positions, walls and land mask. Convert to map projection if
   ;; appropriate

   o = [xra0,era0]  &  c = [xran,eran]

   x_rho = ohis->VarGet(self.lonlat ? 'lon_rho' : 'x_rho', OFFSET=o, COUNT=c)
   y_rho = ohis->VarGet(self.lonlat ? 'lat_rho' : 'y_rho', OFFSET=o, COUNT=c)

   if use_map_structure then begin
      xy = map_proj_forward(x_rho, y_rho, MAP_STRUCTURE=*self.map_structure)
      x_rho = reform(xy[0,*], size(x_rho, /DIMENSIONS))
      y_rho = reform(xy[1,*], size(y_rho, /DIMENSIONS))
      mgh_undefine, xy
   endif

   ;; Establish time variable name, dimension name and units associated with
   ;; the float data

   case !true of
      oflt->HasVar('Xgrid') && oflt->HasAtt('Xgrid','time'): $
           time_var = oflt->AttGet('Xgrid','time')
      oflt->HasVar('Ygrid') && oflt->HasAtt('Ygrid','time'): $
           time_var = oflt->AttGet('Ygrid','time')
      else: $
           time_var = 'ocean_time'
   endcase

   time_dim = (oflt->VarDimNames(time_var))[0]

   if oflt->HasAtt(time_var, 'units') then begin
      time_units = mgh_dt_units(oflt->AttGet(time_var, 'units'))
   endif else begin
      time_units = {scale: 1}
   endelse

   ;; Get number of records in float file and resolve record-selection
   ;; arguments

   n_time = oflt->DimInfo(time_dim, /DIMSIZE)

   mgh_resolve_indices, n_time, record_range, record_stride, records

   ;; Default graph aspect ratio, can be overridden via
   ;; GRAPH_PROPERTIES

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(x_rho)
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(y_rho)

   aspect = mgh_aspect(x_range, y_range, LONLAT=self.lonlat && (~ use_map_structure))
   aspect = (aspect > 0.4) < 1.5

   ;; Create base graph

   ograph = obj_new('MGHgrGraph2D', ASPECT=aspect, NAME='ROMS float animation', $
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

   if keyword_set(show_logo) then begin

      logo = read_png(!mgh_logo_file)

      ograph->GEtProperty, VIEWPLANE_RECT=vrect

      ograph->NewAtom, 'IDLgrImage', logo, $
              XAXIS=0, YAXIS=0, $
              LOCATION=[vrect[0]+0.01*vrect[2],vrect[1]+0.93*vrect[3]], $
              DIMENSIONS=[0.24*vrect[2],0.06*vrect[3]]

   endif

   ;; Add walls. For each wall we extract the 2 rows/columns on each
   ;; side of the physical boundary from the x_rho & y_rho arrays into
   ;; variables xr & yr, call mgh_stagger to get xw & yw along the
   ;; wall, then trim the interior points off xr and yr.

   walls = ohis->GetWalls()

   if walls[0] && (xra0 eq 0) then begin
      ;; Western wall
      xw = mgh_stagger(x_rho[0:1,*], DELTA=[1,1])
      yw = mgh_stagger(y_rho[0:1,*], DELTA=[1,1])
      xw = xw[0:1,*]
      yw = yw[0:1,*]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', DATAX=xw, DATAY=yw, DATAZ=zw, $
           STYLE=2, COLOR=[127,127,127], NAME='Wall west'
   endif

   if walls[1] && (era0 eq 0) then begin
      ;; Southern wall
      xw = mgh_stagger(x_rho[*,0:1], DELTA=[1,1])
      yw = mgh_stagger(y_rho[*,0:1], DELTA=[1,1])
      xw = xw[*,0:1]
      yw = yw[*,0:1]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', DATAX=xw, DATAY=yw, DATAZ=zw, $
           STYLE=2, COLOR=[127,127,127], NAME='Wall south'
   endif

   if walls[2] && (xra1 eq dim_rho[0]-1) then begin
      ;; Eastern wall
      xw = mgh_stagger(x_rho[xran-2:xran-1,*], DELTA=[1,1])
      yw = mgh_stagger(y_rho[xran-2:xran-1,*], DELTA=[1,1])
      xw = xw[1:2,*]
      yw = yw[1:2,*]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', DATAX=xw, DATAY=yw, DATAZ=zw, $
           STYLE=2, COLOR=[127,127,127], NAME='Wall east'
   endif

   if walls[3] && (era1 eq dim_rho[1]-1) then begin
      ;; Northern wall
      xw = mgh_stagger(x_rho[*,eran-2:eran-1], DELTA=[1,1])
      yw = mgh_stagger(y_rho[*,eran-2:eran-1], DELTA=[1,1])
      xw = xw[*,1:2]
      yw = yw[*,1:2]
      zw = make_array(DIMENSION=size(xw, /DIMENSIONS))
      ograph->NewAtom, 'IDLgrSurface', DATAX=xw, DATAY=yw, DATAZ=zw, $
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

   ;; Create symbols to represent each float. Make them reasonably simple to
   ;; avoid slowing the animation. Make two sets: red, to represent floats
   ;; in the water column and blue, to represent floats on the bottom (Zgrid eq 0).

   if n_elements(symbol_properties) gt 0 then begin
     n_sym = n_elements(symbol_properties)
     osym_free = objarr(n_sym)
     osym_stuck = objarr(n_sym)
     for i=0,n_sym-1 do begin
       osym_free[i] = ograph->NewSymbol(0, N_VERTICES=6, COLOR=mgh_color('red'), NORM_SIZE=0.01, _STRICT_EXTRA=symbol_properties[i])
       osym_stuck[i] = ograph->NewSymbol(0, N_VERTICES=6, COLOR=mgh_color('dark green'), NORM_SIZE=0.01, _STRICT_EXTRA=symbol_properties[i])
     endfor
   endif else begin
     osym_free = ograph->NewSymbol(0, N_VERTICES=6, COLOR=mgh_color('red'), NORM_SIZE=0.01)
     osym_stuck = ograph->NewSymbol(0, N_VERTICES=6, COLOR=mgh_color('dark green'), NORM_SIZE=0.01)
   endelse

   ;; Add a title and two invisible polylines to be animated

   ograph->NewTitle, '', RESULT=otitle

   ograph->NewAtom, 'IDLgrPlot', LINESTYLE=6, SYMBOL=osym_free, NAME='Free floats', RESULT=oplt_free
   ograph->NewAtom, 'IDLgrPlot', LINESTYLE=6, SYMBOL=osym_stuck, NAME='Stuck floats', RESULT=oplt_stuck

   ;; Create an animator window to display and manage the movie.

   ok = self->MGH_Datamator::Init(CHANGEABLE=0, GRAPHICS_TREE=ograph, $
                                  _STRICT_EXTRA=extra)
   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Datamator'

   ;; Step through times in the netCDF file, generating new frames &
   ;; plotting data

   oframe = objarr(3)

   for r=0,n_elements(records)-1 do begin

      rec = records[r]

      ;; Check to see if the user has selected the "Finish Loading" menu item

      if self->Finished() then break

      ;; Get time

      time = oflt->VarGet(time_var, OFFSET=records[r], $
                          COUNT=[1], AUTOSCALE=0)*time_units.scale

      ;; Get float position data

      xgrid = oflt->VarGetFloat('Xgrid', FLOATS=floats, RECORDS=rec)
      ygrid = oflt->VarGetFloat('Ygrid', FLOATS=floats, RECORDS=rec)

      ;; Float grid-relative positions should be in the range
      ;; [0.5,dim_rho-1.5]; values outside this range imply the float
      ;; has not been released yet or has become unbounded. (The xgrid
      ;; & ygrid values in this case are usually 1.E37, but for an MPI
      ;; run they may be zero.)

      l_bound = where(xgrid ge 0.5 and xgrid le dim_rho[0]-1.5 and $
                      ygrid ge 0.5 and ygrid le dim_rho[1]-1.5, n_bound, $
                      COMPLEMENT=l_unbound, NCOMPLEMENT=n_unbound)

      x = replicate(!values.f_nan, n_elements(xgrid))
      y = replicate(!values.f_nan, n_elements(xgrid))
      
      if n_bound gt 0 then begin
         x[l_bound] = $
              interpolate(x_rho, xgrid[l_bound]-xra0, ygrid[l_bound]-era0)
         y[l_bound] = $
              interpolate(y_rho, xgrid[l_bound]-xra0, ygrid[l_bound]-era0)
      endif

      ;; Determine which floats are sitting on the bottom

      stuck = oflt->VarGetFloat('Zgrid', FLOATS=floats, RECORDS=rec) eq 0

      ;; Optionally select for floats in a specified depth range. This check
      ;; is independent of the "stuck" check.

      if n_elements(depth_select) gt 0 then begin
         z = oflt->VarGetFloat('depth', FLOATS=floats, RECORDS=rec, AUTOSCALE=0)
         l_good = where(z ge depth_select[0] and z le depth_select[1], n_good)
         if n_good gt 0 then begin
            x = x[l_good]
            y = y[l_good]
            stuck = stuck[l_good]
         endif else begin
            x = [!values.f_nan]
            y = [!values.f_nan]
            stuck = [!values.f_nan]
         endelse
      endif
      
      ;; Plot stuck and free floats separately

      l_stuck = where(stuck, n_stuck, $
         COMPLEMENT=l_free, NCOMPLEMENT=n_free)

      if n_free gt 0 then begin
         oframe[0] = obj_new('MGH_Command', OBJECT=oplt_free, 'SetProperty', DATAX=x[l_free], DATAY=y[l_free])
      endif else begin
         oframe[0] = obj_new('MGH_Command', OBJECT=oplt_free, 'SetProperty', DATAX=[!values.f_nan], DATAY=[!values.f_nan])
      endelse

      if n_stuck gt 0 then begin
         oframe[1] = obj_new('MGH_Command', OBJECT=oplt_stuck, 'SetProperty', DATAX=x[l_stuck], DATAY=y[l_stuck])
      endif else begin
         oframe[1] = obj_new('MGH_Command', OBJECT=oplt_stuck, 'SetProperty', DATAX=[!values.f_nan], DATAY=[!values.f_nan])
      endelse

      ;; Update title with time or date, if applicable

      if show_time gt 0 then begin
         case show_time of
            1: ttt = string(FORMAT='(%"%0.3f days")', time)
            2: ttt = string(FORMAT='(%"%s (%0.3f days)")', mgh_dt_string(time+time_units.offset), time)
         endcase
         if strlen(title) gt 0 then $
              ttt = string(FORMAT='(%"%s: %s")', title, ttt)
         oframe[2] = obj_new('MGH_Command', OBJECT=otitle, $
                             'SetProperty', STRINGS=temporary(ttt))
      endif

      ;; Add frame to animator & display

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; MGH_ROMS_Movie_Hswarm::Cleanup
;
pro MGH_ROMS_Movie_Hswarm::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Datamator::Cleanup

end

; MGH_ROMS_Movie_Hswarm::GetProperty
;
pro MGH_ROMS_Movie_Hswarm::GetProperty, $
     ALL=all, FLOAT_FILE=float_file, HISTORY_FILE=history_file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::GetProperty, ALL=all, _STRICT_EXTRA=extra

   float_file = self.float_file
   history_file = self.history_file

   if arg_present(all) then $
        all = create_struct(all, $
                            'float_file', float_file, 'history_file', history_file)

end

; MGH_ROMS_Movie_Hswarm::SetProperty
;
pro MGH_ROMS_Movie_Hswarm::SetProperty, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=extra

end

; MGH_ROMS_Movie_Hswarm::About
;
;   Print information about the window and its contents
;
pro MGH_ROMS_Movie_Hswarm::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::About, lun

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

pro MGH_ROMS_Movie_Hswarm::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::BuildMenuBar

end


function MGH_ROMS_Movie_Hswarm::Event, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case widget_info(event.id, /UNAME) of

      else: return, self->MGH_Datamator::Event(event)

   endcase

end

pro MGH_ROMS_Movie_Hswarm::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::ExportData, values, labels

   self->GetProperty, $
        ANIMATION=animation, FLOAT_FILE=float_file, $
        HISTORY_FILE=history_file, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[1]->GetProperty, KEYWORDS=keywords

   labels = [labels, 'Float Object', 'History Object', 'X,Y Data']
   values = [values, ptr_new(float_file), ptr_new(history_file), $
             ptr_new(complex(keywords.datax,keywords.datay))]

end

pro MGH_ROMS_Movie_Hswarm__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {MGH_ROMS_Movie_Hswarm, inherits MGH_Datamator, $
         float_file: obj_new(), history_file: obj_new(), $
         lonlat: !false, map_structure: ptr_new()}

end
