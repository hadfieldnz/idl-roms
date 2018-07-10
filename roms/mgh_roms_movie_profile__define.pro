;+
; CLASS NAME:
;   Mgh_Roms_Movie_Profile
;
; PURPOSE:
;   This class implements an animated sequence of graphs showing
;   a vertical profile through a SCRUM/ROMS field from a history OR station
;   file
;
; OBJECT CREATION SEQUENCE
;   mgh_new, 'Mgh_Roms_Movie_Profile', File, Position
;
; POSITIONAL ARGUMENTS:
;   file
;     The name of a ROMS history file OR a reference to an
;     MGHromsHistory or MGHromsStation object.
;
;   position
;     For a history file, a 2-element integer vector specifying the
;     [xi,eta] location of the profile; for a station file, a scalar
;     integer specifying the station number.
;
; KEYWORD ARGUMENTS:
;   DATA_MULTIPLIER (input, scalar numeric)
;     Number by which data values are multiplied before they are loaded into
;     the density surface. Default depends on the variable and is calculated by
;     MGH_ROMS_RESOLVE_DATA
;
;   DATA_RANGE (input, 2-element numeric)
;     Data range for the density surface. Default depends on the variable and
;     is calculated by MGH_ROMS_RESOLVE_DATA
;
;   DEPTH_RANGE
;     If supplied, this keyword specifies the height range of the plot
;
;   VARIABLE
;     The name of a variable in the file.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 1999-02:
;     Written.
;   Mark Hadfield, 2001-07:
;     Extensive modifications. Updated for IDL 5.5. Now can plot
;     profiles from history OR station files. Supports multiple profiles.
;   Mark Hadfield, 2011-09:
;     Removed the file_destroy property.
;   Mark Hadfield, 2011-11:
;     - Variables without a time dimension can be plotted. This
;       functionality has not been tested throughly.
;     - SHOW_TIME keyword added giving support for ISO date/times.
;     - RECORD_AVERAGE keyword giving support for...Oh, go on, guess!
;   Mark Hadfield, 2012-10:
;     - Variable name now specified by keyword argument.
;   Mark Hadfield, 2013-07:
;     - Added [XY]AXIS_PROPERTIES keywords.
;     - Fixed bug: DATA_MULTIPLIER not used.
;-
function Mgh_Roms_Movie_Profile::Init, file, position, $
     DATA_MULTIPLIER=data_multiplier, DATA_RANGE=data_range, $
     DEPTH_RANGE=depth_range, $
     PLOT_COLORS=plot_colors, $
     RECORD_AVERAGE=record_average, RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, RECORDS=records, $
     SHOW_TIME=show_time, TITLE=title, $
     USE_BATH=use_bath, VARIaBLE=variable, $
     GRAPH_PROPERTIES=graph_properties, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process file argument

   if n_elements(file) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'ofile'

   case size(file, /TNAME) of
      'STRING': begin
         is_station = strmatch(mgh_roms_file_type(file), 'ROMS*station file', /FOLD_CASE)
         self.file = is_station ? obj_new('MGHromsStation', file) : obj_new('MGHromsHistory', file)
      end
      'OBJREF': begin
         self.file = file
      end
      else: message, 'The argument is of the wrong data type'
   endcase
   ofile = self.file

   ;; Process variable argument

   if n_elements(variable) ne 1 then $
        message, 'The name of a variable must be supplied'

   if size(variable, /TNAME) ne 'STRING' then $
        message, 'The name of a variable must be supplied'

   self.variable = variable

   ;; Default data range depends on variable name

   mgh_roms_resolve_data, self.variable, $
        DATA_RANGE=data_range, DATA_MULTIPLIER=data_multiplier

   self.data_range = data_range

   if n_elements(title) eq 0 then title = ''
   if n_elements(show_time) eq 0 then show_time = 1B

   ;; Process plot colours

   if n_elements(plot_colors) eq 0 then $
        plot_colors = mgh_color('red')

   if (size(plot_colors, /DIMENSIONS))[0] ne 3 then $
        message, 'PLOT_COLORS must be dimensioned [3,n_colors]'

   n_colors = n_elements(plot_colors)/3

   ;; Handling of dimensions & variables depends on whether
   ;; the file is a history or stations file

   case !true of
      obj_isa(ofile, 'MGHromsHistory'): ftype = 'history'
      obj_isa(ofile, 'MGHromsStation'): ftype = 'station'
      else: message, 'Unknown file type'
   endcase

   ;; Does the file have (lon,lat) or (x,y) data? Remaining code is
   ;; not fully modified to use this yet.

   self.lonlat = ofile->HasVar('lon_rho') && ofile->HasVar('lat_rho')

   ;; Do we use dynamic bathymetry? Currently the bath variable is read
   ;; but assumed static.

   if n_elements(use_bath) eq 0 then use_bath = ofile->HasVar('bath')

   ;; Handling of dimensions & variables depends on whether
   ;; the file is a history or stations file

   case ftype of

      'history': begin

         ;; Check the variable dimensions are appropriate.

         dims = ofile->VarDims(self.variable)

         if min([strlen(dims.horizontal),strlen(dims.vertical)]) eq 0 then begin
            fmt = '(%"The variable %s must have horizontal & vertical dimensions")'
            message, string(FORMAT=fmt, self.variable)
         endif

         ;; Establish position & depth of point(s) to be plotted.

         if n_elements(position) eq 0 then begin
            dim_rho = ofile->DimRho()
            position = [dim_rho[0:1]-1]/2
         endif

         npdims = size(position, /N_DIMENSIONS)
         pdims = size(position, /DIMENSIONS)

         if npdims lt 1 || npdims gt 2 || pdims[0] ne 2 then $
              message, 'Position must be dimensioned [2] or [2,n_pos]'

         self.n_pos = npdims eq 2 ? pdims[1] : 1

         hpos = fltarr(self.n_pos)

         for p=0,self.n_pos-1 do begin

            if use_bath then begin
               hpos[p] = ofile->VarGet('bath', COUNT=[1,1,1], OFFSET=[position[*,p],0])
            endif else begin
               hpos[p] = ofile->VarGet('h', COUNT=[1,1], OFFSET=position[*,p])
            endelse

            if self.lonlat then begin
               xpos = ofile->VarGet('lon_rho', COUNT=[1,1], OFFSET=position[*,p])
               ypos = ofile->VarGet('lat_rho', COUNT=[1,1], OFFSET=position[*,p])
            endif else begin
               xpos = ofile->VarGet('x_rho', COUNT=[1,1], OFFSET=position[*,p])
               ypos = ofile->VarGet('y_rho', COUNT=[1,1], OFFSET=position[*,p])
            endelse

            if use_bath then begin
               fmt = '(%"Point %d, %d at %f, %f in an initial water depth of %f")'
            endif else begin
               fmt = '(%"Point %d, %d at %f, %f in a water depth of %f")'
            endelse
            print, FORMAT=fmt, position[*,p],xpos,ypos,hpos[p]

         endfor

      end

      'station': begin

         ;; Check the variable dimensions are appropriate for plotting.

         dims = ofile->VarDims(self.variable)

         if min([strlen(dims.station),strlen(dims.vertical)]) eq 0 then begin
            fmt = '(%"The variable %s must have vertical & station dimensions")'
            message, string(FORMAT=fmt, self.variable)
         endif

         ;; Establish position & depth of point to be plotted.

         if n_elements(position) eq 0 then position = 0

         self.n_pos = n_elements(position)

         hpos = fltarr(self.n_pos)

         for p=0,self.n_pos-1 do begin

            if use_bath then begin
               hpos[p] = ofile->VarGet('bath', COUNT=[1,1], OFFSET=position[p,0])
            endif else begin
               hpos[p] = ofile->VarGet('h', COUNT=[1], OFFSET=position[p])
            endelse

            ipos = ofile->VarGet('Ipos', COUNT=[1], OFFSET=position[p])
            jpos = ofile->VarGet('Jpos', COUNT=[1], OFFSET=position[p])

            if self.lonlat then begin
               xpos = ofile->VarGet('lon_rho', COUNT=[1], OFFSET=position[p])
               ypos = ofile->VarGet('lat_rho', COUNT=[1], OFFSET=position[p])
            endif else begin
               xpos = ofile->VarGet('x_rho', COUNT=[1], OFFSET=position[p])
               ypos = ofile->VarGet('y_rho', COUNT=[1], OFFSET=position[p])
            endelse
            if use_bath then begin
               fmt = '(%"Station %d position %d %d at %f, %f in an initial water depth of %f")'
            endif else begin
               fmt = '(%"Station %d position %d %d at %f, %f in a water depth of %f")'
            endelse
            print, FORMAT=fmt, position[p],ipos,jpos,xpos,ypos,hpos[p]

         endfor
      end

   endcase

   ;; Determine time variable. Should also check the variable's
   ;; "time" attribute

   if n_elements(record_average) eq 0 then record_average = 1

   has_time = strlen(dims.time) gt 0

   if has_time then begin
      n_time = ofile->DimInfo(dims.time, /DIMSIZE)
      mgh_resolve_indices, n_time, record_range, record_stride, records
      n_records = n_elements(records)
      time_var = ofile->TimeVarName(dims.time)
      if isa(time_var, /NULL) then message, 'Time variable not found'
   endif else begin
      n_records = 1
      if record_average gt 1 then $
           message, 'Cannot average over records for time-independent variable'
   endelse

   n_frames = long(n_records)/long(record_average)

   ;; Establish depth range

   if n_elements(depth_range) eq 0 then depth_range = [-0.02,1.02]*max(hpos)

   ;; Retrieve s-coordinate parameters and s values for the variable
   ;; The following code is complicated by a couple of changes in ROMS
   ;; output files:
   ;;  - In older files, the bottom value of s_w is omitted.
   ;;  - In ROMS 2.1 and earlier, the names of the s-coordinate variables did
   ;;    not match the corresponding dimensions--heavens knows why. In ROMS 2.2
   ;;    this was fixed.

   theta_s = ofile->VarGet('theta_s')
   theta_b = ofile->VarGet('theta_b')
   hc = ofile->VarGet('hc')
   vstretch = ofile->HasVar('Vstretching') ? ofile->VarGet('Vstretching') : 1
   vtransform = ofile->HasVar('Vtransform') ? ofile->VarGet('Vtransform') : 1

   if ofile->HasVar(dims.vertical) then begin
      svar = ofile->VarGet(dims.vertical)
   endif else begin
      case dims.vertical of
         's_rho': svar = ofile->VarGet('sc_r')
         's_w': begin
            svar = ofile->VarGet('sc_w')
            if n_elements(svar) eq ofile->DimInfo('s_rho', /DIMSIZE) then $
               svar = [0,svar]
         end
      endcase
   endelse

   n_s = n_elements(svar)

   ;; For now, calculate the vertical positions associated with
   ;; each profile once only, with the assumption that zeta = 0
   ;; and the initial bathymetry is preserved.

   pro_z = fltarr(n_s, self.n_pos)
   for p=0,self.n_pos-1 do begin
      pro_z[*,p] = mgh_roms_s_to_z(svar, hpos[p], ZETA=0, HC=hc, $
                                   THETA_S=theta_s, THETA_B=theta_b, $
                                   VSTRETCH=vstretch, VTRANSFORM=vtransform)
   endfor

   ;; Get time units from the input file

   if has_time then begin
      if ofile->HasAtt(time_var, 'units') then begin
         time_units = mgh_dt_units(ofile->AttGet(time_var, 'units'))
      endif else begin
         time_units = {scale: 1}
      endelse
   endif

   ;; Create base graph

   mgh_new, 'MGHgrGraph2D', NAME='ROMS profile animation', $
            _STRICT_EXTRA=graph_properties, RESULT=ograph

   ograph->NewMask

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Draw axes

   ograph->NewAxis, 0, RANGE=self.data_range, /EXACT, TITLE=self.variable, $
        _STRICT_EXTRA=xaxis_properties
   ograph->NewAxis, 1, RANGE=-reverse(depth_range), /EXACT, $
        TITLE='Depth (m)', TICKFORMAT='mgh_tf_negative', $
        _STRICT_EXTRA=yaxis_properties

   ;; Add a background to act as a selection target

   ograph->NewBackground

   ;; Add symbols, one per plot color

   for c=0,n_colors-1 do $
        ograph->NewSymbol, CLASS='mghgrsymbol', 0, COLOR=plot_colors[*,c]

   ;; Add various graphics objects that will (or might) be animated...

   ;; ...title

   otitle = ograph->NewTitle(title)

   ;; ...a series of plot objects

   oplot = objarr(self.n_pos)
   for p=0,self.n_pos-1 do begin
      col = p mod n_colors
      oplot[p] = ograph->NewAtom('IDLgrPlot', COLOR=plot_colors[*,col], $
                                 SYMBOL=ograph->GetSymbol(POSITION=col))
   endfor

   ;; Create an animator window to display and manage the movie.

   ok = self->MGH_Datamator::Init(CHANGEABLE=0, GRAPHICS_TREE=ograph, $
                                  MOUSE_ACTION=['Magnify','Translate', 'Context'], $
                                  _STRICT_EXTRA=extra)

   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Datamator'

   ;; Step through times in the netCDF file, generating new frames &
   ;; plotting data

   oframe = objarr(self.n_pos+1)

   ra = record_average
   dm = data_multiplier

   for f=0,n_frames-1 do begin

      if self->Finished() then break

      rec0 = ra*f

      time = 0  &  pro_data = dblarr(n_s, self.n_pos)

      for r=rec0,rec0+ra-1 do begin

         if has_time then begin
            t = ofile->VarGet(time_var, OFFSET=records[r], $
                             COUNT=[1], AUTOSCALE=0)*time_units.scale
            time += temporary(t)/double(ra)
         endif

         for p=0,self.n_pos-1 do begin

            ;; Get profile data

            case ftype of
               'history': begin
                  case strjoin(dims.horizontal, ' ') of
                     'xi_rho eta_rho' : begin
                        c = [1,1,0]
                        o = [position[*,p],0]
                        if has_time then begin
                           c = [c,1]
                           o = [o,records[r]]
                        endif
                        data = ofile->VarGet(self.variable, COUNT=c, OFFSET=o)
                        data = reform(data)
                     end
                     'xi_u eta_u' : begin
                        c = [2,1,0]
                        o = [position[*,p],0]
                        if has_time then begin
                           c = [c,1]
                           o = [o,records[r]]
                        endif
                        data = ofile->VarGet(self.variable, COUNT=c, OFFSET=o)
                        data = reform(mgh_avg(data,1))
                     end
                     'xi_v eta_v' : begin
                        c = [1,2,0]
                        o = [position[*,p],0]
                        if has_time then begin
                           c = [c,1]
                           o = [o,records[r]]
                        endif
                        data = ofile->VarGet(self.variable, COUNT=c, OFFSET=o)
                        data = reform(mgh_avg(data,2))
                     end
                  endcase
               end
               'station': begin
                  c = [0,1]
                  o = [0,position[p]]
                  if has_time then begin
                     c = [c,1]
                     o = [o,records[r]]
                  endif
                  data = ofile->VarGet(self.variable, COUNT=c, OFFSET=o)
                  data = reform(data)
               end
            endcase

	         pro_data[*,p] += dm*temporary(data)/double(ra)

	      endfor

      endfor

      ;; Display title

      if has_time && show_time gt 0 then begin
         case show_time of
            1: ttt = string(FORMAT='(%"%0.3f days")', time)
            2: ttt = mgh_dt_string(time+time_units.offset)
         endcase
         if strlen(title) gt 0 then $
              ttt = string(FORMAT='(%"%s: %s")', title, ttt)
         oframe[self.n_pos] = obj_new('MGH_Command', OBJECT=otitle, 'SetProperty', STRINGS=temporary(ttt))
      endif

      ;; Display profile

      for p=0,self.n_pos-1 do $
           oframe[p] = obj_new('MGH_Command', OBJECT=oplot[p], 'SetProperty', DATAX=pro_data[*,p], DATAY=pro_z[*,p])

      ;; Add frame to animator & display

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; Mgh_Roms_Movie_Profile::GetProperty
;
PRO Mgh_Roms_Movie_Profile::GetProperty, $
     DATA_RANGE=data_range, FILE=file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   data_range = self.data_range

   file = self.file

   self->MGH_Datamator::GetProperty, _STRICT_EXTRA=extra

END

; Mgh_Roms_Movie_Profile::SetProperty
;
PRO Mgh_Roms_Movie_Profile::SetProperty, $
     DATA_RANGE=data_range, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->GetProperty, GRAPHICS_TREE=graph

   if n_elements(data_range) gt 0 then begin
      self.data_range = data_range
      if obj_valid(graph) then begin
         xaxis = graph->GetAxis(DIRECTION=0)
         xaxis->SetProperty, RANGE=self.data_range
      endif
   endif

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=extra

end

; Mgh_Roms_Movie_Profile::About
;
pro Mgh_Roms_Movie_Profile::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::About, lun

   if obj_valid(self.file) then begin
      printf, lun, self, ': the history/station file object is ', self.file
      self.file->GetProperty, FILE_NAME=file_name
      printf, lun, self, ': the files are:', file_name
   endif

end

; Mgh_Roms_Movie_Profile::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Movie_Profile::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::BuildMenuBar

   ombar = mgh_widget_self(self.menu_bar)

   if obj_valid(ombar) then begin
      ombar->NewItem, PARENT='Tools', SEPARATOR=[1,1], $
        ['Set Data Range...','View Profile...']
   endif

end

; Mgh_Roms_Movie_Profile::EventMenubar
;
function Mgh_Roms_Movie_Profile::EventMenubar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      'TOOLS.SET DATA RANGE': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Range', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, /IMMEDIATE, $
                  N_ELEMENTS=2, PROPERTY_NAME='DATA_RANGE'
         return, 0
      end

      'TOOLS.VIEW PROFILE': begin
         self->GetProperty, POSITION=position
         oframe = self.animation->GetFrame(POSITION=position)
         oframe[0]->GetProperty, KEYWORDS=keywords
         data_dims = size(keywords.datax, /DIMENSIONS)
         xvaredit, keywords.datax, GROUP=self.base, $
                   X_SCROLL_SIZE=(data_dims[0] < 12), $
                   Y_SCROLL_SIZE=(data_dims[1] < 30)
         return, 0
      end

      else: return, self->MGH_Datamator::EventMenubar(event)

   endcase

end

; Mgh_Roms_Movie_Profile::ExportData
;
pro Mgh_Roms_Movie_Profile::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::ExportData, values, labels

   self->GetProperty, ANIMATION=animation, FILE=file, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[0]->GetProperty, KEYWORDS=keywords

   labels = [labels, 'File Object', 'Profile Depth', 'Profile Data']
   values = [values, ptr_new(file), ptr_new(keywords.datay),ptr_new(keywords.datax)]

end

pro Mgh_Roms_Movie_Profile__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {Mgh_Roms_Movie_Profile, inherits MGH_Datamator, $
         file: obj_new(), variable: '', n_pos: 0UL, $
         lonlat: !false, data_range: fltarr(2)}

end
