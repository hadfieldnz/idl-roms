;+
; CLASS NAME:
;   Mgh_Roms_Movie_Vector_Profile
;
; PURPOSE:
;   This class implements an animated sequence of graphs showing
;   current vector profiles from a SCRUM/ROMS history OR station
;   file
;
; OBJECT CREATION SEQUENCE
;    mgh_new, 'mgh_roms_movie_vector_profile', history, position
;
; POSITIONAL PARAMETERS:
;   file (input, object reference or string)
;     The name of a ROMS history file OR a reference to an
;     MGHromsHistory or MGHromsStation object.
;
;   position (input, numeric)
;     For a history file, a 2-element vector specifying the
;     [xi,eta] location of the profile; for a station file, a scalar
;     integer specifying the station number.
;
; KEYWORD PARAMETERS:
;   DEPTH_RANGE
;     If supplied, this keyword specifies the height range of the plot
;
;   SPEED_MAX
;     This keyword specifies the range of the speed axes. Default is
;     1.
;
;   VARNAME
;     A 2-element string array specifying the names of a pair of u & v
;     type variables. Default is ['u','v'].
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2001-05:
;       Written.
;   Mark Hadfield, 2009-02:
;       Generalised history-file code to handle fraction positions.
;-

function Mgh_Roms_Movie_Vector_Profile::Init, file, position, $
     DATA_RANGE=data_range, $
     DEPTH_RANGE=depth_range, $
     GRAPH_PROPERTIES=graph_properties, $
     SPEED_MAX=speed_max, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     RECORDS=records, $
     VARNAME=varname, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process file argument

   case size(file, /TNAME) of
      'STRING': begin
         case 1B of
            strmatch(mgh_roms_file_type(file), '*station file', /FOLD_CASE): $
                 self.file = obj_new('MGHromsStation', file)
            else: $
                 self.file = obj_new('MGHromsHistory', file)
         endcase
         self.fdestroy = 1B
      end
      'OBJREF': begin
         self.file = file
         self.fdestroy = 0B
      end
      else: message, 'The argument is of the wrong data type'
   endcase
   ofile = self.file

   ;; Process arguments

   self.varname = n_elements(varname) gt 0 ? varname : ['u','v']

   self.speed_max = n_elements(speed_max) gt 0 ? speed_max : 0.2

   if n_elements(visible) eq 0 then visible = 1

   ;; Establish file type.

   case 1B of
      obj_isa(ofile, 'MGHromsHistory'): ftype = 'history'
      obj_isa(ofile, 'MGHromsStation'): ftype = 'station'
      else: message, 'Unknown file type'
   endcase

   ;; Handling of dimensions & variables depends on whether the file
   ;; is a history or stations file

   case ftype of

      'history': begin

         ;; Read grid dimensions

         dim_rho = ofile->DimRho()

         ;; Check the variable dimensions are appropriate

         dim_u = ofile->VarDims(self.varname[0])
         dim_v = ofile->VarDims(self.varname[1])

         fmt = '(%"Variable %s must have horizontal, vertical & time dimensions")'
         if min([strlen(dim_u.horizontal),strlen(dim_u.vertical), $
                 strlen(dim_u.time)]) eq 0 then begin
            message, string(FORMAT=fmt, self.varname[0])
         endif
         if min([strlen(dim_v.horizontal),strlen(dim_v.vertical), $
                 strlen(dim_v.time)]) eq 0 then begin
            message, string(FORMAT=fmt, self.varname[1])
         endif

         if ~ strmatch(dim_u.time, dim_v.time) then $
              message, 'Dimension mismatch'
         if ~ strmatch(dim_u.vertical, dim_v.vertical) then $
              message, 'Dimension mismatch'

         ;; Establish position of point to be plotted and calculate depth,
         ;; grid angle, location, etc.

         if n_elements(position) eq 0 then position = 0.5*(dim_rho[0:1]-1)

         ipos = floor(position)
         rpos = position - ipos

         hvar = interpolate(ofile->VarGet('h', COUNT=[2,2], OFFSET=ipos), $
                            rpos[0], rpos[1])

         avar = 0
         if ofile->HasVar('angle') then $
              avar = interpolate(ofile->VarGet('angle', COUNT=[2,2], OFFSET=ipos), $
                                 rpos[0], rpos[1])

         self.lonlat = ofile->HasVar('lon_rho') && ofile->HasVar('lat_rho')

         case self.lonlat of
            0B: begin
               self.x = interpolate(ofile->VarGet('x_rho', COUNT=[2,2], OFFSET=ipos), $
                                  rpos[0], rpos[1])
               self.y = interpolate(ofile->VarGet('y_rho', COUNT=[2,2], OFFSET=ipos), $
                                  rpos[0], rpos[1])
            end
            1B: begin
               self.x = interpolate(ofile->VarGet('lon_rho', COUNT=[2,2], OFFSET=ipos), $
                                  rpos[0], rpos[1])
               self.y = interpolate(ofile->VarGet('lat_rho', COUNT=[2,2], OFFSET=ipos), $
                                  rpos[0], rpos[1])
            end
         endcase

         fmt = '(%"Point %s, %s at %s, %s, angle %s deg, depth %s m")'
         message, /INFORM, $
                  string(FORMAT=fmt, mgh_format_float([position,self.x,self.y,avar*!radeg, $
                                                       hvar]))

      end

      'station': begin

         ;; Check the variable dimensions are appropriate

         dim_u = ofile->VarDims(self.varname[0])
         dim_v = ofile->VarDims(self.varname[1])

         fmt = '(%"Variable %s must have vertical, station & time dimensions")'
         if min([strlen(dim_u.station),strlen(dim_u.vertical), $
                 strlen(dim_u.time)]) eq 0 then begin
            message, string(FORMAT=fmt, self.varname[0])
         endif
         if min([strlen(dim_v.station),strlen(dim_v.vertical), $
                 strlen(dim_v.time)]) eq 0 then begin
            message, string(FORMAT=fmt, self.varname[1])
         endif

         if ~ strmatch(dim_u.time, dim_v.time) then $
              message, 'Dimension mismatch'
          if ~ strmatch(dim_u.station, dim_v.station) then $
              message, 'Dimension mismatch'
         if ~ strmatch(dim_u.vertical, dim_v.vertical) then $
              message, 'Dimension mismatch'

         ;; Establish position & depth of point to be plotted.

         if n_elements(position) eq 0 then position = 0

         hvar = ofile->VarGet('h', COUNT=[1], OFFSET=position)

         avar = 0
         if ofile->HasVar('angle') then $
              avar = ofile->VarGet('angle', COUNT=[1], OFFSET=position)

         ivar = ofile->VarGet('Ipos', COUNT=[1], OFFSET=position)
         jvar = ofile->VarGet('Jpos', COUNT=[1], OFFSET=position)

         self.lonlat = ofile->HasVar('lon_rho') && ofile->HasVar('lat_rho')

         case self.lonlat of
            0B: begin
               self.x = ofile->VarGet('x_rho', COUNT=[1], OFFSET=position)
               self.y = ofile->VarGet('y_rho', COUNT=[1], OFFSET=position)
            end
            1B: begin
               self.x = ofile->VarGet('lon_rho', COUNT=[1], OFFSET=position)
               self.y = ofile->VarGet('lat_rho', COUNT=[1], OFFSET=position)
            end
         endcase


         fmt = '(%"Station %d, position %s, %s at %s, %s, angle %s deg, depth %s m")'
         message, /INFORM, $
                  string(FORMAT=fmt, position, mgh_format_float([ivar,jvar,self.x,self.y, $
                                                                 avar*!radeg,hvar]))

      end

   endcase

   ;; Set the plot's depth range

   if n_elements(depth_range) eq 0 then depth_range = [-0.02*hvar,1.02*hvar]

   ;; Establish plotting time

   ;; Should also check the variable's "time" attribute
   case 1B of
      ofile->HasVar(dim_u.time): $
           time_var = dim_u.time
      ofile->HasVar('ocean_time'): $
           time_var = 'ocean_time'
      ofile->HasVar('scrum_time'): $
           time_var = 'scrum_time'
      else: $
           message, 'Time variable not found'
   endcase

   n_time = ofile->DimInfo(dim_u.time, /DIMSIZE)
   mgh_resolve_indices, n_time, record_range, record_stride, records
   n_records = n_elements(records)

   self.time = ptr_new(dblarr(n_records))

   ;; Retrieve depth information for the variable

   theta_s = ofile->VarGet('theta_s')
   theta_b = ofile->VarGet('theta_b')
   hc = ofile->VarGet('hc')
   vstretch = ofile->HasVar('Vstretching') ? ofile->VarGet('Vstretching') : 1
   vtransform = ofile->HasVar('Vtransform') ? ofile->VarGet('Vtransform') : 1

   case dim_u.vertical of
      's_rho': svar = ofile->VarGet('s_rho')
      's_w'  : svar = ofile->VarGet('s_w')
   endcase

   n_s = n_elements(svar)

   ;; Calculate depth with the assumption that zeta = 0.

   hh = mgh_roms_s_to_z(svar, hvar, ZETA=0, $
                        HC=hc, THETA_S=theta_s, THETA_B=theta_b, $
                        VSTRETCH=vstretch, VTRANSFORM=vtransform)

   ;; For history files, set up arrays for horizontal interpolation of (u,v)
   ;; data

   if ftype eq 'history'then begin
      ii = replicate(rpos[0], n_s)
      jj = replicate(rpos[1], n_s)
      kk = findgen(n_s)
   endif

   ;; Create base graph

   ograph = obj_new('MGHgrGraph3D', $
                    NAME='ROMS velocity profile animation', $
                    COLOR=replicate(225B,3), $
                    _STRICT_EXTRA=graph_properties)

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Draw axes

   ograph->NewAxis, 0, RANGE=[-1,1]*self.speed_max, /EXACT, $
        TITLE='X', _STRICT_EXTRA=xaxis_properties
   ograph->NewAxis, 1, RANGE=[-1,1]*self.speed_max, /EXACT, $
        TITLE='Y', _STRICT_EXTRA=yaxis_properties
   ograph->NewAxis, 2, RANGE=-reverse(depth_range), /EXACT, $
        TITLE='Depth (m)', TICKFORMAT='mgh_tf_negative', $
        _STRICT_EXTRA=zaxis_properties

   ;; Barb plot object & title to be animated

   ograph->NewAtom, 'MGHgrBarb', DATAZ=hh, COLOR=mgh_color('red'), RESULT=obarb
   self.barb = obarb

   otitle = ograph->NewTitle()

   ;; Create an animator window to display and manage the movie.

   ok = self->MGH_Datamator::Init(CHANGEABLE=0, GRAPHICS_TREE=ograph, $
                                  MOUSE_ACTION=['Rotate','Pick','Context'], $
                                  _STRICT_EXTRA=extra)
   if ~ ok then message, 'MGH_Datamator initialisation failed'

   ;; Step through times in the netCDF file, generating new frames &
   ;; plotting data

   oframe = objarr(2)

   for r=0,n_records-1 do begin

      ;; Check to see if the user has selected the "Finish Loading" menu item

      if self->Finished() then break

      ;; Get profile

      f = records[r]

      case ftype of
         'history': begin
            uv = ofile->VectorGet(self.varname, COUNT=[2,2,0,1], OFFSET=[ipos,0,f])
            uv = interpolate(temporary(uv), ii, jj, kk)
         end
         'station': begin
            uu = ofile->VarGet(self.varname[0], COUNT=[0,1,1], OFFSET=[0,position,f])
            vv = ofile->VarGet(self.varname[1], COUNT=[0,1,1], OFFSET=[0,position,f])
            uv = complex(temporary(uu), temporary(vv))
         end
      endcase

      uv *= exp(complex(0, 1)*avar)

      oframe[0] = obj_new('MGH_Command', OBJECT=self.barb, $
                          'SetProperty', DATAU=real_part(uv), DATAV=imaginary(uv))

      mgh_undefine, uv

      t = ofile->VarGet(time_var, OFFSET=[f], COUNT=[1], /AUTOSCALE)
      oframe[1] = obj_new('MGH_Command', OBJECT=otitle,'SetProperty', $
                          STRINGS=string(t, FORMAT='(%"%0.3f days")'))
      (*self.time)[r] = t

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; Mgh_Roms_Plot_Hstats::Cleanup
;
pro MGH_ROMS_Movie_Vector_Profile::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if self.fdestroy then obj_destroy, self.file

   ptr_free, self.time

   self->MGH_Datamator::Cleanup

end

; Mgh_Roms_Movie_Vector_Profile::GetProperty
;
pro MGH_ROMS_Movie_Vector_Profile::GetProperty, $
     FILE=file, LONLAT=lonlat, $
     SPEED_MAX=speed_max, $
     VARNAME=varname, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   file = self.file

   lonlat = self.lonlat

   speed_max = self.speed_max

   varname = self.varname

   self->MGH_Datamator::GetProperty, _STRICT_EXTRA=extra

END

; Mgh_Roms_Movie_Vector_Profile::SetProperty
;
pro Mgh_Roms_Movie_Vector_Profile::SetProperty, SPEED_MAX=speed_max, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->GetProperty, GRAPHICS_TREE=graph

   if n_elements(speed_max) gt 0 then begin
      self.speed_max = speed_max
      if obj_valid(graph) then begin
         xaxis = graph->GetAxis(DIRECTION=0)
         xaxis->SetProperty, RANGE=self.speed_max*[-1,1]
         yaxis = graph->GetAxis(DIRECTION=1)
         yaxis->SetProperty, RANGE=self.speed_max*[-1,1]
         graph->GetScaling, XCOORD_CONV=xcoord, YCOORD_CONV=ycoord
         self.barb->SetProperty, SCALE=[xcoord[1],ycoord[1],1]
      endif
   endif

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=extra

end

; Mgh_Roms_Movie_Vector_Profile::About
;
pro Mgh_Roms_Movie_Vector_Profile::About, lun

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

; Mgh_Roms_Movie_Vector_Profile::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Movie_Vector_Profile::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::BuildMenuBar

   obar = mgh_widget_self(self.menu_bar)

   if obj_valid(obar) then begin

      obar->NewItem, PARENT='File.Export Animation', 'NetCDF...'

      obar->NewItem, PARENT='File.Export Frame', 'NetCDF...'

      obar->NewItem, PARENT='Tools', ['Set Max Speed...'], SEPARATOR=[1]

   endif

end


; Mgh_Roms_Movie_Vector_Profile::EventMenuBar
;
function Mgh_Roms_Movie_Vector_Profile::EventMenuBar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      'FILE.EXPORT ANIMATION.NETCDF': begin
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
            self->ExportToNcFile, filename, /ALL
         endif
         return, 0
      end

      'FILE.EXPORT FRAME.NETCDF': begin
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

      'TOOLS.SET MAX SPEED': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Max speed', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, /IMMEDIATE, $
                  N_ELEMENTS=1, PROPERTY_NAME='SPEED_MAX'
         return, 1
      end

      else: return, self->MGH_Datamator::EventMenuBar(event)

   endcase

end

; Mgh_Roms_Movie_Vector_Profile::ExportData
;
pro Mgh_Roms_Movie_Vector_Profile::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::ExportData, values, labels

   self->GetProperty, ANIMATION=animation, FILE=file, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[0]->GetProperty, KEYWORDS=keywords

   labels = [labels, 'File Object', 'U,V Profile']
   values = [values, ptr_new(file), $
             ptr_new(complex(keywords.datau,keywords.datav))]

end

; MGH_Roms_Movie_Vector_Profile::ExportToNcFile
;
pro MGH_Roms_Movie_Vector_Profile::ExportToNcFile, file, position, ALL=all

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->GetProperty, $
        ANIMATION=animation, LONLAT=lonlat, VARNAME=varname

   if n_elements(position) eq 0 then begin
      case keyword_set(all) of
         0B: begin
            self->GetProperty, POSITION=position
         end
         1B: begin
            animation->GetProperty, N_FRAMES=n_frames
            position = lindgen(n_frames)
         end
      endcase
   endif

   self.barb->GetProperty, DATAZ=dataz

   fmt = '(%"Saving %d records velocity data to %s")'
   message, /INFORM, string(FORMaT=fmt, n_elements(position), file)

   onc = obj_new('MGHncFile', file, /CREATE, /CLOBBER)

   onc->AttAdd, /GLOBAL, 'title', 'ROMS velocity profile data'

   onc->AttAdd, /GLOBAL, 'history', $
                'Generated by routine MGH_ROMS_Movie_Velocity_Profile::ExportToNcFile at '+ $
                mgh_dt_string(mgh_dt_now())

   case self.lonlat of
      0B: begin
         onc->AttAdd, /GLOBAL, 'x', self.x
         onc->AttAdd, /GLOBAL, 'y', self.y
      end
      1B: begin
         onc->AttAdd, /GLOBAL, 'lon', self.x
         onc->AttAdd, /GLOBAL, 'lat', self.y
      end
   endcase

   onc->DimAdd, 's', n_elements(dataz)
   onc->DimAdd, 'time'

   onc->VarAdd, 'z', ['s']

   onc->VarAdd, 'time', ['time'], /DOUBLE
   onc->AttAdd, 'time', 'units', 'days'

   onc->VarAdd, self.varname[0], ['s','time']
   onc->VarAdd, self.varname[1], ['s','time']

   onc->VarPut, 'z', dataz

   for p=0,n_elements(position)-1 do begin

      onc->VarPut, 'time', (*self.time)[position[p]], OFFSET=[p]

      oframe = animation->GetFrame(POSITION=position[p])
      oframe[0]->GetProperty, KEYWORDS=keywords

      onc->VarPut, self.varname[0], keywords.datau, OFFSET=[0,p]
      onc->VarPut, self.varname[1], keywords.datav, OFFSET=[0,p]

   endfor

   obj_destroy, onc


end

pro MGH_Roms_Movie_Vector_Profile__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {Mgh_Roms_Movie_Vector_Profile, inherits MGH_Datamator, $
         file: obj_new(), fdestroy: 0B, lonlat: 0B, $
         x: 0.D0, y: 0.D0, varname: strarr(2), $
         barb: obj_new(), speed_max: 0.D0, time: ptr_new()}

end
