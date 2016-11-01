;+
; CLASS NAME:
;   Mgh_Roms_Movie_TS
;
; PURPOSE:
;   This class implements an animated sequence of graphs showing
;   temperature-salinity (or more general property-property) plots from a ROMS
;   history OR station file
;
; OBJECT CREATION SEQUENCE
;   mgh_new, 'Mgh_Roms_Movie_TS', File, Var, Position
;
; POSITIONAL PARAMETERS:
;   file
;     The name of a ROMS history file OR a reference to an
;     MGHromsHistory or MGHromsStation object.
;
;   position
;     For a history file, a 2-element integer vector specifying the
;     [xi,eta] location of the profile; for a station file, a scalar
;     integer specifying the station number.
;
; KEYWORD PARAMETERS:
;   DATA_RANGE_VAR0, DATA_RANGE_VAR1
;     If supplied, these keywords specify the x & y-axis range.
;
;   VARNAME
;     A 2-element string array specifying the names of a pair of
;     rho-type variables. Default is ['salt','temp'].
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2009-11:
;     Written.
;-

function Mgh_Roms_Movie_TS::Init, file, position, $
     DATA_RANGE_VAR0=data_range_var0, DATA_RANGE_VAR1=data_range_var1, $
     GRAPH_PROPERTIES=graph_properties, PLOT_COLORS=plot_colors, $
     RECORD_RANGE=record_range, RECORD_STRIDE=record_stride, RECORDS=records, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process file argument

   case size(file, /TNAME) of
      'STRING': begin
         case 1B of
            strmatch(mgh_roms_file_type(file), 'ROMS*station file', /FOLD_CASE): $
                 self.file = obj_new('MGHromsStation', file)
            else: $
                 self.file = obj_new('MGHromsHistory', file)
         endcase
         self.file_destroy = 1
      end
      'OBJREF': begin
         self.file = file
         self.file_destroy = 0
      end
      else: message, 'The argument is of the wrong data type'
   endcase
   ofile = self.file

   ;; Process varname argument

   if n_elements(varname) eq 0 then $
        varname = ['salt','temp']

   self.varname = varname

   ;; Default data range depends on variable name

   mgh_roms_resolve_data, self.varname[0], DATA_RANGE=data_range_var0
   mgh_roms_resolve_data, self.varname[1], DATA_RANGE=data_range_var1

   self.data_range_var0 = data_range_var0
   self.data_range_var1 = data_range_var1

   ;; Process plot colours

   if n_elements(plot_colors) eq 0 then $
        plot_colors = mgh_color('red')

   if (size(plot_colors, /DIMENSIONS))[0] ne 3 then $
        message, 'PLOT_COLORS must be dimensioned [3,n_colors]'

   n_colors = n_elements(plot_colors)/3

   ;; Handling of dimensions & variables depends on whether
   ;; the file is a history or stations file

   case 1 of
      obj_isa(ofile, 'MGHromsHistory'): ftype = 'history'
      obj_isa(ofile, 'MGHromsStation'): ftype = 'station'
      else: message, 'Unknown file type'
   endcase

   ;; Does the file have (lon,lat) or (x,y) data? Remaining code is
   ;; not fully modified to use this yet.

   self.lonlat = ofile->HasVar('lon_rho') && ofile->HasVar('lat_rho')

   ;; Handling of dimensions & variables depends on whether
   ;; the file is a history or stations file

   case ftype of

      'history': begin

         ;; Check the variable dimensions.

         dim0 = ofile->VarDims(self.varname[0])
         dim1 = ofile->VarDims(self.varname[1])

         fmt = '(%"Variable %s must have horizontal, vertical & time dimensions")'
         if min([strlen(dim0.horizontal),strlen(dim0.vertical), $
                 strlen(dim0.time)]) eq 0 then begin
            message, string(FORMAT=fmt, self.varname[0])
         endif
         if min([strlen(dim1.horizontal),strlen(dim1.vertical), $
                 strlen(dim1.time)]) eq 0 then begin
            message, string(FORMAT=fmt, self.varname[1])
         endif

         if ~ strmatch(dim0.time, dim1.time) then $
              message, 'Dimension mismatch'
         if ~ array_equal(dim0.horizontal, dim1.horizontal) then $
              message, 'Dimension mismatch'
         if ~ strmatch(dim0.vertical, dim1.vertical) then $
              message, 'Dimension mismatch'

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

            hpos[p] = ofile->VarGet('h', COUNT=[1,1], OFFSET=position[*,p])

            case self.lonlat of

               0: begin
                  xpos = ofile->VarGet('x_rho', COUNT=[1,1], OFFSET=position[*,p])
                  ypos = ofile->VarGet('y_rho', COUNT=[1,1], OFFSET=position[*,p])
               end

               1: begin
                  xpos = ofile->VarGet('lon_rho', COUNT=[1,1], OFFSET=position[*,p])
                  ypos = ofile->VarGet('lat_rho', COUNT=[1,1], OFFSET=position[*,p])
               end

            endcase

            fmt = '(%"Point %d, %d at %f, %f in a water depth of %f")'
            print, FORMAT=fmt, position[*,p],xpos,ypos,hpos[p]

         endfor

      end

      'station': begin

         ;; Check the variable dimensions.

         dim0 = ofile->VarDims(self.varname[0])
         dim1 = ofile->VarDims(self.varname[1])

         fmt = '(%"Variable %s must have vertical, station & time dimensions")'
         if min([strlen(dim0.station),strlen(dim0.vertical), $
                 strlen(dim0.time)]) eq 0 then begin
            message, string(FORMAT=fmt, self.varname[0])
         endif
         if min([strlen(dim1.station),strlen(dim1.vertical), $
                 strlen(dim1.time)]) eq 0 then begin
            message, string(FORMAT=fmt, self.varname[1])
         endif

         if ~ strmatch(dim0.time, dim1.time) then $
              message, 'Dimension mismatch'
         if ~ strmatch(dim0.station, dim1.station) then $
              message, 'Dimension mismatch'
         if ~ strmatch(dim0.vertical, dim1.vertical) then $
              message, 'Dimension mismatch'

         ;; Establish position & depth of point to be plotted.

         if n_elements(position) eq 0 then position = 0

         self.n_pos = n_elements(position)

         hpos = fltarr(self.n_pos)

         for p=0,self.n_pos-1 do begin

            hpos[p] = ofile->VarGet('h', COUNT=[1], OFFSET=position[p])

            ipos = ofile->VarGet('Ipos', COUNT=[1], OFFSET=position[p])
            jpos = ofile->VarGet('Jpos', COUNT=[1], OFFSET=position[p])

            case self.lonlat of

               0: begin
                  xpos = ofile->VarGet('x_rho', COUNT=[1], OFFSET=position[p])
                  ypos = ofile->VarGet('y_rho', COUNT=[1], OFFSET=position[p])

                  fmt = '(%"Station %d position %d, %d at %f, %f in a water depth of %f")'
                  print, FORMAT=fmt, position[p],ipos,jpos,xpos,ypos,hpos[p]
               end

               1: begin
                  xpos = ofile->VarGet('lon_rho', COUNT=[1], OFFSET=position[p])
                  ypos = ofile->VarGet('lat_rho', COUNT=[1], OFFSET=position[p])

                  fmt = '(%"Station %d position %d %d at lon %f, lat %f in ' + $
                        'a water depth of %f")'
                  print, FORMAT=fmt, position[p],ipos,jpos,xpos,ypos,hpos[p]
               end

            endcase

         endfor
      end

   endcase

   ;; Should also check the variables "time" attribute
   time_var = ohis->TimeVarName(dim0.time)
   if isa(time_var, /NULL) then message, 'Time variable not found'

   n_time = ofile->DimInfo(dim0.time, /DIMSIZE)
   mgh_resolve_indices, n_time, record_range, record_stride, records
   n_records = n_elements(records)

   ;; Create base graph

   mgh_new, 'MGHgrGraph2D', NAME='ROMS profile animation', $
            _STRICT_EXTRA=graph_properties, RESULT=ograph

   ograph->NewMask

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Draw axes

   ograph->NewAxis, 0, RANGE=self.data_range_var0, /EXACT, TITLE=self.varname[0]
   ograph->NewAxis, 1, RANGE=self.data_range_var1, /EXACT, TITLE=self.varname[1]

   ;; Add a background to act as a selection target

   ograph->NewBackground

   ;; Add symbols, one per plot color

   for c=0,n_colors-1 do $
         ograph->NewSymbol, CLASS='mghgrsymbol', 0, COLOR=plot_colors[*,c]

   ;; Add a text object and a set of plot objects to be animated

   otitle = ograph->NewTitle()

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

   for r=0,n_records-1 do begin

      f = records[r]

      if self->Finished() then break

      for p=0,self.n_pos-1 do begin

         ;; Get profile

         case ftype of
            'history': begin
               case strjoin(dim0.horizontal, ' ') of
                  'xi_rho eta_rho' : begin
                     data0 = ofile->VarGet(self.varname[0], COUNT=[1,1,0,1], $
                                           OFFSET=[position[*,p],0,f])
                     data1 = ofile->VarGet(self.varname[1], COUNT=[1,1,0,1], $
                                           OFFSET=[position[*,p],0,f])
                     data0 = reform(data0)
                     data1 = reform(data1)
                  end
                  'xi_u eta_u' : begin
                     data0 = ofile->VarGet(self.varname[0], COUNT=[2,1,0,1], $
                                           OFFSET=[position[*,p],0,f])
                     data1 = ofile->VarGet(self.varname[1], COUNT=[2,1,0,1], $
                                           OFFSET=[position[*,p],0,f])
                     data0 = reform(mgh_avg(data0, 1))
                     data1 = reform(mgh_avg(data1, 1))
                  end
                  'xi_v eta_v' : begin
                     data0 = ofile->VarGet(self.varname[0], COUNT=[1,2,0,1], $
                                           OFFSET=[position[*,p],0,f])
                     data1 = ofile->VarGet(self.varname[1], COUNT=[1,2,0,1], $
                                           OFFSET=[position[*,p],0,f])
                     data0 = reform(mgh_avg(data0, 2))
                     data1 = reform(mgh_avg(data1, 2))
                  end
               endcase
            end
            'station': begin
               data0 = ofile->VarGet(self.varname[0], COUNT=[0,1,1], $
                                       OFFSET=[0,position[p],f])
               data1 = ofile->VarGet(self.varname[1], COUNT=[0,1,1], $
                                       OFFSET=[0,position[p],f])
               data0 = reform(data0)
               data1 = reform(data1)
            end
         endcase

         ;; Command to set plot properties

         oframe[p] = obj_new('MGH_Command', OBJECT=oplot[p], 'SetProperty', $
                             DATAX=data0, DATAY=data1)

      endfor

      ;; Title

      if ofile->HasVar(time_var) then begin
         t = ofile->VarGet(time_var, OFFSET=[f], COUNT=[1], /AUTOSCALE)
         oframe[self.n_pos] = $
              obj_new('MGH_Command', OBJECT=otitle, 'SetProperty', $
                      STRINGS=string(t,FORMAT='(F0.3)')+' days')
      endif

      ;; Add frame to animator & display

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; Mgh_Roms_Movie_TS::GetProperty
;
PRO Mgh_Roms_Movie_TS::GetProperty, $
     DATA_RANGE_VAR0=data_range_var0, DATA_RANGE_VAR1=data_range_var1, $
     FILE=file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   data_range_var0 = self.data_range_var0
   data_range_var1 = self.data_range_var1

   file = self.file

   self->MGH_Datamator::GetProperty, _STRICT_EXTRA=extra

END

; Mgh_Roms_Movie_TS::SetProperty
;
PRO Mgh_Roms_Movie_TS::SetProperty, $
     DATA_RANGE_VAR0=data_range_var0, $
     DATA_RANGE_VAR1=data_range_var1, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->GetProperty, GRAPHICS_TREE=graph

   if n_elements(data_range_var0) gt 0 then begin
      self.data_range_var0 = data_range_var0
      if obj_valid(graph) then begin
         xaxis = graph->GetAxis(DIRECTION=0)
         xaxis->SetProperty, RANGE=self.data_range_var0
      endif
   endif

   if n_elements(data_range_var1) gt 0 then begin
      self.data_range_var1 = data_range_var1
      if obj_valid(graph) then begin
         yaxis = graph->GetAxis(DIRECTION=1)
         yaxis->SetProperty, RANGE=self.data_range_var1
      endif
   endif

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=extra

end

; Mgh_Roms_Movie_TS::About
;
pro Mgh_Roms_Movie_TS::About, lun

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

; Mgh_Roms_Movie_TS::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Movie_TS::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::BuildMenuBar

   ombar = mgh_widget_self(self.menu_bar)

   if obj_valid(ombar) then begin
      ombar->NewItem, PARENT='Tools', SEPARATOR=[1,0], $
        ['Set X Range...','Set Y Range...']
   endif

end

; Mgh_Roms_Movie_TS::EventMenubar
;
function Mgh_Roms_Movie_TS::EventMenubar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      'TOOLS.SET X RANGE': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Range', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, /IMMEDIATE, $
                  N_ELEMENTS=2, PROPERTY_NAME='DATA_RANGE_VAR0'
         return, 0
      end

      'TOOLS.SET Y RANGE': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Range', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, /IMMEDIATE, $
                  N_ELEMENTS=2, PROPERTY_NAME='DATA_RANGE_VAR1'
         return, 0
      end

      else: return, self->MGH_Datamator::EventMenubar(event)

   endcase

end

; Mgh_Roms_Movie_TS::ExportData
;
pro Mgh_Roms_Movie_TS::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::ExportData, values, labels

   self->GetProperty, ANIMATION=animation, FILE=file, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[0]->GetProperty, KEYWORDS=keywords

   labels = [labels, 'File Object', 'Data X', 'Data Y']
   values = [values, ptr_new(file), ptr_new(keywords.datax), $
             ptr_new(keywords.datay)]

end

pro Mgh_Roms_Movie_TS__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {Mgh_Roms_Movie_TS, inherits MGH_Datamator, $
         file: obj_new(), file_destroy: 0B, n_pos: 0L, $
         lonlat: 0B, varname: strarr(2), $
         data_range_var0: fltarr(2), data_range_var1: fltarr(2)}

end
