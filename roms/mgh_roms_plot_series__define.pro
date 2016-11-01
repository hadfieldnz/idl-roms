;+
; CLASS NAME:
;   Mgh_Roms_Plot_Series
;
; PURPOSE:
;   This class implements a time-series plot from a history OR station
;   file
;
; OBJECT CREATION SEQUENCE
;   mgh_new, 'Mgh_Roms_Plot_Series', File, Var, Position
;
; INPUTS:
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
;   DATA_MULTIPLIER (input, numeric scalar)
;     Number by which data values are multiplied before they are loaded into
;     the density surface. Default depends on the variable and is calculated by
;     MGH_ROMS_RESOLVE_DATA
;
;   DATA_RANGE
;     If supplied, this keyword specifies the data (y-axis) range. Default depends
;     on the variable and is calculated by MGH_ROMS_RESOLVE_DATA
;
;   DEPTH
;     Set this keyword to a numeric scalar or vector to specify the
;     depth(s) at which data are to be extracted. This
;     keyword should be specified only for variables having a depth
;     coordinate and it cannot be used together with LEVEL.
;
;   LAYER
;     Set this keyword to an integer scalar or vector to specify the
;     bed layer(s) at which data are to be plotted. This keyword should
;     be specified only for variables having a bed-layer dimension. It is
;     currently supported only for history files.
;
;   LEVEL
;     Set this keyword to an integer scalar or vector to specify the
;     s-coordinate level(s) to be plotted.  This keyword should be
;     specified only for variables having a depth coordinate and it
;     cannot be used together with DEPTH.
;
;   VARIABLE
;     A variable descriptor, must be a string or a structure that can
;     be interpreted by the MGHromsHistory or MGHromsStation object's
;     methods.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 1999-02:
;     Written.
;   Mark Hadfield, 2001-07:
;     Extensive modifications. Updated for IDL 5.5. Now can plot
;     series from history OR station files. Supports multiple series.
;   Mark Hadfield, 2009-02:
;     Most of thw code to extract the data has now been moved into
;     MGH_ROMS_SERIES_SCALAR.
;   Mark Hadfield, 2011-08:
;     No longer destroys any mghromsfile or mghromsstation objects that
;     are created.
;   Mark Hadfield, 2013-06:
;     Added DATA_MULTIPLIER keyword.
;-
function Mgh_Roms_Plot_Series::Init, file, position, $
     DATA_MULTIPLIER=data_multiplier, $
     DATA_RANGE=data_range, $
     DEPTH=depth, LAYER=layer, LEVEL=level, SIGMA=sigma, $
     NEAREST=nearest, $
     GRAPH_PROPERTIES=graph_properties, PLOT_COLORS=plot_colors, $
     PLOT_PROPERTIES=plot_properties, $
     RECORD_RANGE=record_range, TIME_RANGE=time_range, $
     VARIaBLE=variable, $
     XAXIS_PROPERTIES=xaxis_properties, YAXIS_PROPERTIES=yaxis_properties, $
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
      end
      'OBJREF': begin
         self.file = file
      end
      else: message, 'The argument is of the wrong data type'
   endcase
   ofile = self.file

   ;; Process VARIABLE argument

   if n_elements(variable) eq 0 then variable = 'zeta'

   self.variable = ptr_new(variable)

   ;; Other defaults

   mgh_roms_resolve_data, self.variable, $
        DATA_RANGE=data_range, DATA_MULTIPLIER=data_multiplier

   self.data_range = data_range

   ;; Process plot colours

   if n_elements(plot_colors) eq 0 then plot_colors = mgh_color('red')

   if (size(plot_colors, /DIMENSIONS))[0] ne 3 then $
        message, 'PLOT_COLORS must be dimensioned [3,n_colors]'

   n_colors = n_elements(plot_colors)/3

   ;; Process position parameter

   case 1B of
      obj_isa(ofile, 'MGHromsHistory'): ftype = 'history'
      obj_isa(ofile, 'MGHromsStation'): ftype = 'station'
      else: message, 'Unknown file type'
   endcase

   case ftype of
      'history': begin
         if n_elements(position) eq 0 then begin
            dim_rho = ofile->DimRho()
            position = 0.5*(dim_rho[0:1]-1)
         endif
         if ~ isa(position, /NUMBER) then $
            message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'position'
         n_pos = n_elements(position)/2
      end
      'station': begin
         if n_elements(position) eq 0 then position = 0
         if ~ isa(position, /NUMBER) then $
            message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'position'
         n_pos = n_elements(position)
      end
   endcase

   ;; Calculate number of series

   n_depth = n_elements(depth)
   n_layer = n_elements(layer)
   n_level = n_elements(level)
   n_sigma = n_elements(sigma)

   self.n_series =  n_pos > n_depth > n_layer > n_level > n_sigma

   ;; Create base graph

   mgh_new, 'MGHgrGraph2D', NAME='ROMS time series', ASPECT=0.4, $
            _STRICT_EXTRA=graph_properties, RESULT=ograph

   ograph->NewMask

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ograph->NewTitle

   ;; Draw vertical axis. (The time axis will be drawn after data for
   ;; the first series has been extracted, so we know the range.)

   ograph->NewAxis, 1, RANGE=data_range, /EXACT, $
        _STRICT_EXTRA=yaxis_properties

   ograph->NewBackground

   ;; Work through the time series, extracting & plotting

   for s=0,self.n_series-1 do begin

      case ftype of
         'history': pos = position[*,s < (n_pos-1)]
         'station': pos = position[s < (n_pos-1)]
      endcase

      case 1B of
         n_level gt 0: begin
            data = mgh_roms_series_scalar(ofile, pos, VARIaBLE=*self.variable, $
                                          RECORD_RANGE=record_range, TIME_RANGE=time_range, $
                                          LEVEL=level[s < (n_level-1)], NEAREST=nearest)
         end
         n_depth gt 0: begin
            data = mgh_roms_series_scalar(ofile, pos, VARIaBLE=*self.variable, $
                                          RECORD_RANGE=record_range, TIME_RANGE=time_range, $
                                          DEPTH=depth[s < (n_depth-1)], NEAREST=nearest)
         end
         n_sigma gt 0: begin
            data = mgh_roms_series_scalar(ofile, pos, VARIaBLE=*self.variable, $
                                          RECORD_RANGE=record_range, TIME_RANGE=time_range, $
                                          SIGMA=sigma[s < (n_sigma-1)], NEAREST=nearest)
         end
         n_layer gt 0: begin
            data = mgh_roms_series_scalar(ofile, pos, VARIaBLE=*self.variable, $
                                          RECORD_RANGE=record_range, TIME_RANGE=time_range, $
                                          LAYER=layer[s < (n_layer-1)], NEAREST=nearest)
         end
         else: begin
            data = mgh_roms_series_scalar(ofile, pos, VARIaBLE=*self.variable, $
                                          RECORD_RANGE=record_range, TIME_RANGE=time_range, $
                                          NEAREST=nearest)
         end
      endcase

      data.value *= data_multiplier

      ;; Draw time axis

      if s eq 0 then $
           ograph->NewAxis, 0, RANGE=mgh_minmax(data.time), /EXACT, TITLE='Time', $
                _STRICT_EXTRA=xaxis_properties


      ;; Plot

      ograph->NewAtom, 'IDLgrPlot', data.time, data.value, $
           COLOR=plot_colors[*,s mod n_colors], _STRICT_EXTRA=plot_properties

      mgh_undefine, data

   endfor

   ok = self->MGH_Window::Init(GRAPHICS_TREE=ograph,  _STRICT_EXTRA=extra)

   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   self->Finalize, 'Mgh_Roms_Plot_Series'

   return, 1

end

; Mgh_Roms_Plot_Series::Cleanup
;
pro Mgh_Roms_Plot_Series::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.variable

   self->MGH_Window::Cleanup

end

; Mgh_Roms_Plot_Series::GetProperty
;
PRO Mgh_Roms_Plot_Series::GetProperty, $
     DATA_RANGE=data_range, FILE=file, VARIABLE=variable, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   data_range = self.data_range

   file = self.file

   variable = *self.variable

   self->MGH_Window::GetProperty, _STRICT_EXTRA=extra

END

; Mgh_Roms_Plot_Series::SetProperty
;
PRO Mgh_Roms_Plot_Series::SetProperty, $
     DATA_RANGE=data_range, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->GetProperty, GRAPHICS_TREE=graph

   if n_elements(data_range) gt 0 then begin
      self.data_range = data_range
      if obj_valid(graph) then begin
         oyaxis = graph->GetAxis(DIRECTION=1)
         oyaxis->SetProperty, RANGE=self.data_range
      endif
   endif

   self->MGH_Window::SetProperty, _STRICT_EXTRA=extra

end

; Mgh_Roms_Plot_Series::About
;
pro Mgh_Roms_Plot_Series::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::About, lun

   if obj_valid(self.file) then begin
      printf, lun, self, ': the history/station file object is ', self.file
      self.file->GetProperty, FILE_NAME=file_name
      printf, lun, self, ': the files are:', file_name
   endif

end

; Mgh_Roms_Plot_Series::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Plot_Series::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::BuildMenuBar

   ombar = mgh_widget_self(self.menu_bar)

   if obj_valid(ombar) then begin
      ombar->NewItem, PARENT='Tools', SEPARATOR=[1,1], $
        ['Set Data Range...']
   endif

end

; Mgh_Roms_Plot_Series::EventMenubar
;
function Mgh_Roms_Plot_Series::EventMenubar, event

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

      else: return, self->MGH_Window::EventMenubar(event)

   endcase

end

; Mgh_Roms_Plot_Series::ExportData
;
pro Mgh_Roms_Plot_Series::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::ExportData, values, labels

end

pro Mgh_Roms_Plot_Series__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {Mgh_Roms_Plot_Series, inherits MGH_Window, $
         file: obj_new(), variable: ptr_new(), n_series: 0L, $
         data_range: fltarr(2)}

end
