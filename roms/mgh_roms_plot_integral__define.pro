;+
; CLASS NAME:
;   Mgh_Roms_Plot_Integral
;
; PURPOSE:
;   This class implements a time-series plot of a volume integral
;   from a history file
;
; OBJECT CREATION SEQUENCE
;   mgh_new, 'Mgh_Roms_Plot_integral', file, VARIABLE=variable
;
; INPUTS:
;   history
;     The name of a ROMS history file OR a reference to an
;     MGHromsHistory object.
;
;   var
;     Variable name
;
; KEYWORD PARAMETERS:
;   DATA_RANGE
;     If supplied, this keyword specifies the data (y-axis) range.
;
;   VARIABLE
;     A variable descriptor, must be a string or a structure that can
;     be interpreted by the MGHromsHistory or object's methods.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2011-08:
;     Written.
;-
function Mgh_Roms_Plot_Integral::Init, file, $
     DATA_RANGE=data_range, $
     GRAPH_PROPERTIES=graph_properties, $
     PLOT_PROPERTIES=plot_properties, $
     RECORD_RANGE=record_range, VARIaBLE=variable, $
     XAXIS_PROPERTIES=xaxis_properties, YAXIS_PROPERTIES=yaxis_properties, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process file argument

   ;; Process history argument.

   case size(history, /TYPE) of
      7: begin  ;;; String
         ohis = obj_new('MGHromsHistory', history)
         self.history_file = ohis
         self.history_destroy = 1B
      end
      11: begin   ;;; Object reference
         ohis = history
         self.history_file = history
         self.history_destroy = 0B
      end
      0: $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'history'
      else: $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'history'
   endcase

   ;; Process VARIABLE argument

   if n_elements(variable) eq 0 then variable = 'zeta'

   self.variable = ptr_new(variable)

   ;; Default data range depends on variable name

   mgh_roms_resolve_data, *self.variable, DATA_RANGE=data_range

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
         n_pos = n_elements(position)/2
      end
      'station': begin
         if n_elements(position) eq 0 then position = 0
         n_pos = n_elements(position)
      end
   endcase

   ;; Calculate number of series

   n_level = n_elements(level)
   n_depth = n_elements(depth)

   self.n_series =  n_pos > n_level > n_depth

   ;; Create base graph

   mgh_new, 'MGHgrGraph2D', NAME='ROMS time series', ASPECT=0.4, $
            _STRICT_EXTRA=graph_properties, RESULT=ograph

   ograph->NewMask

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

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
                                          LEVEL=level[s < (n_level-1)])
         end
         n_depth gt 0: begin
            data = mgh_roms_series_scalar(ofile, pos, VARIaBLE=*self.variable, $
                                          DEPTH=depth[s < (n_depth-1)])
         end
         else: begin
            data = mgh_roms_series_scalar(ofile, pos, VARIaBLE=*self.variable)
         end
      endcase

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

   self->Finalize, 'Mgh_Roms_Plot_integral'

   return, 1

end

; Mgh_Roms_Plot_integral::Cleanup
;
pro Mgh_Roms_Plot_integral::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.variable

   self->MGH_Window::Cleanup

end

; Mgh_Roms_Plot_integral::GetProperty
;
PRO Mgh_Roms_Plot_integral::GetProperty, $
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

; Mgh_Roms_Plot_integral::SetProperty
;
PRO Mgh_Roms_Plot_integral::SetProperty, $
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

; Mgh_Roms_Plot_integral::About
;
pro Mgh_Roms_Plot_integral::About, lun

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

; Mgh_Roms_Plot_integral::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Plot_integral::BuildMenuBar

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

; Mgh_Roms_Plot_integral::EventMenubar
;
function Mgh_Roms_Plot_integral::EventMenubar, event

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

; Mgh_Roms_Plot_integral::ExportData
;
pro Mgh_Roms_Plot_integral::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::ExportData, values, labels

   self->GetProperty, ANIMATION=animation, FILE=file, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[0]->GetProperty, KEYWORDS=keywords

   labels = [labels, 'File Object', 'Profile Depth', 'Profile Data']
   values = [values, ptr_new(file), ptr_new(keywords.datay), $
             ptr_new(keywords.datax)]

end

pro Mgh_Roms_Plot_integral__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {Mgh_Roms_Plot_integral, inherits MGH_Window, $
         history: obj_new(), variable: ptr_new(), $
         data_range: fltarr(2)}

end
