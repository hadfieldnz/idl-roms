;+
; CLASS NAME:
;   Mgh_Roms_Plot_Scatter
;
; PURPOSE:
;   This class implements a scatter plot of a pair of vector components
;   from a history OR station file
;
; OBJECT CREATION SEQUENCE
;   mgh_new, 'Mgh_Roms_Plot_Scatter', File, Position
;
; INPUTS:
;   file (input, string or object)
;     The name of a ROMS history file sequence OR a reference to an
;     MGHromsHistory or MGHromsStation object.
;
;   position (input, numeric)
;     For a history file, a [2,n] array specifying the
;     [xi,eta] location of the profiles; for a station file, an integer
;     vector specifying the station number.
;
; KEYWORD PARAMETERS:
;   DEPTH (numeric scalar or vector)
;     Set this keyword to specify the depth at which data are to be
;     extracted. This keyword should be specified only for variables
;     having a depth coordinate and it cannot be used together with
;     LEVEL.
;
;   LEVEL (integer scalar or vector)
;     Set this keyword to specify the s-coordinate level to be
;     plotted.  This keyword should be specified only for variables
;     having a depth coordinate and it cannot be used together with
;     DEPTH.
;
;   VARIABLE (input, string or structure vector)
;     A 2-element vector specifying the names of a pair of u & v
;     type variables. Default is ['u','v'] if available, otherwise
;     ['ubar','vbar'].
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2009-02:
;     Written.
;   Mark Hadfield, 2010-07:
;     Miscellaneous enhancements.
;-
function Mgh_Roms_Plot_Scatter::Init, file, position, $
     DEPTH=depth, LEVEL=level, RECORD_RANGE=record_range, VARIABLE=variable, $
     ELLIPSE_COLORS=ellipse_colors, SCATTER_COLORS=scatter_colors, $
     GRAPH_PROPERTIES=graph_properties, PLOT_PROPERTIES=plot_properties, $
     SHOW_ELLIPSE=show_ellipse, $
     SYMBOL_PROPERTIES=symbol_properties, $
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
         self.file_destroy = 1
      end
      'OBJREF': begin
         self.file = file
         self.file_destroy = 0
      end
      else: message, 'The argument is of the wrong data type'
   endcase
   ofile = self.file

   ;; Process VARIABLE argument

   if n_elements(variable) eq 0 then begin
      has_uv = ofile->HasVar('u') && ofile->HasVar('v')
      variable = has_uv ? ['u','v'] : ['ubar','vbar']
   endif

   self.variable = ptr_new(variable)

   ;; Process appearance keywords

   if n_elements(show_ellipse) eq 0 then show_ellipse = 1B

   if keyword_set(show_ellipse) then begin
      if n_elements(ellipse_colors) eq 0 then $
           ellipse_colors = mgh_color('blue')
      if (size(ellipse_colors, /DIMENSIONS))[0] ne 3 then $
           message, 'ELLIPSE_COLORS must be dimensioned [3,n_colors]'
   endif

   if n_elements(scatter_colors) eq 0 then scatter_colors = mgh_color('red')

   if (size(scatter_colors, /DIMENSIONS))[0] ne 3 then $
        message, 'SCATTER_COLORS must be dimensioned [3,n_colors]'

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

   ;; For each site, extract and save U & V data

   uv = ptrarr(self.n_series)

   for s=0,self.n_series-1 do begin

      case ftype of
         'history': pos = position[*,s < (n_pos-1)]
         'station': pos = position[s < (n_pos-1)]
      endcase

      case 1B of
         n_level gt 0: begin
            data = mgh_roms_series_vector(ofile, pos, RECORD_RANGE=record_range, $
                                          VARIABLE=*self.variable, $
                                          LEVEL=level[s < (n_level-1)])
         end
         n_depth gt 0: begin
            data = mgh_roms_series_vector(ofile, pos, RECORD_RANGE=record_range, $
                                          VARIABLE=*self.variable, $
                                          DEPTH=depth[s < (n_depth-1)])
         end
         else: begin
            data = mgh_roms_series_vector(ofile, pos, RECORD_RANGE=record_range, $
                                          VARIABLE=*self.variable)
         end
      endcase

      ;; Rotate

      data.uv *= exp(complex(0, 1)*data.angle)

      ;; Save uv data

      uv[s] = ptr_new(data.uv)

      mgh_undefine, data

   endfor

   ;; Establish x & y ranges

   xrange = mgh_minmax(real_part(*uv[0]))
   yrange = mgh_minmax(imaginary(*uv[0]))

   for s=1,self.n_series-1 do begin
      xrange[0] = xrange[0] < min(real_part(*uv[s]))
      xrange[1] = xrange[1] > max(real_part(*uv[s]))
      yrange[0] = yrange[0] < min(imaginary(*uv[s]))
      yrange[1] = yrange[1] > max(imaginary(*uv[s]))
   endfor


   ;; Create base graph

   mgh_new, 'MGHgrGraph2D', NAME='ROMS time series', ASPECT=1, $
            _STRICT_EXTRA=graph_properties, RESULT=ograph

   ograph->NewMask

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Draw axes

   ograph->NewAxis, 0, RANGE=xrange, /EXACT, $
        _STRICT_EXTRA=xaxis_properties
   ograph->NewAxis, 1, RANGE=yrange, /EXACT, $
        _STRICT_EXTRA=yaxis_properties

   ograph->NewBackground

   ;; Plot

   for s=0,self.n_series-1 do begin

      ;; Plot points

      c = scatter_colors[*,s mod (n_elements(scatter_colors)/3)]

      ograph->NewSymbol, 0, RESULT=osym, $
           NORM_SIZE=0.005, FILL=0, N_VERT=6, COLOR=c, $
           _STRIcT_EXTRA=symbol_properties

      ograph->NewAtom, 'IDLgrPlot', real_part(*uv[s]), imaginary(*uv[s]), $
           SYMBOL=osym, LINESTYLE=6, COLOR=c, $
           _STRICT_EXTRA=plot_properties

      ;; Optionally plot variance ellipse.

      if keyword_set(show_ellipse) then begin

         c = ellipse_colors[*,s mod (n_elements(ellipse_colors)/3)]

         ell = mgh_variance_ellipse(*uv[s])

         case ftype of
            'history': begin
               fmt = '(%" Position %d %d ellipse parameters: mean %f %f, ' + $
                     'semi-major %f, semi-minor %f, inclination %f deg")'
               message, /INFORM, string(FORMAT=temporary(fmt), pos, $
                                        real_part(ell.mean), imaginary(ell.mean), $
                                        ell.major, ell.minor, ell.angle*!radeg)
            end
            'station': begin
               fmt = '(%" Station %d ellipse parameters: mean %f %f, ' + $
                     'semi-major %f, semi-minor %f, inclination %f deg")'
               message, /INFORM, string(FORMAT=temporary(fmt), pos, $
                                        real_part(ell.mean), imaginary(ell.mean), $
                                        ell.major, ell.minor, ell.angle*!radeg)
            end
         endcase

         cj = complex(0, 1)
         ang = mgh_range([0,2*!pi], N_ELEMENTS=21)
         uve = ell.mean + complex(ell.major*cos(ang), ell.minor*sin(ang))*exp(cj*ell.angle)

         ograph->NewAtom, 'IDLgrPlot', real_part(uve), imaginary(uve), $
              THICK=2, COLOR=c, /USE_ZVALUE, ZVALUE=10*deltaz

      endif

   endfor

   ok = self->MGH_Window::Init(GRAPHICS_TREE=ograph,  _STRICT_EXTRA=extra)

   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   self->Finalize, 'Mgh_Roms_Plot_Scatter'

   return, 1

end

; Mgh_Roms_Plot_Scatter::Cleanup
;
pro Mgh_Roms_Plot_Scatter::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if self.file_destroy then $
        obj_destroy, self.history_file

   self->MGH_Window::Cleanup

end

; Mgh_Roms_Plot_Scatter::GetProperty
;
pro Mgh_Roms_Plot_Scatter::GetProperty, $
     FILE=file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   file = self.file

   self->MGH_Window::GetProperty, _STRICT_EXTRA=extra

end

; Mgh_Roms_Plot_Scatter::SetProperty
;
PRO Mgh_Roms_Plot_Scatter::SetProperty, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->GetProperty, GRAPHICS_TREE=graph

;   if n_elements(data_range) gt 0 then begin
;      self.data_range = data_range
;      if obj_valid(graph) then begin
;         oxax = graph->GetAxis(DIRECTION=0)
;         oxax->SetProperty, RANGE=self.data_range
;         oyax = graph->GetAxis(DIRECTION=1)
;         oyax->SetProperty, RANGE=self.data_range
;      endif
;   endif

   self->MGH_Window::SetProperty, _STRICT_EXTRA=extra

end

; Mgh_Roms_Plot_Scatter::About
;
pro Mgh_Roms_Plot_Scatter::About, lun

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

; Mgh_Roms_Plot_Scatter::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Plot_Scatter::BuildMenuBar

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

; Mgh_Roms_Plot_Scatter::EventMenubar
;
function Mgh_Roms_Plot_Scatter::EventMenubar, event

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

; Mgh_Roms_Plot_Scatter::ExportData
;
pro Mgh_Roms_Plot_Scatter::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::ExportData, values, labels

   self->GetProperty, FILE=file

   labels = [labels, 'File Object']
   values = [values, ptr_new(file)]

end

pro Mgh_Roms_Plot_Scatter__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {Mgh_Roms_Plot_Scatter, inherits MGH_Window, $
         file: obj_new(), file_destroy: 0B, n_series: 0L, $
         variable: ptr_new()}

end
