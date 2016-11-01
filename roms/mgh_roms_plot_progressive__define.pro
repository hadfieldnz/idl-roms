;+
; CLASS NAME:
;   Mgh_Roms_Plot_Progressive
;
; PURPOSE:
;   This class implements a progressive-vector plot from a history OR
;   station file
;
; OBJECT CREATION SEQUENCE
;   mgh_new, 'Mgh_Roms_Plot_Progressive', File, Var, Position
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
;     Set this keyword to specify the
;     depth at which data are to be extracted. This
;     keyword should be specified only for variables having a depth
;     coordinate and it cannot be used together with LEVEL.
;
;   LEVEL (integer scalar or vector)
;     Set this keyword to specify the
;     s-coordinate level to be plotted.  This keyword should be
;     specified only for variables having a depth coordinate and it
;     cannot be used together with DEPTH.
;
;   VARIABLE (input, string or structure vector)
;     A 2-element vector specifying the names of a pair of u & v
;     type variables. Default is ['u','v'] if available, otherwise
;     ['ubar','vbar'].
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2008-06:
;     Written.
;   Mark Hadfield, 2009-01:
;     Cleaned up documentation.
;   Mark Hadfield, 2009-02:
;     - Added support for fractional positions for history files.
;     - Now gets velocity data from the file using MGH_ROMS_SERIES_VECTOR.
;   Mark Hadfield, 2012-09:
;     - Removed FILE_DESTROY functionality.
;-
function mgh_roms_plot_progressive::Init, file, position, $
     DEPTH=depth, LEVEL=level, SIGMA=sigma, $
     ASPECT=aspect, PLOT_COLORS=plot_colors, $
     RECORD_RANGE=record_range, VARIABLE=variable, $
     GRAPH_PROPERTIES=graph_properties, $
     PLOT_PROPERTIES=plot_properties, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
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

   if n_elements(variable) eq 0 then begin
      has_uv = ofile->HasVar('u') && ofile->HasVar('v')
      variable = has_uv ? ['u','v'] : ['ubar','vbar']
   endif

   self.variable = ptr_new(variable)

   ;; Process plot colours

   if n_elements(plot_colors) eq 0 then $
        plot_colors = mgh_color('red')

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
   n_sigma = n_elements(sigma)

   self.n_series =  n_pos > n_depth > n_level > n_sigma

   ;; Extract U & V data and calculate pseudo-trajectories (in kilometres)

   xprog = ptrarr(self.n_series)
   yprog = ptrarr(self.n_series)

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
         n_sigma gt 0: begin
            data = mgh_roms_series_vector(ofile, pos, RECORD_RANGE=record_range, $
                                          VARIABLE=*self.variable, $
                                          SIGMA=sigma[s < (n_sigma-1)])
         end
         else: begin
            data = mgh_roms_series_vector(ofile, pos, RECORD_RANGE=record_range, $
                                          VARIABLE=*self.variable)
         end
      endcase

      ;; Rotate

      data.uv *= exp(complex(0, 1)*data.angle)

      ;; Calculate pseudo-trajectory

      n_time = n_elements(data.time)
      dt = mgh_diff(mgh_stagger(data.time, DELTA=1))

      u = real_part(data.uv)
      v = imaginary(data.uv)

      xx = fltarr(n_time+1)
      yy = fltarr(n_time+1)
      for i=0,n_time-1 do begin
         xx[i+1] = xx[i] + 1.E-3*dt[i]*u[i]*(24*3600)
         yy[i+1] = yy[i] + 1.E-3*dt[i]*v[i]*(24*3600)
      endfor

      xprog[s] = ptr_new(xx, /NO_COPY)
      yprog[s] = ptr_new(yy, /NO_COPY)

   endfor

   ;; Set up default axis ranges so that the data aspect ratio
   ;; matches the graph aspect ratio...

   ;; ...establish ranges of x & y data.

   xrange = mgh_minmax(*xprog[0])
   yrange = mgh_minmax(*yprog[0])

   for s=1,self.n_series-1 do begin
      xrange[0] = xrange[0] < min(*xprog[s])
      xrange[1] = xrange[1] > max(*xprog[s])
      yrange[0] = yrange[0] < min(*yprog[s])
      yrange[1] = yrange[1] > max(*yprog[s])
   endfor

   ;; ...tweak x & y ranges

   if xrange[0] eq xrange[1] then xrange = [-1,1]
   if yrange[0] eq yrange[1] then yrange = [-1,1]

   xrange += [-1,1]*0.05*(xrange[1]-xrange[0])
   yrange += [-1,1]*0.05*(yrange[1]-yrange[0])

   ;; ...adjust x & y ranges

   if n_elements(aspect) eq 0 then aspect = 1

   dx = xrange[1]-xrange[0]
   dy = yrange[1]-yrange[0]

   daspect = dy/dx

   if daspect gt aspect then begin
      xrange += [-1,1]*0.5*(dy/aspect-dx)
   endif else begin
      yrange += [-1,1]*0.5*(dx*aspect-dy)
   endelse

   ;; Create base graph

   mgh_new, 'MGHgrGraph2D', NAME='ROMS progressive vectors', ASPECT=aspect, $
            _STRICT_EXTRA=graph_properties, RESULT=ograph

   ograph->NewMask

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Draw axes.

   ograph->NewAxis, 0, RANGE=xrange, /EXACT, TITLE='X (km)', $
        _STRICT_EXTRA=xaxis_properties
   ograph->NewAxis, 1, RANGE=yrange, /EXACT, TITLE='Y (km)', $
        _STRICT_EXTRA=yaxis_properties

   ;; Add a background to act as a selection target

   ograph->NewBackground

   ;; Draw plots

   for s=0,self.n_series-1 do begin

      ograph->NewAtom, 'IDLgrPlot', *xprog[s], *yprog[s], $
           COLOR=plot_colors[*,s mod n_colors], _STRICT_EXTRA=plot_properties

   endfor

   ptr_free, xprog, yprog

   ;; Open window and finalise

   ok = self->MGH_Window::Init(GRAPHICS_TREE=ograph,  _STRICT_EXTRA=extra)

   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   self->Finalize, 'Mgh_Roms_Plot_Progressive'

   return, 1

end

; Mgh_Roms_Plot_Scatter::Cleanup
;
pro Mgh_Roms_Plot_Scatter::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.variable

   self->MGH_Window::Cleanup

end

; Mgh_Roms_Plot_Progressive::GetProperty
;
PRO Mgh_Roms_Plot_Progressive::GetProperty, $
     FILE=file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   file = self.file

   self->MGH_Window::GetProperty, _STRICT_EXTRA=extra

END

; Mgh_Roms_Plot_Progressive::SetProperty
;
PRO Mgh_Roms_Plot_Progressive::SetProperty, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::SetProperty, _STRICT_EXTRA=extra

end

; Mgh_Roms_Plot_Progressive::About
;
pro Mgh_Roms_Plot_Progressive::About, lun

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

; Mgh_Roms_Plot_Progressive::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Plot_Progressive::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::BuildMenuBar

end

; Mgh_Roms_Plot_Progressive::EventMenubar
;
function Mgh_Roms_Plot_Progressive::EventMenubar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      else: return, self->MGH_Window::EventMenubar(event)

   endcase

end

; Mgh_Roms_Plot_Progressive::ExportData
;
pro Mgh_Roms_Plot_Progressive::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::ExportData, values, labels

   self->GetProperty, FILE=file

   labels = [labels, 'File Object']
   values = [values, ptr_new(file)]

end

pro Mgh_Roms_Plot_Progressive__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {Mgh_Roms_Plot_Progressive, inherits MGH_Window, $
         file: obj_new(), n_series: 0L, $
         variable: ptr_new()}

end
