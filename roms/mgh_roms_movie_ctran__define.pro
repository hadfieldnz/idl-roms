;+
; CLASS NAME:
;   Mgh_Roms_Movie_Ctran
;
; PURPOSE:
;   This class generates and displays a an animated sequence of graphs
;   showing a C-transect through a ROMS 2D or 3D output field as a line
;   plot
;
; CATEGORY:
;   Ocean models.
;   Object graphics
;
; CALLING SEQUENCE:
;   mgh_new, 'mgh_roms_movie_ctran', history, variable, slice
;
; POSITIONAL PARAMETERS:
;   history (input, optional?)
;     A synonym for the HISORY_FILE keyword.
;
;   variable (input, optional)
;     The name of a 2-D or 3-D variable in the netCDF file.
;
;   index
;     An integer scalar or array representing slice position(s) in the
;     perpendicular direction. Default is ny/2, where ny is the number
;     of rho points in the perpendicular direction.
;
; KEYWORD PARAMETERS:
;   DATA_RANGE
;     Range of data to be plotted on the y axis. Default is [-1,1]
;
;   DEPTH
;     Set this keyword to a scalar numeric value to specify the depth
;     of a z surface on which data are to be plotted. This keyword
;     should be specified only for variables having a depth coordinate
;     and it cannot be used together with LEVEL.
;
;   DIRECTION
;     Set this equal to 0 for an xi transect and 1 for an eta
;     transect.  Default is 0.
;
;   LEVEL
;     Set this keyword to a scalar integer to specify the s-coordinate
;     level to be plotted.  This keyword should be specified only for
;     variables having a depth coordinate and it cannot be used
;     together with DEPTH.
;
;   SYMBOL
;     Symbols specified via this argument are sized and added to the
;     graph's symbols container.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, Feb 1999:
;     Written as ROMS_TRANSECT_MOVIE.
;   Mark Hadfield, 2000-01:
;     Converted from a procedure (ROMS_TRANSECT_MOVIE) to an object. The front-end
;     procedure (now MGH_ROMS_TRANSECT_MOVIE) can be used to create the object
;   Mark Hadfield, 2000-08:
;     Class name is now MGH_Roms_Transect_Movie and it can be created easily using
;     MGH_NEW. The history file sequence now *must* be in the form of an object.
;   Mark Hadfield, 2001-07:
;     The object is now named Mgh_Roms_Movie_Ctran and uses the C-slice, C-transect
;     functionality built into MGHromsHistory.
;   Mark Hadfield, 2001-10:
;     Updated for IDL 5.5.
;   Mark Hadfield, 2010-07:
;     Minor updates.
;   Mark Hadfield, 2013-05:
;     Modified to accommodate recent changes in the ROMS C-slice code.
;-
function Mgh_Roms_Movie_Ctran::Init, history, variable, index, $
     ALONG_RANGE=along_range, $
     ALL_SLICES=all_slices, $
     AVERAGE_SLICES=average_slices, $
     DATA_RANGE=data_range, DEPTH=depth, $
     DIRECTION=direction, $
     GRAPH_PROPERTIES=graph_properties, $
     HISTORY_FILE=history_file, $
     LEVEL=level, PLOT_COLORS=plot_colors, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     SYMBOL_PROPERTIES=symbol_properties, $
     _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process variable name argument

   if n_elements(variable) eq 0 then variable = 'zeta'

   if size(variable, /TNAME) ne 'STRING' then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', variable

   if n_elements(variable) ne 1 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', variable

   self.variable = variable

   ;; Process history argument.

   case size(history, /TNAME) of
      'STRING': begin
         ohis = obj_new('MGHromsHistory', history)
         self.history_file = ohis
         self.history_destroy = 1
      end
      'OBJREF': begin
         ohis = history
         self.history_file = history
         self.history_destroy = 0
      end
      else: $
         message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', history
   endcase

   ;; Other defaults

   ;; Default data range depends on variable name. This code is replicated
   ;; in several routines. Should it be moved into a single place?

   if n_elements(data_range) eq 0 then begin
      case self.variable of
         'temp'  : data_range=[0,25]
         'salt'  : data_range=[34,36]
         'zeta'  : data_range=[-1,1]
         'u'     : data_range=[-0.3,0.3]
         'v'     : data_range=[-0.3,0.3]
         'ubar'  : data_range=[-0.3,0.3]
         'vbar'  : data_range=[-0.3,0.3]
         'h'     : data_range=[0,4000]
         else    : data_range = [-1,1]
      endcase
   endif

   self.data_range = data_range

   if keyword_set(average_slices) then $
        message, "Sorry I'm not doing AVERAGE_SLICES right now"

   if n_elements(direction) ne 1 then direction = 0

   if n_elements(graph_aspect) eq 0 then graph_aspect = 0.6

   if n_elements(plot_colors) eq 0 then $
        plot_colors = mgh_color(['red','blue','black'])
   numc = n_elements(plot_colors)/3

   ;; Select slice(s) and get grid data for each one

   if n_elements(index) eq 0 then begin
      case direction of
         0: n = ohis->DimInfo('eta_rho', /DIMSIZE)
         1: n = ohis->DimInfo('xi_rho', /DIMSIZE)
      endcase
      index = keyword_set(all_slices) ? lindgen(n) : round(0.5*(n-2))
   endif

   n_slice = n_elements(index)
   grid = ptrarr(n_slice)
   for i=0,n_slice-1 do begin
      grid[i] = ptr_new(ohis->CsliceGrid(ALONG_RANGE=along_range, $
                                         DIRECTION=direction, INDEX=index[i]))
   endfor
   grid0 = *grid[0]

   ;; Establish records to be plotted (if applicable). Establish name of
   ;; time variable--the logic here is a bit shaky, as conventions
   ;; for defining time dimensions and variables in ROMS files are
   ;; not entirely consistent.

   dims = ohis->VarDims(variable)

   has_time = strlen(dims.time) gt 0

   if has_time then begin
      n_time = ohis->DimInfo(dims.time, /DIMSIZE)
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

   ;; Create base graph

   mgh_new, 'MGHgrGraph2D', NAME='ROMS C-transect animation', $
            RESULT=ograph, _STRICT_EXTRA=graph_properties

   ograph->NewFont

   ograph->NewSymbol, 0, /FILL, _STRICT_EXTRA=symbol_properties

   ;; Get minimum and maximum of arc distance

   amin = min(grid0.arc)
   amax = max(grid0.arc)
   for i=1,n_slice-1 do begin
      g = *grid[i]
      amin = amin < min(g.arc)
      amax = amax > max(g.arc)
   endfor

   ;; Draw axes

   ograph->NewAxis, 0, $
        RANGE=[amin,amax], /EXACT, $
        TICKFORMAT='mgh_tf_linear', TICKFRMTDATA={scale:1.E-3, format:'(F10.0)'}
   ograph->NewAxis, 1, $
        RANGE=self.data_range, /EXACT, TITLE=self.variable

   ;; Add atoms to be animated

   otitle = ograph->NewTitle('')

   oplot = objarr(n_slice)
   for i=0,n_slice-1 do begin
      g = *grid[i]
      oplot[i] = ograph->NewAtom('IDLgrPlot', DATAX=g.arc, $
                                 DATAY=fltarr(n_elements(g.arc)), $
                                 COLOR=plot_colors[*,i mod numc], $
                                 SYMBOL=symbol)
   endfor

   ;; Create an animator window to display and manage the movie.

   ok = self->MGH_Datamator::Init(CHANGEABLE=0, GRAPHICS_TREE=ograph, $
                                  _STRICT_EXTRA=_extra)
   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Datamator'

   ;; Create an object array in which to collect command objects

   oframe = objarr(n_slice+1)

   ;; Step through frames, retrieving data & plotting

   for r=0,n_records-1 do begin

      ;; Check to see if the user has selected the "Finish Loading"
      ;; menu item

      if self->Finished() then break

      ;; Label

      if has_time then begin
         t = ohis->VarGet(time_var, OFFSET=records[r], COUNT=[1], /AUTOSCALE)
         oframe[0] = obj_new('MGH_Command', OBJECT=otitle, 'SetProperty', $
                             STRINGS=string(t, FORMAT='(%"%0.3f days")'))
      endif

      ;; Retrieve and plot data

      if keyword_set(average_slices) then begin
        message, "I don't do AVERAGE_SLICES"
      endif else begin
        for i=0,n_slice-1 do begin
          transect = ohis->CtranData(variable, GRID=*grid[i], RECORD=records[r], $
                                      DEPTHS=depth, LEVELS=level)
          oframe[i+1] = obj_new('MGH_Command', OBJECT=oplot[i], 'SetProperty', DATAY=transect)
        endfor
      endelse

      self->AddFrame, oframe

   endfor                       ;; End of loop over records

   self->Finish

   return, 1

end

; Mgh_Roms_Movie_Ctran::Cleanup
;
pro Mgh_Roms_Movie_Ctran::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if self.history_destroy then obj_destroy, self.history_file

   self->MGH_Datamator::Cleanup

end

; Mgh_Roms_Movie_Ctran::GetProperty
;
pro Mgh_Roms_Movie_Ctran::GetProperty, $
     DATA_RANGE=data_range, HISTORY_FILE=history_file, _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   data_range = self.data_range

   history_file = self.history_file

   self->MGH_Datamator::GetProperty, _STRICT_EXTRA=_extra

end

; Mgh_Roms_Movie_Ctran::SetProperty
;
pro Mgh_Roms_Movie_Ctran::SetProperty, $
     DATA_RANGE=data_range, _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

    self->GetProperty, GRAPHICS_TREE=graph

    if n_elements(data_range) gt 0 then $
        if obj_valid(graph) then begin
            self.data_range = data_range
            yaxis = graph->GetAxis(DIRECTION=1)
            yaxis->SetProperty, RANGE=self.data_range
        endif

    self->MGH_Datamator::SetProperty, _STRICT_EXTRA=_extra

end

; Mgh_Roms_Movie_Ctran::About
;
;   Print information about the window and its contents
;
pro Mgh_Roms_Movie_Ctran::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

    self->MGH_Datamator::About, lun

    if obj_valid(self.history_file) then begin
        printf, lun, self, ': the history file sequence is ',self.history_file
        self.history_file->GetProperty, FILE_NAME=file_name
        printf, lun, self, ': the files are:', file_name
    endif

end

; Mgh_Roms_Movie_Ctran::EventMenuBar
;
function Mgh_Roms_Movie_Ctran::EventMenuBar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

    case event.value of

        'TOOLS.SET DATA RANGE': begin
            mgh_new, 'MGH_GUI_SetArray', CAPTION='Data range', CLIENT=self $
                , /FLOATING, GROUP_LEADER=self.base, N_ELEMENTS=2 $
                , PROPERTY_NAME='DATA_RANGE'
            return, 0
        end

        else: return, self->MGH_Datamator::EventMenuBar(event)

    endcase

end

; Mgh_Roms_Movie_Ctran::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Movie_Ctran::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::BuildMenuBar

   ombar = mgh_widget_self(self.menu_bar)

   ombar->NewItem, PARENT='Tools', ['Set Data Range...'], SEPARATOR=[1]

end


pro Mgh_Roms_Movie_Ctran__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, {Mgh_Roms_Movie_Ctran, inherits MGH_Datamator, $
                 history_file: obj_new(), history_destroy: 0B, $
                 variable: '', data_range: fltarr(2)}

end
