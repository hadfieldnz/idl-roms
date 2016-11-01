;+
; CLASS NAME:
;   MGH_ROMS_Movie_Tport
;
; PURPOSE:
;   This class generates and displays a an animated sequence of graphs
;   showing volume transport vs distance around a box enclosing a ROMS
;   simulation.
;
; CALLING SEQUENCE:
;   mgh_new, 'mgh_roms_Movie_Tport', history
;
; POSITIONAL PARAMETERS:
;   history
;     A reference to a ROMS history sequence object or a string arrary
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
; KEYWORD PARAMETERS:
;   DATA_RANGE
;     Data range for the transport plot. Default is [-30,30]
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2002-10:
;     Written.
;-
function MGH_ROMS_Movie_Tport::Init, history, $
     ANIMATION_PROPERTIES=animation_properties, $
     DATA_RANGE=data_range, $
     ETA_RANGE=eta_range, $
     GRAPH_PROPERTIES=graph_properties, $
     PALETTE=palette, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     RECORDS=records, $
     XI_RANGE=xi_range, $
     USE_ZETA=use_zeta, $
     VAR_UBAR=var_ubar, VAR_VBAR=var_vbar, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process history argument.

   case size(history, /TNAME) of
      'STRING': begin
         ohis = obj_new('MGHromsHistory', history)
         self.history_file = ohis
      end
      'OBJREF': begin
         ohis = history
         self.history_file = history
      end
      else: message, 'The argument is of the wrong data type'
   endcase

   ;; Set other keyword defaults

   if n_elements(data_range) eq 0 then data_range = [-100,100]

   ;; Check variable names

   if n_elements(var_ubar) eq 0 then var_ubar = 'ubar'
   if n_elements(var_vbar) eq 0 then var_vbar = 'vbar'

   ;; Establish records to be plotted (if applicable)

   dim_u = ohis->VarDims(var_ubar)

   has_time = strlen(dim_u.time) gt 0

   if has_time then begin
      n_time = ohis->DimInfo(dim_u.time, /DIMSIZE)
      mgh_resolve_indices, n_time, record_range, record_stride, records
      n_records = n_elements(records)
      time_var = ohis->TimeVarName(dim_u.time)
      if isa(time_var, /NULL) then message, 'Time variable not found'
   endif else begin
      n_records = 1
      ;; Check other record keywords for validity here?
      if record_average gt 1 then $
           message, 'Cannot average over records for time-independent variable'
   endelse

   ;; Get transport data for RECORD 0 so we can set up the graph

   tport = ohis->GetTransportBox(RECORD=0, VAR_UBAR=var_ubar, VAR_VBAR=var_vbar)

   land = where(~ finite(tport.depth), n_land)
   if n_land gt 0 then tport.depth[land] = 0
   mgh_undefine, land, n_land

   ;; Create base graph

   ograph = obj_new('MGHgrGraph2D', ASPECT=0.7, $
                    NAME='ROMS transport '+mgh_get_property(ohis, /NAME), $
                    _STRICT_EXTRA=graph_properties)

   ograph->GetProperty, FONTSIZE=fontsize, PLOT_RECT=prect

   ograph->NewFont, SIZE=0.9*fontsize

   del = 0.025

   ;; Show boundary number & depth

   rect = ograph->Rect([0,0.8+del,1,0.2-2*del])

   ograph->NewAxis, 0, RANGE=mgh_minmax(tport.distance*1.E-3), /EXACT, $
        NOTEXT=1, RECT=rect, RESULT=oxaxis
   ograph->NewAxis, 1, RANGE=[0.5,4.5], /EXACT, $
        TICKVALUES=[1,2,3,4], MINOR=0, RECT=rect, $
        RESULT=oyaxis

   ograph->NewAtom, 'IDLgrPlot', COLOR=mgh_color('red'), $
        tport.distance*1.E-3, tport.boundary+1, $
        XAXIS=oxaxis, YAXIS=oyaxis

   rect = ograph->Rect([0,0.6+del,1,0.2-del])

   ograph->NewAxis, 0, RANGE=mgh_minmax(tport.distance*1.E-3), /EXACT, $
        NOTEXT=1, RECT=rect, RESULT=oxaxis
   ograph->NewAxis, 1, RANGE=[-6000,0], /EXACT, /EXTEND, $
        TICKVALUES=[-6000,-4000,-2000,0], MINOR=1, RECT=rect, $
        RESULT=oyaxis

   ograph->NewAtom, 'IDLgrPlot', COLOR=mgh_color('blue'), $
        tport.distance*1.E-3, -tport.depth, $
        XAXIS=oxaxis, YAXIS=oyaxis

   ;; Set up transport axes and create a plot object to be
   ;; animated.

   rect = [prect[0],prect[1]+del*prect[3], $
           prect[2],(0.6-del)*prect[3]]

   ograph->NewAxis, 0, RANGE=mgh_minmax(tport.distance*1.E-3), /EXACT, $
        RECT=rect, RESULT=oxaxis
   ograph->NewAxis, RESULT=oyaxis, 1, RANGE=data_range, /EXACT, RECT=rect, $
        _STRICT_EXTRA=yaxis_properties
   self.yaxis = oyaxis[0]

   ograph->NewAtom, 'IDLgrPlot', THICK=2, $
        tport.distance*1.E-3, replicate(0,n_elements(tport.distance)), $
        XAXIS=oxaxis, YAXIS=oyaxis, RESULT=oplot
   self.plot = oplot

   ;; Add a title object to be animated

   otitle = ograph->NewTitle('')

   ;; Create an animator window to display and manage the movie.

   ok = self->MGH_Datamator::Init(CHANGEABLE=0, GRAPHICS_TREE=ograph, $
                                  _STRICT_EXTRA=extra)
   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Datamator'

   ;; Step through times in the history file, generating new frames &
   ;; plotting data

   oframe = objarr(2)

   for r=0,n_records-1 do begin

      ;; Check to see if the user has selected the "Finish Loading"
      ;; menu item

      if self->Finished() then break

      ;; Retrieve & plot slice.

      if has_time then rec = records[r]

      tport = ohis->GetTransportBox(RECORD=rec, VAR_UBAR=var_ubar, VAR_VBAR=var_vbar)

      oframe[0] = obj_new('MGH_Command', OBJECT=oplot, 'SetProperty', $
                          DATAY=tport.transport*1.E-6)
      ;; Title

      if has_time then begin
         t = ohis->VarGet(time_var, OFFSET=rec, COUNT=[1], /AUTOSCALE)
         oframe[1] = obj_new('MGH_Command', OBJECT=otitle, 'SetProperty', $
                             STRINGS=string(t, FORMAT='(F0.3)')+' days')
      endif

      ;; Add frame to animator & display

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; MGH_ROMS_Movie_Tport::Cleanup
;
pro MGH_ROMS_Movie_Tport::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR

   self->MGH_Datamator::Cleanup

end

; MGH_ROMS_Movie_Tport::GetProperty
;
pro MGH_ROMS_Movie_Tport::GetProperty, $
     ALL=all, DATA_RANGE=data_range, HISTORY_FILE=history_file, PLOT=plot, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR

   self->MGH_Datamator::GetProperty, ALL=all, _STRICT_EXTRA=extra

   if arg_present(all) or arg_present(data_range) then $
        self.yaxis->GetProperty, RANGE=data_range

   history_file = self.history_file

   plot = self.plot

   if arg_present(all) then $
        all = create_struct(all, 'data_range', data_range, $
                            'history_file', history_file, $
                            'plot', plot)

end

; MGH_ROMS_Movie_Tport::SetProperty
;
pro MGH_ROMS_Movie_Tport::SetProperty, $
     DATA_RANGE=data_range, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR

   if n_elements(data_range) gt 0 then $
        self.yaxis->SetProperty, RANGE=data_range

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=extra

end

; MGH_ROMS_Movie_Tport::About
;
;   Print information about the window and its contents
;
pro MGH_ROMS_Movie_Tport::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR

   self->MGH_Datamator::About, lun

   self->GetProperty, HISTORY_FILE=history_file

   if obj_valid(history_file) then begin
      printf, lun, self, ': the history file sequence is ', history_file
      history_file->GetProperty, FILE_NAME=file_name
      printf, lun, self, ': the files are:', file_name
   endif

end

; MGH_ROMS_Movie_Tport::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro MGH_ROMS_Movie_Tport::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR

    self->MGH_Datamator::BuildMenuBar

    obar = mgh_widget_self(self.menu_bar)

    obar->NewItem, PARENT='Tools', $
      ['Set Data Range...','View Data Values...']

end


; MGH_ROMS_Movie_Tport::EventMenuBar
;
function MGH_ROMS_Movie_Tport::EventMenuBar, event

   compile_opt DEFINT32
   compile_opt STRICTARR

   case event.value of

      'TOOLS.SET DATA RANGE': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Range', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, IMMEDIATE=0, $
                  N_ELEMENTS=2, PROPERTY_NAME='DATA_RANGE'
         return, 0
      end

      'TOOLS.VIEW DATA VALUES': begin
         self->GetProperty, POSITION=position
         oframe = self.animation->GetFrame(POSITION=position)
         oframe[0]->GetProperty, KEYWORDS=keywords
         n_data = size(keywords.datay, /N_ELEMENTS)
         ;; Call REFORM so that XVAREDIT cannot modify values
         xvaredit, reform(keywords.datay, 1, n_data), GROUP=self.base, $
                   Y_SCROLL_SIZE=(n_data < 30)
         return, 0
      end

      else: return, self->MGH_Datamator::EventMenuBar(event)

   endcase

end

; MGH_ROMS_Movie_Tport::ExportData
;
pro MGH_ROMS_Movie_Tport::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR

   self->MGH_Player::ExportData, values, labels

   self->GetProperty, ANIMATION=animation, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[0]->GetProperty, KEYWORDS=keywords

   labels = [labels, 'Transport Data']
   values = [values, ptr_new(keywords.datay)]

end

pro MGH_ROMS_Movie_Tport__Define

   compile_opt DEFINT32
   compile_opt STRICTARR

   struct_hide, {MGH_ROMS_Movie_Tport, inherits MGH_Datamator, $
                 history_file: obj_new(), $
                 plot: obj_new(), yaxis: obj_new()}

end
