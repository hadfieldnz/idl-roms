; svn $Id$
;+
; CLASS NAME:
;   MGH_ROMS_Ztraj_Plot
;
; PURPOSE:
;   This procedure displays ROMS vertital float trajectories in (z,t) space
;
; CALLING SEQUENCE:
;   mgh_new, 'MGH_ROMS_Ztraj_Plot', ffile, hfile
;
; SUPERCLASS:
;   MGH_Window
;
; INIT ARGUMENTS:
;  ffile
;    Synonym for FLOAT_FILE property.
;
;  hfile
;    Synonym for HISTORY_FILE property.
;
; PROPERTIES:
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
;  SYMBOL_PROPERTIES (Init)
;    A structure containing keywords to be passed to the plotting
;    symbol.
;
;  ZAXIS_PROPERTIES (Init)
;  TAXIS_PROPERTIES (Init)
;    Structures containing keywords to be passed to the z (vertical) &
;    t (horizontal) axes.
;
;  Z_RANGE (Init)
;    A 2-element integer vector specifying the range of the vertical
;    (z) axis. The default is max(h)*[-1,0] and the axis has the
;    EXTEND property set by default.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2002-02:
;     Written.
;   Mark Hadfield, 2009-10:
;     Removed calls to widget_event(/NOWAIT).
;-

function MGH_ROMS_Ztraj_Plot::Init, ffile, hfile, $
     FLOAT_FILE=float_file, $
     FLOAT_RANGE=float_range, $
     FLOAT_STRIDE=float_stride, $
     FLOATS=floats, $
     GRAPH_PROPERTIES=graph_properties, $
     HISTORY_FILE=history_file, $
     PLOT_PROPERTIES=plot_properties, $
     SHOW_ORIGIN=show_origin, $
     TAXIS_PROPERTIES=taxis_properties, $
     TIME_RANGE=time_range, $
     Z_RANGE=xi_range, $
     ZAXIS_PROPERTIES=zaxis_properties, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process float-file argument

   if n_elements(float_file) eq 0 && n_elements(ffile) gt 0 then $
        float_file = ffile

   if n_elements(float_file) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'float_file'

   case size(float_file, /TNAME) of
      'STRING': begin
         self.float_file = obj_new('MGHromsFloat', float_file)
         self.float_destroy = 1
      end
      'OBJREF': begin
         self.float_file = float_file
         self.float_destroy = 0
      end
   endcase
   oflt = self.float_file

   ;; Process history-file argument

   if n_elements(history_file) eq 0 && n_elements(hfile) gt 0 then $
        history_file = hfile

   if n_elements(history_file) eq 0 then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_UNDEFVAR', 'history_file'

   case size(history_file, /TNAME) of
      'STRING': begin
         self.history_file = obj_new('MGHromsHistory', history_file)
         self.history_destroy = 1B
      end
      'OBJREF': begin
         self.history_file = history_file
         self.history_destroy = 0B
      end
   endcase
   ohis = self.history_file

   ;; Set defaults

   if n_elements(z_range) eq 0 then $
        z_range = [-1.02,0.02]*max(ohis->VarGet('h'))

   ;; Establish variable and dimension name associated with float data

   case 1B of
      oflt->HasVar('depth') && oflt->HasAtt('depth','time'): $
           time_var = oflt->AttGet('depth','time')
      else: $
           time_var = 'ocean_time'
   endcase

   time_dim = (oflt->VarDimNames(time_var))[0]

   ;; Get time and set default time range. TO DO: scale result
   ;; according to "units" property

   time = oflt->VarGet(time_var)/(24.D0*3600)

   if n_elements(time_range) eq 0 then time_range = mgh_minmax(time)

   ;; Get number of floats and resolve float-selection arguments

   n_float = oflt->DimInfo('drifter', /DIMSIZE)

   mgh_resolve_indices, n_float, float_range, float_stride, floats

   ;; Get grid dimensions

   dim_rho = ohis->DimRho()

   ;; Create base graph

   ograph = obj_new('MGHgrGraph2D', ASPECT=0.6, NAME='ROMS float trajectory', $
                    _STRICT_EXTRA=graph_properties)

   ograph->GetProperty, DELTAZ=deltaz

   ograph->NewFont, SIZE=10

   ograph->NewMask

   ;; Draw axes

   ograph->NewAxis, DIRECTION=0, RANGE=time_range, /EXACT, $
        TITLE='Time (days)', _STRICT_EXTRA=taxis_properties

   ograph->NewAxis, DIRECTION=1, RANGE=z_range, /EXACT, /EXTEND, $
        TITLE='Depth (m)', TICKFORMAT='mgh_tf_negative', _STRICT_EXTRA=taxis_properties

   ;; If SHOW_ORIGIN keyword is set, create an appropriate symbol

   if keyword_set(show_origin) then $
        osym = ograph->NewSymbol(0, /FILL, NORM_SIZE=0.01, COLOR=mgh_color('dark green'))


   ;; Step through floats, extracting and plotting trajectories

   for f=0,n_elements(floats)-1 do begin

      ;; Get float vertical position data

      z = reform(oflt->VarGet('depth', OFFSET=[floats[f],0], COUNT=[1,0], $
                              AUTOSCALE=0))

      ;; To robustly establish if floats are unbounded, we need to use xgrid and ygrid data.

      xgrid = reform(oflt->VarGet('Xgrid', OFFSET=[floats[f],0], COUNT=[1,0], $
                                  AUTOSCALE=0))
      ygrid = reform(oflt->VarGet('Ygrid', OFFSET=[floats[f],0], COUNT=[1,0], $
                                  AUTOSCALE=0))

      l_bound = where(xgrid ge 0.5 and xgrid le dim_rho[0]-1.5 and $
                      ygrid ge 0.5 and ygrid le dim_rho[1]-1.5, n_bound, $
                      COMPLEMENT=l_unbound, NCOMPLEMENT=n_unbound)

      mgh_undefine, xgrid, ygrid

      if n_unbound gt 0 then z[l_unbound] = !values.f_nan

      ograph->NewAtom, 'IDLgrPlot', time, z, COLOR=mgh_color('red'), $
           _STRICT_EXTRA=plot_properties

      if keyword_set(show_origin) then begin
         r0 = min(l_bound)
         if r0 ge 0 then $
              ograph->NewAtom, 'IDLgrPlot', z[[r0]], time[[r0]], SYMBOL=osym
      endif

   endfor

   ma = ['Magnify','Translate','Context']

   ok = self->MGH_Window::Init(GRAPHICS_TREE=ograph, $
                               MOUSE_ACTION=ma,  _STRICT_EXTRA=extra)

   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Window'

   ;; Finalise plot & return

   self->Finalize, 'MGH_ROMS_Ztraj_Plot'

   return, 1

end

; MGH_ROMS_Ztraj_Plot::Cleanup
;
pro MGH_ROMS_Ztraj_Plot::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if self.float_destroy then obj_destroy, self.float_file
   if self.history_destroy then obj_destroy, self.history_file

   self->MGH_Window::Cleanup

end

; MGH_ROMS_Ztraj_Plot::GetProperty
;
pro MGH_ROMS_Ztraj_Plot::GetProperty, $
     ALL=all, FLOAT_FILE=float_file, HISTORY_FILE=history_file, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::GetProperty, ALL=all, _STRICT_EXTRA=extra

   float_file = self.float_file
   history_file = self.history_file

   if arg_present(all) then $
        all = create_struct(all, $
                            'float_file', float_file, 'history_file', history_file)

end

; MGH_ROMS_Ztraj_Plot::SetProperty
;
pro MGH_ROMS_Ztraj_Plot::SetProperty, _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::SetProperty, _STRICT_EXTRA=extra

end

; MGH_ROMS_Ztraj_Plot::About
;
;   Print information about the window and its contents
;
pro MGH_ROMS_Ztraj_Plot::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::About, lun

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

pro MGH_ROMS_Ztraj_Plot::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::BuildMenuBar

end


function MGH_ROMS_Ztraj_Plot::Event, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case widget_info(event.id, /UNAME) of

      else: return, self->MGH_Window::Event(event)

   endcase

end

pro MGH_ROMS_Ztraj_Plot::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Window::ExportData, values, labels

   self->GetProperty, $
        FLOAT_FILE=float_file, HISTORY_FILE=history_file

   labels = [labels, 'Float Object', 'History Object']
   values = [values, ptr_new(float_file), ptr_new(history_file)]

end

pro MGH_ROMS_Ztraj_Plot__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
        {MGH_ROMS_Ztraj_Plot, inherits MGH_Window, $
         float_file: obj_new(), float_destroy: 0B, $
         history_file: obj_new(), history_destroy: 0B}

end
