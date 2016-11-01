;+
; CLASS NAME:
;   MGH_ROMS_Movie_Cslice
;
; PURPOSE:
;   This procedure generates and displays a an animated sequence of
;   graphs showing a ROMS Cslice (in the sense described in
;   mghromshistory__define.pro)
;
; CATEGORY:
;   Ocean models.
;   Object graphics
;
; CALLING SEQUENCE:
;   mgh_new, 'MGH_ROMS_Movie_Cslice', File, Var
;
; ARGUMENTS:
;   File (input)
;     The name of a history or restart file produced by ROMS.
;
;   Var (input)
;     The name of a 3-D variable in the netCDF file.
;
; KEYWORD PARAMETERS:
;   ALONG_RANGE
;     Range of indices along the , relative to the rho grid.
;
;   DATA_MULTIPLIER (input, scalar numeric)
;     Number by which data values are multiplied before they are loaded into
;     the density surface. Default depends on the variable and is calculated by
;     MGH_ROMS_RESOLVE_DATA
;
;   DATA_RANGE (input, 2-element numeric)
;     Data range for the density surface. Default depends on the variable and
;     is calculated by MGH_ROMS_RESOLVE_DATA
;
;   DATA_TRANSFORMATION (input, numeric string)
;     Function applied to data values before they are loaded into the density
;     surface.
;
;   DEPTH_RANGE
;     Range for the vertical (depth) axis. Default is calculated from
;     the bathymetry and is [-1.02,0.02]*max(h)
;
;   DIRECTION
;     Set this equal to 0 for xi-z slice and 1 for a eta-z
;     slice. Default is 0.
;
;   SLICE
;     The position in the perpendicular direction. Default is nyp/2,
;     where nyp is the number of rho points in the perpendicular direction.
;
;   USE_BATH (input, switch)
;     Controls whether the vertical grid changes with time using bath (dynamic
;     bathymetry) data. Default is now to use data if available.
;
;   USE_ZETA (input, switch)
;     Controls whether the vertical grid changes with time using zeta data.
;     Default is now to use data if available.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 1998-09:
;     Written.
;   Mark Hadfield, 2000-09:
;     Modernised.
;   Mark Hadfield, 2001-10:
;     - The object is now named MGH_Roms_Cslice_Movie and uses the
;       C-slice functionality built into MGHromsHistory.
;     - Updated for IDL 5.5.
;   Mark Hadfield, 2002-01
;     * Reformatted with IDLWAVE.
;   Mark Hadfield, 2010-??:
;     Renamed MGH_Roms_Movie_Cslice
;   Mark Hadfield, 2010-12:
;     The default values for ALONG_RANGE were broken and have been
;     set to something more sensible.
;   Mark Hadfield, 2011-05:
;     Now supports a dynamic vertical grid affected by bathymetry (USE_BATH)
;     and sea surface height. Cleaned up the grid handling.
;   Mark Hadfield, 2011-07:
;     Removed history_destroy field: not necessary with automatic garbage
;     collection.
;   Mark Hadfield, 2011-11:
;     Default aspect ratio for the plot is now 0.5.
;   Mark Hadfield, 2013-02:
;     - Ported several items of functionality from mgh_roms_movie_hslice:
;       RECORD_AVERAGE, COLORBAR_PROPERTIES, DATA_TRANFORMATION, LOGARITHMIC.
;     - Renamed BAR_VISIBLE to SHOW_COLORBAR.
;-
function MGH_ROMS_Movie_Cslice::Init, history, variable, $
     DIRECTION=direction, INDEX=index, ALONG_RANGE=along_range, $
     DATA_MULTIPLIER=data_multiplier, $
     DATA_RANGE=data_range, LOGARITHMIC=logarithmic, $
     DATA_TRANSFORMATION=data_transformation, $
     DEPTH_RANGE=depth_range, $
     EXAGGERATE_ZETA=exaggerate_zeta, $
     RECORD_AVERAGE=record_average, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     RECORDS=records, $
     SHOW_COLORBAR=show_colorbar, $
     SHOW_CONTOUR=show_contour, $
     SHOW_MESH=show_mesh, SHOW_TIME=show_time, TITLE=title, $
     USE_BATH=use_bath, USE_ZETA=use_zeta, $
     COLORBAR_PROPERTIES=colorbar_properties, $
     CONTOUR_PROPERTIES=contour_properties, $
     GRAPH_PROPERTIES=graph_properties, $
     PALETTE_PROPERTIES=palette_properties, $
     _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process history argument

   if n_elements(history) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'history'

   case size(history, /TNAME) of
      'STRING': begin
         ohis = obj_new('MGHromsHistory', history)
         self.history_file = ohis
      end
      'OBJREF': begin
         ohis = history
         self.history_file = history
      end
      else: $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'history'
   endcase
   ohis = self.history_file

   ;; Check variable name

   if n_elements(variable) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'variable'

   if n_elements(variable) gt 1 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'variable'

   if size(variable, /TNAME) ne 'STRING' then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'variable'

   self.variable = variable

   ;; Other defaults

   mgh_roms_resolve_data, self.variable, $
        DATA_RANGE=data_range, DATA_MULTIPLIER=data_multiplier

   if n_elements(exaggerate_zeta) eq 0 then exaggerate_zeta = 1.

   if n_elements(show_colorbar) eq 0 then show_colorbar = 1B

   if n_elements(show_mesh) eq 0 then show_mesh = 0B

   if n_elements(use_bath) eq 0 then use_bath = ohis->HasVar('bath')

   if n_elements(use_zeta) eq 0 then use_zeta = ohis->HasVar('zeta')

   if n_elements(title) eq 0 then title = ''
   if n_elements(show_time) eq 0 then show_time = 1B

   ;; Determine the Cslice grid

   grid = ohis->CsliceGrid(DIRECTION=direction, INDEX=index, ALONG_RANGE=along_range, $
                           LONLAT=lonlat)
   self.grid = ptr_new(grid)

   n_arc = n_elements(grid.arc)

   ;; Establish records to be plotted (if applicable). Establish name of
   ;; time variable--the logic here is a bit shaky, as conventions
   ;; for defining time dimensions and variables in ROMS files are
   ;; not entirely consistent.

   if n_elements(record_average) eq 0 then record_average = 1

   dims = ohis->VarDims(variable)

   has_time = strlen(dims.time) gt 0

   if has_time then begin
      n_time = ohis->DimInfo(dims.time, /DIMSIZE)
      mgh_resolve_indices, n_time, record_range, record_stride, records
      n_records = n_elements(records)
      time_var = ohis->TimeVarName(dims.time)
      if isa(time_var, /NULL) then message, 'Time variable not found'
   endif else begin
      n_records = 1
      ;; Check other record keywords for validity here?
      if record_average gt 1 then $
         message, 'Cannot average over records for time-independent variable'
   endelse

   n_frames = long(n_records)/long(record_average)

   ;; Get time units from the history file

   if has_time then begin
      if ohis->HasAtt(time_var, 'units') then begin
         time_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
      endif else begin
         time_units = {scale: 1, offset: 0}
      endelse
   endif

   ;; Get (initial or static) bathymetry along the slice. This will be
   ;; recalculated for each record if USE_BATH is set.

   if use_bath then begin
      bath = ohis->CtranData('bath', GRID=grid, RECORD=records[0])
   endif else begin
      bath = grid.h
   endelse

   ;; Specify initial zeta. This will be recalculated for each record
   ;; if USE_ZETA is set.

   if use_zeta then begin
      zeta = ohis->CtranData('zeta', GRID=grid, RECORD=records[0])
   endif else begin
      zeta = fltarr(n_elements(grid.arc))
   endelse

   ;; Set default for DEPTH_RANGE

   if n_elements(depth_range) eq 0 then $
        depth_range = [-1.05,0.05]*max(bath)

   ;; Create base graph

   xmargin = show_colorbar ? [0.375,0.4] : [0.375,0.15]

   ograph = obj_new('MGHgrGraph2D', NAME='ROMS C-slice animation', $
                    ASPECT=0.5, XMaRGIN=xmargin, _STRICT_EXTRA=graph_properties)

   ograph->NewMask

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fsize

   ograph->NewFont, SIZE=fsize
   ograph->NewFont, SIZE=0.9*fsize

   ;; Draw axes

   ograph->NewAxis, 0, RANGE=mgh_minmax(grid.arc), /EXACT, $
        TITLE='Distance (km)', TICKFORMAT='mgh_tf_linear', $
        TICKFRMTDATA={scale: 1.E-3, format: '(F10.1)'}
   ograph->NewAxis, 1, RANGE=depth_range, /EXACT, $
        TITLE='Depth (m)', TICKFORMAT='mgh_tf_negative', $
        TICKFRMTDATA={format: '(F10.1)'}

   ;; Plot bathymetry. This will be animated if USE_BATH is set,

   xb = grid.arc#replicate(1, 2)
   yb = [[replicate(-1.05*max(bath), n_arc)],[-bath]]
   zb = make_array(DIMENSION=size(xb, /DIMENSIONS))
   ograph->NewAtom, 'IDLgrSurface', DATAX=xb, DATAY=yb, DATAZ=zb, $
        STYLE=2, COLOR=[127,127,127], RESULT=obath

   ;; Plot land mask. We have to inflate the mask vector into
   ;; an [n,2] array and the vertex arrasy to [n+1,3] to get the
   ;; color-plane code working.

   xm = mgh_stagger(grid.arc, DELTA=1)#replicate(1, 3)
   ym = [[replicate(-1.05*max(bath), n_arc+1)], $
         [replicate(-0.55*max(bath), n_arc+1)], $
         [replicate(0, n_arc+1)]]
   ograph->NewAtom, 'MGHgrColorPlane', RESULT=oland, $
        STYLE=0, DEFAULT_COLOR=[190,190,190], ZVALUE=0, $
        MISSING_POINTS=grid.mask#replicate(1, 2), $
        DATAX=xm, DATAY=ym, $
        NAME='Land mask'

   ;; Create the palette and add a colour bar

   ograph->NewPalette, 'MGH Special 6', RESULT=palette, $
        _STRICT_EXTRA=palette_properties

   ograph->NewColorBar, RESULT=obar, $
        FONT=ograph->GetFont(POS=1), $
        DATA_RANGE=data_range, HIDE=(~ show_colorbar), $
        LOGARITHMIC=logarithmic, PALETTE=palette, $
        CONTOUR_PROPERTIES=contour_properties, $
        SHOW_CONTOUR=show_contour, $
        _STRICT_EXTRA=colorbar_properties
   self.bar = obar

   ;; Add various graphics objects that will (or might) be animated...

   ;; ...title

   otitle = ograph->NewTitle(title)

   ;; ...density plane showing data values.  Position will be
   ;; recalculated for each record if USE_BATH or USE_ZETA is set.

   grid_z = ohis->CSliceZ(variable, GRID=grid, BATH=bath, ZETA=zeta)

   dim_z = size(grid_z, /DiMENSIONS)

   ograph->NewAtom, 'MGHgrDensityPlane', $
        STYLE=1, DEPTH_OFFSET=1, ZVALUE=-2*deltaz, $
        DATA_VALUES=fltarr(dim_z), $
        DATAX=cmreplicate(grid.arc, dim_z[1]), DATAY=grid_z, $
        DATA_RANGE=data_range, COLORSCALE=self.bar, RESULT=oplane
   self.plane = oplane

   ;; ...mesh or ruled surface showing levels.  Position will be
   ;; recalculated for each record if USE_BATH or USE_ZETA is set.

   ograph->NewAtom, 'IDLgrSurface', STYLE=3, COLOR=mgh_color('blue'), $
        DATAZ=fltarr(dim_z), $
        DATAX=cmreplicate(grid.arc, dim_z[1]), DATAY=grid_z, $
        HIDE=(1-show_mesh), RESULT=omesh
   self.mesh = omesh

   ;; ...plot showing surface height. If USE_ZETA is set,
   ;; this will be recalculated at each frame.

   ograph->NewAtom, 'IDLgrPlot', grid.arc, replicate(0, n_arc), RESULT=ozeta

   ;; ...contour showing data values.

   if keyword_set(show_contour) then begin
     ograph->NewAtom, 'IDLgrContour', RESULT=ocont, $
          /PLANAR, GEOMZ=-deltaz, $
          DATA=fltarr(dim_z), $
          GEOMX=cmreplicate(grid.arc, dim_z[1]), GEOMY=grid_z, $
          _STRICT_EXTRA=contour_properties
     self.contour = ocont
   endif

   ;; Create an animator window to display and manage the movie.

   mouse_action = ['Magnify','Translate','Context']
   ok = self->MGH_Datamator::Init(CHANGEABLE=0, GRAPHICS_TREE=ograph, $
                                  MOUSE_ACTION=mouse_action, $
                                  _STRICT_EXTRA=_extra)

   if ~ ok then $
        message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Datamator'

   ;; Step through times in the netCDF file, generating new frames &
   ;; plotting data

   oframe = objarr(6)

   ra = record_average
   dm = data_multiplier

   for f=0,n_frames-1 do begin

      if self->Finished() then break

      rec0 = ra*f

      time = 0
      slice = 0

      if use_bath then bath = 0
      if use_zeta then zeta = 0

      for r=rec0,rec0+ra-1 do begin

         if has_time then begin

            time += ohis->VarGet(time_var, OFFSET=records[r], $
                                 COUNT=[1], AUTOSCALE=0)*time_units.scale

            slice += ohis->CSliceData(variable, GRID=grid, RECORD=records[r])

            if use_bath then begin
               bath += ohis->CtranData('bath', GRID=grid, RECORD=records[r])
            endif

            if use_zeta then begin
               zeta += ohis->CtranData('zeta', GRID=grid, RECORD=records[r])
            endif

         endif else begin

            slice += ohis->CsliceData(variable, GRID=grid)

         endelse

      endfor

      if has_time then time = time/double(ra)

      slice = dm*slice/float(ra)

      if n_elements(data_transformation) gt 0 then $
           slice = call_function(data_transformation, slice)

      if use_zeta then begin
         zeta = zeta*exaggerate_zeta/float(ra)
      endif

      if use_bath then begin
         bath = bath/float(ra)
      endif

      if use_bath || use_zeta then $
           grid_z = ohis->CSliceZ(variable, GRID=grid, BATH=bath, ZETA=zeta)

      ;; Update title

      if has_time && show_time gt 0 then begin
         case show_time of
            1: ttt = string(FORMAT='(%"%0.3f days")', time)
            2: ttt = mgh_dt_string(time+time_units.offset)
         endcase
         if strlen(title) gt 0 then $
              ttt = string(FORMAT='(%"%s: %s")', title, ttt)
         oframe[0] = obj_new('MGH_Command', OBJECT=otitle, $
                             'SetProperty', STRINGS=temporary(ttt))
      endif

      ;; Set properties for remaining animated atoms. Note that some
      ;; of the GUI event-handling code relies on the order of the
      ;; command objects in the oframe array.

      ;; Note that setting DATAY alone for the plane object seems to
      ;; discard any existing DATA_VALUES. I'm not sure that this
      ;; should be happening, but the following code works around this.

      if use_bath || use_zeta then begin
         oframe[1] = obj_new('MGH_Command', OBJECT=oplane, $
                             'SetProperty', DATA_VALUES=slice, $
                             DATAY=grid_z)
         oframe[2] = obj_new('MGH_Command', OBJECT=omesh, $
                             'SetProperty', DATAY=grid_z)
      endif else begin
         oframe[1] = obj_new('MGH_Command', OBJECT=oplane, $
                             'SetProperty', DATA_VALUES=slice)
      endelse

      if use_zeta then begin
         oframe[3] = obj_new('MGH_Command', OBJECT=ozeta, $
                             'SetProperty', DATAY=zeta)
      endif

      if use_bath then begin
         yb = [[mgh_reproduce(-1.05*max(bath), bath)],[-bath]]
         oframe[4] = obj_new('MGH_Command', OBJECT=obath, $
                             'SetProperty', DATAY=yb)
      endif

      if keyword_set(show_contour) then begin
        oframe[5] = obj_new('MGH_Command', OBJECT=self.contour, $
                            'SetProperty', DATA=slice)
      endif

      ;; Add frame to animator & display

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; MGH_ROMS_Movie_Cslice::Cleanup
;
pro MGH_ROMS_Movie_Cslice::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.grid

   self->MGH_Datamator::Cleanup

end

; MGH_ROMS_Movie_Cslice::GetProperty
;
pro MGH_ROMS_Movie_Cslice::GetProperty, $
     ALL=all, BAR=bar, BYTE_RANGE=byte_range, DATA_RANGE=data_range, GRID=grid, $
     HISTORY_FILE=history_file, PALETTE=palette, SHOW_MESH=show_mesh, VARIABLE=variable, $
     _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::GetProperty, ALL=all, _STRICT_EXTRA=_extra

   self.mesh->GetProperty, HIDE=show_mesh
   show_mesh = 1 - show_mesh

   bar = self.bar

   if arg_present(all) or arg_present(grid) then grid = *self.grid

   history_file = self.history_file

   variable = self.variable

   if obj_valid(self.bar) then begin
      self.bar->GetProperty, BYTE_RANGE=byte_range, $
           DATA_RANGE=data_range, PALETTE=palette
   endif else begin
      byte_range = [0B,0B]
      data_range = [0,0]
      palette = obj_new()
   endelse

   if arg_present(all) then $
        all = create_struct(all, 'byte_range', byte_range, $
                            'data_range', data_range, 'grid=', grid, $
                            'history', history, 'palette', palette, $
                            'show_mesh', show_mesh, 'variable', variable)

end

; MGH_ROMS_Movie_Cslice::SetProperty
;
pro MGH_ROMS_Movie_Cslice::SetProperty, $
     BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
     SHOW_MESH=show_mesh, _REF_EXTRA=_extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(byte_range) gt 0 || n_elements(data_range) gt 0 then begin

      if n_elements(byte_range) eq 0 then $
           self->GetProperty, BYTE_RANGE=byte_range
      if n_elements(data_range) eq 0 then $
           self->GetProperty, DATA_RANGE=data_range

      if obj_valid(self.bar) then $
           self.bar->SetProperty, BYTE_RANGE=byte_range, DATA_RANGE=data_range

      self.plane->SetProperty, BYTE_RANGE=byte_range, DATA_RANGE=data_range

   endif

   if n_elements(show_mesh) gt 0 then $
        self.mesh->SetProperty, HIDE=(1-show_mesh)

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=_extra

end

; MGH_ROMS_Movie_Cslice::About
;
;   Print information about the window and its contents
;
pro MGH_ROMS_Movie_Cslice::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::About, lun

   self->GetProperty, BAR=bar, GRID=grid, HISTORY=history, $
     PALETTE=palette, VARIABLE=variable

   if obj_valid(history) then begin
      printf, lun, self, ': the history file sequence is ', $
              mgh_obj_string(history)
      history->GetProperty, FILE_NAME=file_name
      printf, lun, self, ': the files are:', file_name
   endif

   printf, lun, self, ': the variable name is '+variable

   printf, lun, self, ': the slice direction is ', strtrim(grid.direction,2)

   printf, lun, self, ': the slice index and range are ', $
           grid.slice, grid.along_range

   if obj_valid(self.plane) then begin
      printf, lun, self, ': the density plane object is ', $
              mgh_obj_string(self.plane, /SHOW_NAME)
   endif

   if obj_valid(bar) then begin
      printf, lun, self, ': the colour bar object is ', $
              mgh_obj_string(self.bar, /SHOW_NAME)
   endif

   if obj_valid(palette) then begin
      printf, lun, self, ': the palette is ', $
              mgh_obj_string(palette, /SHOW_NAME)
   endif

end

; MGH_ROMS_Movie_Cslice::BuildMenuBar

pro MGH_ROMS_Movie_Cslice::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::BuildMenuBar

   obar = mgh_widget_self(self.menu_bar)

   obar->NewItem, PARENT='Tools', $
        SEPARATOR=[1,0,0,0,1,1], MENU=[1,0,0,0,0,0], $
        CHECKED_MENU=[0,0,0,0,1,0], $
        ['Data Range','Edit Palette...','View Colour Scale...', $
         'View Data Values...','Mesh','Plot Grid']

   obar->NewItem, PARENT='Tools.Data Range', ['Set...','Fit this Frame']

end


; MGH_ROMS_Movie_Cslice::EventMenuBar
;
function MGH_ROMS_Movie_Cslice::EventMenuBar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      'TOOLS.DATA RANGE.SET': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Range', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, IMMEDIATE=0, $
                  N_ELEMENTS=2, PROPERTY_NAME='DATA_RANGE'
         return, 0
      end

      'TOOLS.DATA RANGE.FIT THIS FRAME': begin
         self->GetProperty, POSITION=position
         oframe = self.animation->GetFrame(POSITION=position)
         oframe[1]->GetProperty, KEYWORDS=keywords
         data_range = mgh_minmax(keywords.data_values, /NAN)
         if data_range[0] eq data_range[1] then data_range += [-1,1]
         self->SetProperty, DATA_RANGE=data_range
         self->Update
         return, 0
      end

      'TOOLS.EDIT PALETTE': begin
         self->GetProperty, PALETTE=palette
         mgh_new, 'MGH_GUI_Palette_Editor', palette, CLIENT=self, /FLOATING, $
                  GROUP_LEADER=self.base, /IMMEDIATE
         return, 0
      end

      'TOOLS.VIEW COLOUR SCALE': begin
         mgh_new, 'MGH_GUI_ColorScale', CLIENT=self, /FLOATING, $
                  GROUP_LEADER=self.base
         return, 0
      end

      'TOOLS.VIEW DATA VALUES': begin
         self->GetProperty, POSITION=position
         oframe = self.animation->GetFrame(POSITION=position)
         oframe[1]->GetProperty, KEYWORDS=keywords
         data_dims = size(keywords.data_values, /DIMENSIONS)
         xvaredit, keywords.data_values, GROUP=self.base, $
                   X_SCROLL_SIZE=(data_dims[0] < 12), $
                   Y_SCROLL_SIZE=(data_dims[1] < 30)
         return, 0
      end

      'TOOLS.MESH': begin
         self->GetProperty, SHOW_MESH=show_mesh
         self->SetProperty, SHOW_MESH=(1-show_mesh)
         self->Update
         return, 1
      end

      'TOOLS.PLOT GRID': begin
         self->ShowGridPlot
         return, 1
      end

      else: return, self->MGH_Datamator::EventMenuBar(event)

   endcase

end

; MGH_ROMS_Movie_Cslice::ExportData
;
pro MGH_ROMS_Movie_Cslice::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Player::ExportData, values, labels

   self->GetProperty, ANIMATION=animation, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[1]->GetProperty, KEYWORDS=keywords

   labels = [labels, 'Slice Data']
   values = [values, ptr_new(keywords.data_values)]

end

; MGH_ROMS_Movie_Cslice::PickReport
;
pro MGH_ROMS_Movie_Cslice::PickReport, pos, LUN=lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(lun) eq 0 then lun = -1

   if n_elements(pos) ne 2 then $
        message, 'Parameter POS must be a 2-element vector'

   self->GetProperty, GRAPHICS_TREE=ograph

   if ~ obj_valid(ograph) then begin
      printf, lun, FORMAT='(%"%s: no graphics tree")', mgh_obj_string(self)
      return
   endif

   printf, lun, FORMAT='(%"%s: graphics tree %s")', $
           mgh_obj_string(self), mgh_obj_string(ograph, /SHOW_NAME)

   atoms = self->Select(ograph, pos)
   valid = where(obj_valid(atoms), n_atoms)

   if n_atoms eq 0 then begin
      printf, lun, FORMAT='(%"%s: no atoms selected")', mgh_obj_string(self)
      return
   endif

   atoms = atoms[valid]

   for j=0,n_atoms-1 do begin
      atom = atoms[j]
      status = self->PickData(ograph, atom, pos, data)
      case (atom eq self.plane) of
         0: begin
            printf, lun, FORMAT='(%"%s: atom %s, success: %d, value: %f %f %f")', $
                    mgh_obj_string(self), mgh_obj_string(atom,/SHOW_NAME), $
                    status, double(data)
         end
         1: begin
            ;; If the selected atom is the density plane, report the
            ;; data value at the selected location.

            ;; Locate the selection point in the index space of the
            ;; density planes' pixel vertices.
            self.plane->GetProperty, $
                 DATAX=datax, DATAY=datay, DATA_VALUES=data_values, STYLE=style
            xy2d = size(datax, /N_DIMENSIONS) eq 2
            case xy2d of
               0: begin
                  ii = mgh_locate(datax, XOUT=data[0])
                  jj = mgh_locate(datay, XOUT=data[1])
               end
               1: begin
                  loc = mgh_locate2a(datax, datay, XOUT=data[[0]], YOUT=data[[1]], $
                                     MISSING=-1)
                  ii = reform(loc[0,*,*])
                  jj = reform(loc[1,*,*])
                  mgh_undefine, loc
               end
            endcase

            ;; If style is 0, allow for offset of data locations
            ;; (pixel centres) and use nearest-neighbour interpolation
            if style eq 0 then begin
               ii = round(ii-0.5)
               jj = round(jj-0.5)
            endif
            ;; Interpolate & report
            value = mgh_interpolate(data_values, ii, jj, GRID=(xy2d eq 0), $
                                    MISSING=!values.f_nan)
            printf, lun, FORMAT='(%"%s: atom %s, success: %d, value: %f %f %f")', $
                    mgh_obj_string(self), mgh_obj_string(atom,/SHOW_NAME), $
                    status, double(data)
            printf, lun, FORMAT='(%"%s: atom %s, index: %f %f, value: %f")', $
                    mgh_obj_string(self), mgh_obj_string(atom,/SHOW_NAME), $
                    ii,jj,value
         end
      endcase
   endfor

end

; MGH_ROMS_Movie_Cslice::ShowGridPlot
;
; Purpose:
;   Generate a grph showing the slice position(s)

pro MGH_ROMS_Movie_Cslice::ShowGridPlot

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   mgh_graph_default, SCALE=default_scale

   self->GetProperty, GRID=grid

   ohis = mgh_obj_clone(self.history_file)

   oplot = obj_new('mgh_roms_plot_grid', ohis, GRAPH_PROPERTIES={scale: 0.7*default_scale}, $
                   GROUP_LEADER=self.base, VISIBLE=0)

   oplot->GetProperty, GRAPHICS_TREE=ograph

   ograph->NewAtom, 'IDLgrPlot', grid.x, grid.y, COLOR=mgh_color('white')

   oplot->Update

   oplot->SetProperty, /VISIBLE

end


; MGH_ROMS_Movie_Cslice::UpdateMenuBar
;
pro MGH_ROMS_Movie_Cslice::UpdateMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::UpdateMenuBar

   obar = mgh_widget_self(self.menu_bar)

   if obj_valid(obar) then begin
      self->GetProperty, SHOW_MESH=show_mesh
      obar->SetItem, 'Tools.Mesh', SET_BUTTON=show_mesh
   endif

end

pro MGH_ROMS_Movie_Cslice__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, {MGH_ROMS_Movie_Cslice, inherits MGH_Datamator, $
                 history_file: obj_new(), variable: '', $
                 direction: 0S, grid: ptr_new(), bar: obj_new(), $
                 plane: obj_new(), mesh: obj_new(), $
                 contour: obj_new()}

end
