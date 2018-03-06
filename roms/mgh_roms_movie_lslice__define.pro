;+
; CLASS NAME:
;  MGH_ROMS_Movie_Lslice
;
; PURPOSE:
;   This class generates and displays an animated sequence of graphs
;   showing an Lslice (lon-lat slice) through a ROMS 2D forcing field.
;
; CALLING SEQUENCE:
;   mgh_new, 'mgh_roms_movie_lslice', history, variable
;
; POSITIONAL PARAMETERS:
;   history (input)
;     A reference to a ROMS history sequence object or a string arrary
;     specifying a list of ROMS history files or a single string with
;     wildcards specifying a list of ROMS history files.
;
;   variable
;     The name of a 2-D or 3-D variable in the netCDF file.
;
; KEYWORD PARAMETERS:
;   DATA_MULTIPLIER (input, numeric scalar)
;     Number by which data values are multiplied before they are loaded into
;     the density surface. Default depends on the variable and is calculated by
;     MGH_ROMS_RESOLVE_DATA
;
;   DATA_RANGE (input, numeric 2-element vector)
;     Data range for the density surface. Default depends on the variable and
;     is calculated by MGH_ROMS_RESOLVE_DATA
;
;   DATA_TRANSFORMATION (input, numeric string)
;     Function applied to data values before they are loaded into the density
;     surface.
;
;   MAP_STRUCTURE
;     If a map structure is supplied via this keyword, and if there
;     are longitude and latitude data in the file, then the longitudes
;     and latitudes are converted with MAP_PROJ_FORWARD before display.
;
;   SHOW_TIME (input, integer)
;     Show the time or date in the title of each frame.
;
;   I_RANGE
;   J_RANGE
;     Use these keywords to display a subset of the domain in the lon (I)
;     and lat (J) directions.
;
;###########################################################################
; Copyright (c) 2018 NIWA:
;   http://www.niwa.co.nz/
; Licensed under the MIT open source license:
;   http://www.opensource.org/licenses/mit-license.php
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2018-03:
;     Written.
;-
function mgh_roms_movie_lslice::Init, $
     history, variable, $
     BYTE_RANGE=byte_range, $
     DATA_MULTIPLIER=data_multiplier, $
     DATA_RANGE=data_range, LOGARITHMIC=logarithmic, $
     DATA_TRANSFORMATION=data_transformation, $
     MAP_STRUCTURE=map_structure, $
     RECORD_AVERAGE=record_average, $
     RECORD_RANGE=record_range, $
     RECORD_STRIDE=record_stride, $
     RECORDS=records, $
     STYLE=style, $
     SHOW_COLORBAR=show_colorbar, SHOW_CONTOUR=show_contour, $
     SHOW_TITLE=show_title, SHOW_TIME=show_time, $
     TITLE=title, DT_FORMaT=dt_format, $
     I_RANGE=i_range, J_RANGE=j_range, $
     X_RANGE=x_range, Y_RANGE=y_range, $
     COLORBAR_PROPERTIES=colorbar_properties, $
     CONTOUR_PROPERTIES=contour_properties, $
     GRAPH_PROPERTIES=graph_properties, $
     LAND_PROPERTIES=land_properties, $
     PALETTE_PROPERTIES=palette_properties, $
     XAXIS_PROPERTIES=xaxis_properties, $
     YAXIS_PROPERTIES=yaxis_properties, $
     _REF_EXTRA=extra

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

   ;; Check variable name

   if n_elements(variable) eq 0 then variable = 'shflux'

   if n_elements(variable) gt 1 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'variable'

   if size(variable, /TNAME) ne 'STRING' then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'variable'

   self.variable = variable

   ;; Other defaults

   mgh_roms_resolve_data, self.variable, $
        DATA_RANGE=data_range, DATA_MULTIPLIER=data_multiplier

   if n_elements(style) eq 0 then style = 0

   if n_elements(show_colorbar) eq 0 then show_colorbar = !true

   if n_elements(show_contour) eq 0 then show_contour = !false

   if n_elements(map_structure) gt 0 then self.map_structure = ptr_new(map_structure)

   use_map_structure = ptr_valid(self.map_structure)

   if n_elements(show_title) eq 0 then show_title = !true

   if keyword_set(show_title) then begin
      if n_elements(title) eq 0 then title = ''
      if n_elements(show_time) eq 0 then show_time = 1
      if n_elements(dt_format) eq 0 then $
         dt_format = '(C(CYI4.4,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))'
   endif

   ;; Set I_RANGE and J_RANGE. The default is to show all points.

   dim = [ohis->DimInfo('lon', /DIMSIZE),ohis->DimInfo('lat', /DIMSIZE)]

   if n_elements(i_range) eq 0 then i_range = [0,dim[0]-1]
   if n_elements(j_range) eq 0 then j_range = [0,dim[1]-1]

   if i_range[0] lt 0 then i_range[0] += dim[0]
   if i_range[1] lt 0 then i_range[1] += dim[0]

   if j_range[0] lt 0 then j_range[0] += dim[1]
   if j_range[1] lt 0 then j_range[1] += dim[1]

   ;; Get grid data required for horizontal slice retrievals

   grid = ohis->LsliceGrid(I_RANGE=i_range, J_RANGE=j_range)

   ;; Convert 2D lon-lat data and convert to to map projection if appropriate

   x = mgh_inflate(grid.dim, grid.lon, 1)
   y = mgh_inflate(grid.dim, grid.lat, 2)

   if use_map_structure then begin
      xy = map_proj_forward(x, y, MAP_STRUCTURE=*self.map_structure)
      x = reform(xy[0,*], grid.dim)
      y = reform(xy[1,*], grid.dim)
      mgh_undefine, xy
   endif else begin
   endelse

   ;; Establish records to be plotted (if applicable). Establish name of
   ;; time variable--the logic here is a bit shaky, as conventions
   ;; for defining time dimensions and variables in ROMS files are
   ;; not entirely consistent. The following does not allow variables
   ;; without a time dimension.

   if n_elements(record_average) eq 0 then record_average = 1

   case !true of
      ohis->HasAtt(variable, 'time'): $
         time_var = ohis->AttGet(variable,'time')
      else: $
         time_var = 'ocean_time'
   endcase

   time_dim = (ohis->VarDimNames(time_var))[0]

   n_time = ohis->DimInfo(time_dim, /DIMSIZE)
   mgh_resolve_indices, n_time, record_range, record_stride, records
   n_records = n_elements(records)

   n_frames = long(n_records)/long(record_average)

   ;; Get time units from the history file. Also retrieve the time data

   if ohis->HasAtt(time_var, 'units') then begin
      time_units = mgh_dt_units(ohis->AttGet(time_var, 'units'))
   endif else begin
      time_units = {scale: 1, offset: 0}
   endelse

   time = ohis->VarGet(time_var, AUTOSCALE=0)*time_units.scale
   time = time[records]

   ;; Default graph aspect ratio, can be overridden via GRAPH_PROPERTIES

   if n_elements(x_range) eq 0 then x_range = mgh_minmax(x)
   if n_elements(y_range) eq 0 then y_range = mgh_minmax(y)

   aspect = mgh_aspect(x_range, y_range, LONLAT=(~ use_map_structure))
   aspect = (aspect > 0.4) < 1.5

   ;; Create base graph

   xmargin = show_colorbar ? [0.375,0.4] : [0.375,0.15]

   ograph = obj_new('MGHgrGraph2D', ASPECT=aspect, XMARGIN=xmargin, $
                    NAME='ROMS lon-lat slice animation', $
                    _STRICT_EXTRA=graph_properties)

   ograph->GetProperty, DELTAZ=deltaz, FONTSIZE=fontsize, PLOT_RECT=prect

   ograph->NewFont
   ograph->NewFont, SIZE=0.9*fontsize

   ;; Add mask around plot area

   ograph->NewMask

   ;; Draw axes

   if (~ use_map_structure) then begin
      xap = {tickformat: 'mgh_tf_longitude', tickfrmtdata: {format:'(F10.1)'}}
      yap = {tickformat: 'mgh_tf_latitude', tickfrmtdata: {format:'(F10.1)'} }
   endif else begin
      xap = {title: 'X (km)', tickformat: 'mgh_tf_linear', tickfrmtdata: {scale:1.E-3, format:'(F10.1)'}}
      yap = {title: 'Y (km)', tickformat: 'mgh_tf_linear', tickfrmtdata: {scale:1.E-3, format:'(F10.1)'}}
   endelse

   ograph->NewAxis, 0, $
        RANGE=x_range, /EXACT, $
        _STRICT_EXTRA=mgh_struct_merge(xap, xaxis_properties)
   ograph->NewAxis, 1, $
        RANGE=y_range, /EXACT, $
        _STRICT_EXTRA=mgh_struct_merge(yap, yaxis_properties)

   ;; Create the palette and add a colour bar

   ograph->NewPalette, 'Matlab Jet', RESULT=palette, $
        _STRICT_EXTRA=palette_properties

   ograph->NewColorBar, RESULT=obar, FONT=ograph->GetFont(POS=1), $
        DATA_RANGE=data_range, HIDE=(~ show_colorbar), $
        LOGARITHMIC=logarithmic, PALETTE=palette, $
        CONTOUR_PROPERTIES=contour_properties, $
        SHOW_CONTOUR=show_contour, $
        _STRICT_EXTRA=colorbar_properties
   self.bar = obar

   ;; Add various graphics objects that will (or might) be animated...

   ;; ...title

   if keyword_set(show_title) then otitle = ograph->NewTitle(title)

   ;; ...density plane showing data values.

   ograph->NewAtom, 'MGHgrDensityPlane', RESULT=oplane, STYLE=style, ZVALUE=-5*deltaz , $
        DATAX=mgh_stagger(x, DELTA=(style eq 0)), $
        DATAY=mgh_stagger(y, DELTA=(style eq 0)), $
        DATA_VALUES=fltarr(grid.dim), $
        COLORSCALE=self.bar, /STORE_DATA
   self.plane = oplane

   ;; ...contour showing data values.

   if show_contour then begin
      ograph->NewAtom, 'IDLgrContour', RESULT=ocont, /PLANAR, GEOMZ=-2*deltaz , $
           GEOMX=x, GEOMY=y, DATA=fltarr(grid.dim), $
           _STRICT_EXTRA=contour_properties
      self.contour = ocont
   endif

   ;; Create an animator window to display and manage the movie.

   ma = ['Magnify','Translate','Context']
   ok = self->MGH_Datamator::Init(CHANGEABLE=!false, GRAPHICS_TREE=ograph, MOUSE_ACTION=ma, _STRICT_EXTRA=extra)
   if ~ ok then message, BLOCK='MGH_MBLK_MOTLEY', NAME='MGH_M_INITFAIL', 'MGH_Datamator'

   ;; Step through the netCDF file, generating new frames & plotting data

   oframe = objarr(4)

   ra = record_average
   dm = data_multiplier

   for f=0,n_frames-1 do begin

      if self->Finished() then break

      rec0 = ra*f

      frame_time = 0
      frame_slice = 0

      frame_wind_uv = 0

      for r=rec0,rec0+ra-1 do begin

         frame_time += time[r]/double(ra)
         frame_slice += ohis->LsliceData(variable, GRID=grid, RECORD=records[r])

      endfor

      frame_slice = dm*frame_slice/float(ra)

      if n_elements(data_transformation) gt 0 then $
           frame_slice = call_function(data_transformation, frame_slice)

      ;; Note that if the LOGARITHMIC keyword is set, the logarithmic transformation
      ;; is handled within the methods of the MGHgrDensityPlane class.

      oframe[0] = obj_new('MGH_Command', OBJECT=self.plane, 'SetProperty', DATA_VALUES=frame_slice)
      if show_contour then $
           oframe[1] = obj_new('MGH_Command', OBJECT=self.contour, 'SetProperty', DATA=frame_slice)
     mgh_undefine, frame_slice

      ;; Update title with time or date, if applicable

      if keyword_set(show_title) && show_time gt 0 then begin
         case show_time of
            1: ttt = string(FORMAT='(%"%0.3f days")', frame_time)
            2: ttt = string(FORMAT='(%"%s")', mgh_dt_string(frame_time+time_units.offset, FORMAT=dt_format))
            3: ttt = string(FORMAT='(%"%s (%0.3f days)")', mgh_dt_string(frame_time+time_units.offset, FORMAT=dt_format), frame_time)
         endcase
         if strlen(title) gt 0 then $
            ttt = string(FORMAT='(%"%s: %s")', title, ttt)
         oframe[2] = obj_new('MGH_Command', OBJECT=otitle, 'SetProperty', STRINGS=temporary(ttt))
      endif

      if keyword_set(show_nzlam_wind) then begin
         oframe[3] = obj_new('MGH_Command', OBJECT=owind, 'SetProperty', DATAU=real_part(frame_wind_uv), DATAV=imaginary(frame_wind_uv))
      endif

      self->AddFrame, oframe

   endfor

   self->Finish

   return, 1

end

; Mgh_Roms_Movie_lslice::Cleanup
;
pro Mgh_Roms_Movie_lslice::Cleanup

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ptr_free, self.map_structure

   self->MGH_Datamator::Cleanup

end

; Mgh_Roms_Movie_lslice::GetProperty
;
pro Mgh_Roms_Movie_lslice::GetProperty, $
     ALL=all, BAR=bar, BYTE_RANGE=byte_range, DATA_RANGE=data_range, $
     HISTORY_FILE=history_file, LONLAT=lonlat, MAP_STRUCTURE=map_structure, $
     PALETTE=palette, STYLE=style, VARIABLE=variable, $
     _REF_EXTRA=extra

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::GetProperty, ALL=all, _STRICT_EXTRA=extra

   bar = self.bar

   history_file = self.history_file

   if arg_present(all) || arg_present(map_structure) then $
        map_structure = ptr_valid(self.map_structure) ? *self.map_structure : -1

   variable = self.variable

   self.plane->GetProperty, STYLE=style

   self.bar->GetProperty, $
        BYTE_RANGE=byte_range, DATA_RANGE=data_range, PALETTE=palette

   if arg_present(all) then $
        all = create_struct(all, $
                            'bar', bar, 'byte_range', byte_range, $
                            'data_range', data_range, $
                            'history_file', history_file, $
                            'lonlat', lonlat, 'map_structure', map_structure, $
                            'palette', palette, 'style', style, $
                            'variable', variable)

end

; Mgh_Roms_Movie_lslice::SetProperty
;
pro Mgh_Roms_Movie_lslice::SetProperty, $
     BYTE_RANGE=byte_range, DATA_RANGE=data_range, STYLE=style, _REF_EXTRA=extra

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

   self.plane->SetProperty, STYLE=style

   self->MGH_Datamator::SetProperty, _STRICT_EXTRA=extra

end

; Mgh_Roms_Movie_lslice::About
;
pro Mgh_Roms_Movie_lslice::About, lun

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::About, lun

   self->GetProperty, BAR=bar, HISTORY_FILE=history_file, $
        PALETTE=palette, VARIABLE=variable

   if obj_valid(history_file) then begin
      printf, lun, FORMAT='(%"%s: my history file sequence is %s")', $
           mgh_obj_string(self), mgh_obj_string(history_file)
      history_file->GetProperty, FILE_NAME=file_name
      printf, lun, self, ': its files are:', file_name
   endif

   printf, lun, FORMAT='(%"%s: my variable name is %s")', $
        mgh_obj_string(self), variable

   if obj_valid(self.bar) then begin
      printf, lun, self, ': the colour bar is ' , $
              mgh_obj_string(self.bar, /SHOW_NAME)
   endif

   if obj_valid(self.plane) then begin
      printf, lun, self, ': the density plane object is ', $
              mgh_obj_string(self.plane, /SHOW_NAME)
   endif

   if obj_valid(palette) then begin
      printf, lun, self, ': the palette is ', $
              mgh_obj_string(palette, /SHOW_NAME)
   endif

   if obj_valid(self.contour) then begin
      printf, lun, self, ': the contour is ' , $
              mgh_obj_string(self.contour, /SHOW_NAME)
   endif

end

; Mgh_Roms_Movie_lslice::BuildMenuBar
;
; Purpose:
;   Add menus, sub-menus & menu items to the menu bar

pro Mgh_Roms_Movie_lslice::BuildMenuBar

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

    self->MGH_Datamator::BuildMenuBar

    obar = mgh_widget_self(self.menu_bar)

    obar->NewItem, PARENT='File.Export Animation', 'NetCDF...'

   obar->NewItem, PARENT='Tools', SEPARATOR=[1,0,0,0,1], MENU=[1,0,0,0,0], $
        ['Data Range','Edit Palette...','Set Style...', $
         'View Colour Scale...','View Data Values...']

   obar->NewItem, PARENT='Tools.Data Range', ['Set...','Fit this Frame']

end


; Mgh_Roms_Movie_lslice::EventMenuBar
;
function Mgh_Roms_Movie_lslice::EventMenuBar, event

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   case event.value of

      'FILE.EXPORT ANIMATION.NETCDF': begin
         self.graphics_tree->GetProperty, NAME=name
         ext = '.nc'
         default_file = strlen(name) gt 0 ? mgh_str_vanilla(name)+ext : ''
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
            self->ExportAnimationDataToNcFile, filename
         endif
         return, 0
      end

      'TOOLS.DATA RANGE.SET': begin
         mgh_new, 'MGH_GUI_SetArray', CAPTION='Range', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, IMMEDIATE=0, $
                  N_ELEMENTS=2, PROPERTY_NAME='DATA_RANGE'
         return, 0
      end

      'TOOLS.DATA RANGE.FIT THIS FRAME': begin
         self->GetProperty, POSITION=position
         oframe = self.animation->GetFrame(POSITION=position)
         oframe[0]->GetProperty, KEYWORDS=keywords
         data_range = mgh_minmax(keywords.data_values, /NAN)
         if data_range[0] eq data_range[1] then data_range += [-1,1]
         self->SetProperty, DATA_RANGE=data_range
         self->Update
         return, 0
      end

      'TOOLS.EDIT PALETTE': begin
         self->GetProperty, PALETTE=palette
         mgh_new, 'MGH_GUI_Palette_Editor', palette, CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, /IMMEDIATE
         return, 0
      end

      'TOOLS.SET STYLE': begin
         mgh_new, 'MGH_GUI_SetList', CAPTION='Style', CLIENT=self, $
                  /FLOATING, GROUP_LEADER=self.base, IMMEDIATE=0, $
                  ITEM_STRING=['Block','Interpolated'], $
                  PROPERTY_NAME='STYLE'
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
         oframe[0]->GetProperty, KEYWORDS=keywords
         data_dims = size(keywords.data_values, /DIMENSIONS)
         ;; Call REFORM so that XVAREDIT cannot modify values
         xvaredit, reform(keywords.data_values), GROUP=self.base, $
                   X_SCROLL_SIZE=(data_dims[0] < 8), $
                   Y_SCROLL_SIZE=(data_dims[1] < 30)
         return, 0
      end

      else: return, self->MGH_Datamator::EventMenuBar(event)

   endcase

end

; Mgh_Roms_Movie_lslice::ExportData
;
pro Mgh_Roms_Movie_lslice::ExportData, values, labels

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   self->MGH_Datamator::ExportData, values, labels

   self->GetProperty, ANIMATION=animation, POSITION=position

   oframe = animation->GetFrame(POSITION=position)
   oframe[0]->GetProperty, KEYWORDS=keywords

   labels = [labels, 'Slice Data']
   values = [values, ptr_new(keywords.data_values)]

end

; Mgh_Roms_Movie_lslice::PickReport
;
pro Mgh_Roms_Movie_lslice::PickReport, pos, LUN=lun

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

; Mgh_Roms_Movie_lslice::ExportAnimationDataToNcFile
;
pro Mgh_Roms_Movie_lslice::ExportAnimationDataToNcFile, file

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Query the appropriate objects for info

   self->GetProperty, $
      LONLAT=lonlat, MAP_STRUCTURE=map_structure, VARIABLE=variable

   self.animation->GetProperty, $
      GRAPHICS_TREE=ograph, N_FRAMES=n_frames

   self.plane->GetProperty, $
      DATAX=datax, DATAY=datay, STYLE=style

   self.animator->GetPlayBack, $
      RANGE=play_range, USE_RANGE=play_use_range

   ;; Sort out grid data

   if style eq 0 then begin
      datax = mgh_stagger(temporary(datax), DELTA=[-1,-1])
      datay = mgh_stagger(temporary(datay), DELTA=[-1,-1])
   endif

   dim = size(datax, /DIMENSIONS)

   if lonlat && size(map_structure, /TYPE) eq 8 then begin
      ll = map_proj_inverse(datax, datay, MAP_STRUCTURE=map_structure)
      datax = (reform(ll[0,*], dim)+360) mod 360
      datay = reform(ll[1,*], dim)
   endif

   ;; Sort out frames to be read

   if play_use_range then begin
      if n_elements(range) eq 0 then range = play_range[0:1]
      if n_elements(stride) eq 0 then stride = play_range[2]
   endif else begin
      if n_elements(range) eq 0 then range = [0,n_frames-1]
      if n_elements(stride) eq 0 then stride = 1
   endelse

   ;; Set up netCDF file

   onc = obj_new('MGHncFile', file, /CREATE, /CLOBBER)

   onc->AttAdd, /GLOBAL, 'title', 'ROMS lslice scalar data'

   onc->AttAdd, /GLOBAL, 'history', $
      'Generated by routine Mgh_Roms_Movie_lslice::ExportAnimationDataToNcFile at '+ $
      mgh_dt_string(mgh_dt_now())

   onc->AttAdd, /GLOBAL, 'variable', variable

   onc->DimAdd, 'xi', dim[0]
   onc->DimAdd, 'eta', dim[1]
   onc->DimAdd, 'time'

   x_name = lonlat ? 'lon' : 'x'
   y_name = lonlat ? 'lat' : 'y'
   v_name = mgh_str_vanilla(variable)

   onc->VarAdd, x_name, ['xi','eta'], /DOUBLE
   onc->VarAdd, y_name, ['xi','eta'], /DOUBLE
   onc->VarAdd, v_name, ['xi','eta','time'], /FLOAT

   onc->VarPut, x_name, datax
   onc->VarPut, y_name, datay

   ;; Work through frames, retrieving and write data

   n_rec = 1 + (range[1]-range[0])/stride

   fmt ='(%"Writing %d frames of %d x %d data to netCDF file %s")'
   message, /INFORM, string(n_rec, dim, file, FORMAT=fmt)

   for r=0,n_rec-1 do begin

      pos = range[0]+stride*r

      oframe = self.animation->GetFrame(POSITION=pos)
      oframe[0]->GetProperty, KEYWORDS=keywords

      onc->VarPut, v_name, keywords.data_values, OFFSET=[0,0,r]

   endfor

   obj_destroy, onc

   fmt ='(%"Finished saving netCDF file %s")'
   message, /INFORM, string(file, FORMAT=fmt)

end

pro Mgh_Roms_Movie_Lslice__Define

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   struct_hide, $
      {mgh_roms_movie_lslice, inherits MGH_Datamator, $
       variable: '', map_structure: ptr_new(), $
       history_file: obj_new(), $
       bar: obj_new(), plane: obj_new(), contour: obj_new()}

end
