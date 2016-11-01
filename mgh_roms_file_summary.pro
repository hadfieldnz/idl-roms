;+
; NAME:
;   MGH_ROMS_FILE_SUMMARY
;
; PURPOSE:
;   This procedure prints a readable summary of the information in a ROMS netCDF
;   output file.
;
; CATEGORY:
;   Ocean models.
;
; CALLING SEQUENCE:
;   ROMS_FILE_SUMMARY, history
;
; POSITIONAL PARAMETERS:
;   file (input, scalar string or object reference)
;     A reference to a ROMS history-file or station-file object or a list of file names.
;
; KEYWORD PARAMETERS:
;   OUTPUT (input, optional, scalar string)
;     If a string value is supplied, print the output to a file with
;     this name, otherwise to the console.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 1999-01:
;       Written.
;   Mark Hadfield, 2004-10:
;       Added station-file support.
;-
pro mgh_roms_file_summary, file, OUTPUT=output

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Process history argument.

   case size(file, /TNAME) of
      'STRING': begin
         case mgh_roms_file_type(file) of
            'ROMS station file': $
                 ofile = obj_new('MGHromsStation', file)
            else: $
                 ofile = obj_new('MGHromsHistory', file)
         endcase
         ofile_destroy = 1B
      end
      'OBJREF': begin
         ofile = file
         ofile_destroy = 0B
      end
      else: $
           message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'file'
   endcase

   ;; Arrange a unit for output

   case size(output, /TNAME) of
      'STRING': begin
         message, /INFORM, 'Printing output to file '+output
         openw, lun, output, /GET_LUN
      end
      else: begin
         message, /INFORM, 'Printing output to standard output'
         lun = -1
      endelse
   endcase

   ;; Print info

   if ofile->HasAtt(/GLOBAL, 'type') then begin
      printf, lun, 'File type: '
      printf, lun, ofile->AttGet(/GLOBAL, 'type')
   endif

   if ofile->HasAtt(/GLOBAL, 'title') then begin
      printf, lun, 'Title: '
      printf, lun, ofile->AttGet(/GLOBAL, 'title')
   endif

   if ofile->HasAtt(/GLOBAL, 'CPP_options') then begin
      opt = ofile->AttGet(/GLOBAL, 'CPP_options')
      printf, lun, 'CPP options'
      printf, lun, strsplit(opt, ',', /EXTRACT)
      printf, lun, 'Wall parameters: '
      printf, lun, $
              max(strcmp(opt,'WESTERN_WALL' )) gt 0, $
              max(strcmp(opt,'SOUTHERN_WALL')) gt 0, $
              max(strcmp(opt,'EASTERN_WALL' )) gt 0, $
              max(strcmp(opt,'NORTHERN_WALL')) gt 0
   endif

   printf, lun, 'NetCDF dimension names: '
   printf, lun, ofile->DimNames()

   printf, lun, 'NetCDF variable names: '
   printf, lun, ofile->VarNames()

   if ofile->HasVar('dt') then begin
      printf, lun, 'Long time step: '
      printf, lun, ofile->VarGet('dt')
   endif

   if ofile->HasVar('dtfast') then begin
      printf, lun, 'Short time step: '
      printf, lun, ofile->VarGet('dtfast')
   endif

   if ofile->HasDim('time') then begin
      printf, lun, 'Number of time levels: '
      printf, lun, ofile->DimInfo('time', /DIMSIZE)
   endif

   if ofile->HasDim('station') then begin
      printf, lun, 'Number of stations: '
      printf, lun, ofile->DimInfo('station', /DIMSIZE)
   endif

   if ofile->HasVar('xl') && ofile->HasVar('el') then begin
      printf, lun, 'Domain size: '
      printf, lun, ofile->VarGet('xl'), ofile->VarGet('el')
   endif

   if ofile->HasDim('xi_rho') && ofile->HasDim('eta_rho') then begin
      printf, lun, 'Horizontal grid dimensions (interior rho points): '
      printf, lun, ofile->DimInfo('xi_rho', /DIMSIZE)-2, ofile->DimInfo('eta_rho', /DIMSIZE)-2
   endif

   blonlat = ofile->HasVar('lon_rho') && ofile->HasVar('lat_rho')

   x_var = blonlat ? 'lon_rho' : 'x_rho'
   y_var = blonlat ? 'lat_rho' : 'y_rho'

   if ofile->HasVar(x_var) && ofile->HasVar(y_var) then begin
      x_psi = mgh_stagger(ofile->VarGet(x_var), DELTA=[-1,-1])
      y_psi = mgh_stagger(ofile->VarGet(y_var), DELTA=[-1,-1])
      dim_psi = size(x_psi, /DIMENSIONS)
      printf, lun, 'Minimum & maximum x (PSI):'
      printf, lun, mgh_minmax(x_psi)
      printf, lun, 'Minimum & maximum y (PSI):'
      printf, lun, mgh_minmax(y_psi)
      printf, lun, 'SW corner (PSI):'
      printf, lun, x_psi[0,0], y_psi[0,0]
      printf, lun, 'SE corner (PSI):'
      printf, lun, x_psi[dim_psi[0]-1,0], y_psi[dim_psi[0]-1,0]
      printf, lun, 'NE corner (PSI):'
      printf, lun, x_psi[dim_psi[0]-1,dim_psi[1]-1], y_psi[dim_psi[0]-1,dim_psi[1]-1]
      printf, lun, 'NW corner (PSI):'
      printf, lun, x_psi[0,dim_psi[1]-1], y_psi[0,dim_psi[1]-1]
   endif

   if ofile->HasVar('angle') then begin
      angle = ofile->VarGet('angle')
      printf, lun, 'Maximum, minimum & average grid angle (deg):'
      printf, lun, [mgh_minmax(angle), mgh_avg(angle)]*!radeg
   endif

   if ofile->HasVar('f') then begin
      f = ofile->VarGet('f')
      printf, lun, 'Maximum, minimum & average Coriolis frequency (10^4 s-1):'
      printf, lun, [mgh_minmax(f), mgh_avg(f)]*1.D4
   endif

   if ofile->HasVar('pm') && ofile->HasVar('pn') then begin
      pm = ofile->VarGet('pm')
      pn = ofile->VarGet('pn')
      printf, lun, 'Maximum, minimum & average grid spacing (xi):'
      printf, lun, 1./[mgh_minmax(pm), mgh_avg(pm)]
      printf, lun, 'Maximum, minimum & average grid spacing (eta):'
      printf, lun, 1./[mgh_minmax(pn), mgh_avg(pn)]
      printf, lun, 'Maximum along-axis stretching ratio (xi & eta):'
      printf, lun, 10.^max(abs(mgh_diff(alog10(pm), 1))), $
                   10.^max(abs(mgh_diff(alog10(pn), 2)))
      ucrit = [8,10,12]
      for i=0,n_elements(ucrit)-1 do begin
         printf, lun, 'Baroclinic CFL limit (velocity = ',mgh_format_float(ucrit[i]),')'
         printf, lun, 1./(ucrit[i]*max(pm > pn))
      endfor
      if ofile->HasVar('h') then begin
         printf, lun, 'Barotropic CFL limit:'
         printf, lun, 1./max((pm > pn)*sqrt(9.8*ofile->VarGet('h')))
      endif
   endif

   if ofile->HasVar('h') then begin
      h = ofile->VarGet('h')
      if ofile->HasVar('mask_rho') then mask = ofile->VarGet('mask_rho')
      if size(h, /N_DIMENSIONS) eq 2 then begin
         printf, lun, 'Minimum and maximum depth:'
         printf, lun, mgh_minmax(h)
         printf, lun, 'Minimum and maximum r-value:'
         printf, lun, mgh_minmax(mgh_roms_r_value(h, MASK=mask), /NAN)
         printf, lun, 'Minimum and maximum slope:'
         printf, lun, mgh_minmax(abs(mgh_roms_slope(h, pm, pn)))
      endif
      mgh_undefine, h
   endif

   if ofile->HasVar('mask_rho') then begin
      mask = ofile->VarGet('mask_rho')
      printf, lun, 'Fraction of sea points: '
      printf, lun, total(mask)/n_elements(mask)
   endif

   if ofile->HasDim('s_rho') then begin
      printf, lun, 'Number of vertical grid levels: '
      printf, lun, ofile->DimInfo('s_rho', /DIMSIZE)
   endif

   if ofile->HasVar('theta_s') && ofile->HasVar('theta_b') && $
        ofile->HasVar('tcline') && ofile->HasVar('hc') then begin
      printf, lun, 'S coordinate parameters (theta_s, theta_b, tcline, hc): '
      printf, lun, $
              ofile->VarGet('theta_s'), ofile->VarGet('theta_b'), $
              ofile->VarGet('Tcline'), ofile->VarGet('hc')
   endif

   if ofile->HasDim('station') then begin
      ipos = ofile->VarGet('Ipos')
      jpos = ofile->VarGet('Jpos')
      xpos = ofile->VarGet(x_var)
      ypos = ofile->VarGet(y_var)
      hpos = ofile->VarGet('h')
      fmt = '(%"Station %d position %d, %d at %f, %f in a water depth of %f")'
      for p=0,n_elements(hpos)-1 do $
           printf, lun, FORMAT=fmt, p,ipos[p], jpos[p], xpos[p], ypos[p], hpos[p]
   endif

   if keyword_set(ofile_destroy) then obj_destroy, ofile

   if lun ne -1 then free_lun, lun

end

