;+
; FUNCTION NAME:
;   MGH_ROMS_SERIES_INTEGRAL
;
; PURPOSE:
;   This function calculates and returns a time series of the volume, area
;   or bed integral of a ROMS variable, evaluated over the domain or a rectangular
;   subset thereof.
;
;   The variable must be appropriately dimensioned:
;     - For a variable dimensioned ['xi_rho','eta_rho','s_rho'] and optionally time
;       a volume integral is calculated
;     - For a variable dimensioned ['xi_rho','eta_rho','N_bed'] and optionally time
;       a bed integral is calculated
;     - For a variable dimensioned ['xi_rho','eta_rho'] and optionally time
;       an area integral is calculated
;
;   The integral at each time step is calculated with the MGH_ROMS_INTEGRAL_VOLUME,
;   MGH_ROMS_INTEGRAL_AREA or MGH_ROMS_INTEGRAL_BED function, as appropriate.
;
;   Note that I normally call this function via higher-level wrappers like
;   MGH_COOK_SERIES_INTEGRAL, which adds temporary file support.
;
; CALLING SEQUENCE
;   result = mgh_roms_series_integral(ohis)
;
; RETURN VALUE:
;   The function returns the time series and ancillary information in
;   a structure.
;
; POSITIONAL PARAMETERS:
;   ohis (input, object)
;     A reference to an MGHromsHistory object.
;
; KEYWORD PARAMETERS:
;   RECORD_RANGE (integer 2-element vector)
;   TIME_RANGE (numeric 2-element vector)
;     Specify the range of records to be extracted either as either indices
;     (RECORD_RANGE) or times (in days) from the simulation's reference time
;     (TIME_RANGE).
;
;   USE_ZETA (input, switch)
;     Use zeta data? Default is to do so if available.
;
;   VARIABLE (input, string or structure scalar)
;     A variable descriptor, must be a string or a structure that can
;     be interpreted by the MGHromsHistory object's  methods. Default is 'temp'.
;
;   XI_RANGE (integer 2-element vector)
;   ETA_RANGE (integer 2-element vector)
;     Used to defined subsets of the domain in xi or eta directions.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2011-08:
;     Written.
;   Mark Hadfield, 2014-06:
;     Documentation improvements.
;-
function mgh_roms_series_integral, ohis, position, $
     RECORD_RANGE=record_range, TIME_RANGE=time_range, $
     XI_RANGE=xi_range, ETA_RANGE=eta_range, $
     USE_ZETA=use_zeta, VARIABLE=variable

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  ;; Process file argument

  if n_elements(ohis) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'ohis'

  if n_elements(ohis) gt 1 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', 'ohis'

  if ~ obj_valid(ohis) gt 1 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_objref_bad', 'ohis'

  ;; Process VARIABLE argument

  if n_elements(variable) eq 0 then variable = 'temp'

  ;; Get & check variable dimensions.

  dim = ohis->VarDims(variable)

  if ~ array_equal(dim.horizontal, ['xi_rho','eta_rho']) then begin
    fmt = '(%"The variable %s must be dimensioned by xi_rho & eta_rho")'
    message, string(FORMAT=fmt, string(variable, /PRINT))
  endif

  case dim.vertical of
    's_rho': begin
      has_vertical = 1B
    end
    '': begin
      has_vertical = 0B
    end
    else: begin
      fmt = '(%"The vertical dimension of variable %s, if any, must be s_rho")'
      message, string(FORMAT=fmt, string(variable, /PRINT))
    end
  endcase

  case dim.bed of
    'Nbed': begin
      has_bed = 1B
    end
    '': begin
      has_bed = 0B
    end
    else: begin
      fmt = '(%"The bed dimension of variable %s, if any, must be Nbed")'
      message, string(FORMAT=fmt, string(variable, /PRINT))
    end
  endcase

  has_time = strlen(dim.time) gt 0

  fmt = '(%"Calculating integrals of %s data")'
  message, /INFORM, string(FORMAT=fmt, string(variable, /PRINT))

  ;; Use zeta data?

  if ~ has_bed then begin
    if n_elements(use_zeta) eq 0 then $
      use_zeta = has_time ? ohis->HasVar('zeta') : 0B
    if use_zeta && (~ has_time) then $
      message, 'The USE_ZETA keyword is not currently supported for time-invariant variables'
  endif

  ;; Some final checks on keyword consistency

  if has_bed and has_vertical then $
    message, "I wasn't expecting that!"

  if has_bed and keyword_set(use_zeta) then $
    message, "I wasn't expecting that!"

  ;; Read grid dimensions & calculate horizontal ranges. Default is all
  ;; interior points. Protect they keyword arguments from modification

  dim_rho = ohis->DimRho()

  xr = n_elements(xi_range) gt 0 ? xi_range : [1,dim_rho[0]-2]
  er = n_elements(eta_range) gt 0 ? eta_range : [1,dim_rho[1]-2]

  xr0 = xr[0]
  xr1 = xr[1]
  xrn = xr1 - xr0 + 1

  er0 = er[0]
  er1 = er[1]
  ern = er1 - er0 + 1

  mgh_undefine, xr, er

  ;; Get horizontal grid data

  if ohis->HasVar('mask_rho') then $
    mask = ohis->VarGet('mask_rho', OFFSET=[xr0,er0], COUNT=[xrn,ern])

  h = ohis->VarGet('h', OFFSET=[xr0,er0], COUNT=[xrn,ern])

  pm = ohis->VarGet('pm', OFFSET=[xr0,er0], COUNT=[xrn,ern])
  pn = ohis->VarGet('pn', OFFSET=[xr0,er0], COUNT=[xrn,ern])

  ;; Get vertical coordinate data

  if has_vertical then begin
    scoord = {theta_s: ohis->VarGet('theta_s'), $
              theta_b: ohis->VarGet('theta_b'), $
              hc: ohis->VarGet('hc'), $
              vstretch: ohis->HasVar('Vstretching') ? ohis->VarGet('Vstretching'): 1, $
              vtransform: ohis->HasVar('Vtransform') ? ohis->VarGet('Vtransform'): 1}
  endif

  if has_time then begin

    ;; Process variable with time dimension

    n_time = ohis->DimInfo(dim.time, /DIMSIZE)
    mgh_resolve_indices, n_time, record_range

    case 1B of
      ohis->HasVar(dim.time): $
        time_var = dim.time
      ohis->HasVar('ocean_time'): $
        time_var = 'ocean_time'
      ohis->HasVar('scrum_time'): $
        time_var = 'scrum_time'
      else: $
        message, 'Time variable not found'
    endcase

    time = ohis->VarGet(time_var, /AUTOSCALE)

    if n_elements(time_range) gt 0 then begin
      record_range = mgh_subset(time, time_range)
    endif
    if n_elements(record_range) eq 0 then begin
      n_time = ohis->DimInfo(dim_u.time, /DIMSIZE)
      record_range = [0,n_time-1]
    endif

    rr0 = record_range[0]
    rr1 = record_range[1]
    rrn = rr1 - rr0 + 1

    time = time[rr0:rr1]

    integral = dblarr(rrn)

    case 1B of

      has_vertical: begin
        if use_zeta then volume = dblarr(rrn)
        for r=rr0,rr1 do begin
          data = ohis->VarGet(variable, OFFSET=[xr0,er0,0,r], COUNT=[xrn,ern,0,1])
          if use_zeta then $
            zeta = ohis->VarGet('zeta', OFFSET=[xr0,er0,r], COUNT=[xrn,ern,1])
          integral[r-rr0] = mgh_roms_integral_volume(data, pm, pn, h, scoord, $
                                                     MASK=mask, ZETA=zeta, VOLUME=vol)
          if use_zeta then begin
            volume[r-rr0] = vol
          endif else begin
            if r eq rr0 then volume = vol
          endelse
        endfor
        result = {time: temporary(time), integral: temporary(integral), volume: temporary(volume)}
      end

      has_bed: begin
        for r=rr0,rr1 do begin
          data = ohis->VarGet(variable, OFFSET=[xr0,er0,0,r], COUNT=[xrn,ern,0,1])
          integral[r-rr0] = mgh_roms_integral_bed(data, pm, pn, MASK=mask)
        endfor
        result = {time: temporary(time), integral: temporary(integral)}
      end

      else: begin
        ;; No vertical or bed dimension
        for r=rr0,rr1 do begin
          data = ohis->VarGet(variable, OFFSET=[xr0,er0,r], COUNT=[xrn,ern,1])
          if r eq rr0 then begin
            integral[r-rr0] = mgh_roms_integral_area(data, pm, pn, MASK=mask, AREA=area)
          endif else begin
            integral[r-rr0] = mgh_roms_integral_area(data, pm, pn, MASK=mask)
          endelse
        endfor
        result = {time: temporary(time), integral: temporary(integral), area: temporary(area)}
      end

    endcase

  endif else begin

    ;; Process variable with no time dimension

    case 1B of

      has_vertical: begin
        data = ohis->VarGet(var, OFFSET=[xr0,er0,0], COUNT=[xrn,ern,0])
        integral = mgh_roms_integral_volume(data, pm, pn, h, scoord, $
          MASK=mask, ZETA=zeta, VOLUME=volume)
        result = {integral: temporary(integral), volume: temporary(volume)}
      end

      has_bed: begin
        data = ohis->VarGet(var, OFFSET=[xr0,er0,0], COUNT=[xrn,ern,0])
        integral = mgh_roms_integral_bed(data, pm, pn, MASK=mask)
        result = {integral: temporary(integral)}
      end

      else: begin
        ;; No vertical or bed dimension
        data = ohis->VarGet(var, OFFSET=[xr0,er0], COUNT=[xrn,ern])
        integral = mgh_roms_integral_area(data, pm, pn, MASK=mask, AREA=area)
        result = {integral: temporary(integral), area: temporary(area)}
      end
    endcase

  endelse

  return, result

end
