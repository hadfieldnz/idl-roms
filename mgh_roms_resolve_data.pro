;+
; ROUTINE NAME:
;   MGH_ROMS_RESOLVE_DATA
;
; PURPOSE:
;   For a specified ROMS variable, this function supplies default values for the
;   DATA_RANGE and DATA_MULTIPLIER keywords
;
; CALLING SEQUENCE:
;   mgh_roms_resolve_data, var, DATA_RANGE=data_range, DATA_MULTIPLIER=data_multiplier
;
; POSITIONAL PARAMETERS:
;   var (input, scalar string)
;     The name of a ROMS 2-D or 3-D variable.
;
; KEYWORD PARAMETERS:
;   DATA_MULTIPLIER (input/output, scalar numeric)
;     Number by which data values are multiplied before they displayed.
;
;   DATA_RANGE (input/output, 2-element numeric)
;     Data range for the graphics element in which data are displayed.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-12:
;     Written, based on & superseding function MGH_ROMS_DATA_RANGE.
;   Mark Hadfield, 2011-05:
;     Added support for several wave and sediment variables.
;   Mark Hadfield, 2016-03:
;     - Updated code
;     - Added 'sand*'.
;-
pro mgh_roms_resolve_data, var, $
     DATA_RANGE=data_range, DATA_MULTIPLIER=data_multiplier

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(var) ne 1 then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumelem', var

   if n_elements(data_multiplier) eq 0 then begin
      case !true of
         isa(var, /STRING) && strmatch(var, 'psi*'): $
            data_multiplier = 1.E-6
         else: $
            data_multiplier = 1
      endcase
   endif

   if n_elements(data_range) eq 0 then begin
      case !true of
         isa(var, /STRING) && var eq 'h': $
            data_range = [0,4000]
         isa(var, /STRING) && var eq 'bath': $
            data_range = [0,4000]
         isa(var, /STRING) && (var eq 'zeta' && var eq 'zeta_mean'): $
            data_range = [-1,1]
         isa(var, /STRING) && (var eq 'u' || var eq 'v' || var eq 'u_mean' || var eq 'v_mean'): $
            data_range = [-1,1]*0.3
         isa(var, /STRING) && (var eq 'ubar' || var eq 'vbar' || var eq 'ubar_mean' || var eq 'vbar_mean'): $
            data_range = [-1,1]*0.3
         isa(var, /STRING) && (var eq 'spd' || var eq 'spd_mean'): $
            data_range = [0,0.4]
         isa(var, /STRING) && (var eq 'sbar' || var eq 'sbar_mean'): $
            data_range = [0,0.4]
         isa(var, /STRING) && var eq 'w': $
            data_range = [-1,1]*1.E-3
         isa(var, /STRING) && var eq 'omega': $
            data_range = [-1,1]*1.E-3
         isa(var, /STRING) && (var eq 'temp' || var eq 'temp_mean'): $
            data_range = [0,25]
         isa(var, /STRING) && (var eq 'salt' || var eq 'salt_mean'): $
            data_range = [34,36]
         isa(var, /STRING) && var eq 'SST': $
            data_range = [0,25]
         isa(var, /STRING) && var eq 'SSS': $
            data_range = [34,36]
         isa(var, /STRING) && var eq 'AKs': $
            data_range = [0,1]
         isa(var, /STRING) && var eq 'AKt': $
            data_range = [0,1]
         isa(var, /STRING) && var eq 'AKv': $
            data_range = [0,1]
         isa(var, /STRING) && var eq 'gls': $
            data_range = [0,5E-3]
         isa(var, /STRING) && var eq 'rho': $
            data_range = [20,50]
         isa(var, /STRING) && (var eq 'psi' || var eq 'psi_mean'): $
            data_range = [-100,100]
         isa(var, /STRING) && var eq 'NO3': $
            data_range = [0,20]
         isa(var, /STRING) && var eq 'NH4': $
            data_range = [0,2]
         isa(var, /STRING) && var eq 'chlorophyll': $
            data_range = [0,5]
         isa(var, /STRING) && var eq 'phytoplankton': $
            data_range = [0,2]
         isa(var, /STRING) && var eq 'zooplankton': $
            data_range = [0,0.5]
         isa(var, /STRING) && var eq 'LdetritusN': $
            data_range = [0,2]
         isa(var, /STRING) && var eq 'SdetritusN': $
            data_range = [0,4]
         isa(var, /STRING) && var eq 'TIC': $
            data_range = [2000.,2200.]
         isa(var, /STRING) && var eq 'alkalinity': $
            data_range = [2200.,2400.]
         isa(var, /STRING) && var eq 'shflux': $
            data_range = [-200,200]
         isa(var, /STRING) && var eq 'Hsbl': $
            data_range = [-300,0]
         isa(var, /STRING) && var eq 'Hwave': $
            data_range = [0,2]
         isa(var, /STRING) && var eq 'Lwave': $
            data_range = [0,100]
         isa(var, /STRING) && var eq 'bed_thickness': $
            data_range = [0,1]
         isa(var, /STRING) && var eq 'bed_wave_amp': $
            data_range = [0,1]
         isa(var, /STRING) && strmatch(var, 'mud*'): $
            data_range = [0,10]
         isa(var, /STRING) && strmatch(var, 'sand*'): $
            data_range = [0,10]
         else: $
            data_range = [-1,1]
      endcase
   endif

end
