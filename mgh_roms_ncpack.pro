;+
; NAME:
;   MGH_ROMS_NCPACK
;
; PURPOSE:
;   Copy a ROMS netCDF file, applying packing to selected variables.
;   Any packed variables in the input file are unpacked and repacked
;   with no problem.
;
;   This routine has been superseded by Python script rncpack5, which
;   is more suited for embedding in a shell script, though it does not
;   offer fine control over the packing data ranges.
;
; CATEGORY:
;   ROMS
;
; CALLING SEQUENCE:
;   MGH_ROMS_NCPACK, file_in, file_out
;
; POSITIONAL PARAMETERS:
;   file_in (input, string scalar)
;     The name of a ROMS netCDF file to be read.
;
;   file_out (input, string scalar, optional)
;     The name of a ROMS netCDF file to be written with packed
;     data. The default is the input file name with ".pack" inserted
;     before the extension.
;
; KEYWORD PARAMETERS:
;   BSTR_RANGE, HBBL_RANGE, HSBL_RANGE, OMEGA_RANGE, SALT_RANGE,
;   TEMP_RANGE, UV_RANGE, UVBAR_RANGE, W_RANGE, ZETA_RANGE (input,
;   two-element numeric vector)
;     Range of allowed values for the specified variables or sets of
;     variables.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2010-07:
;     Written.
;   Mark Hadfield, 2010-11:
;     - Added packing parameters for bottom stress.
;     - Added support for detided and tide-harmonic variables.
;     - Streamlined variable matching with a dedicated function.
;   Mark Hadfield, 2010-11:
;     - Disabled packing of vertical diffusivity (AKt and AKv)
;       because of the large dynamic range.
;   Mark Hadfield, 2011-01:
;     - Default output file name now ends in "_pack.nc"
;   Mark Hadfield, 2011-02:
;     - Enabled packing of omega.
;   Mark Hadfield, 2011-03:
;     - Disabled packing of omega and w: too many out-of-bound values.
;-
function mgh_roms_ncpack_strmatch, varname, pattern

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE
   compile_opt HIDDEN

   for p=0,n_elements(pattern)-1 do $
         if strmatch(varname, pattern[p]) then return, 1B

   return, 0B

end


pro mgh_roms_ncpack, file_in, file_out, $
     AK_RANGE=ak_range, BSTR_RANGE=bstr_range, HBBL_RANGE=hbbl_range, $
     HSBL_RANGE=hsbl_range, OMEGA_RANGE=omega_range, SALT_RANGE=salt_range, $
     TEMP_RANGE=temp_range, UV_RANGE=uv_range, UVBAR_RANgE=uvbar_range, $
     W_RANGE=w_range, ZETA_RANGE=zeta_range

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(file_in) eq 0 then $
        message, 'Name for input file not supplied'

   if ~ file_test(file_in, /READ) then $
        message, 'Input file cannot be read'

   if n_elements(file_out) eq 0 then begin
      if ~ strmatch(file_in, '*.nc') then $
           message, 'Cannot generate default output name if input name ' + $
                    'does not end with ".nc"'
      file_out = strmid(file_in, 0, strlen(file_in)-3)+'_pack.nc'
   endif

;   if n_elements(ak_range) eq 0 then ak_range = [0,10]
   if n_elements(bstr_range) eq 0 then bstr_range = [-100,100]
   if n_elements(hbbl_range) eq 0 then hbbl_range = [-6000,10]
   if n_elements(hsbl_range) eq 0 then hsbl_range = [-6000,10]
   if n_elements(omega_range) eq 0 then omega_range = [-0.05,0.05]
   if n_elements(salt_range) eq 0 then salt_range = [0,40]
   if n_elements(temp_range) eq 0 then temp_range = [-5,45]
   if n_elements(uv_range) eq 0 then uv_range = [-10,10]
   if n_elements(uvbar_range) eq 0 then uvbar_range = [-10,10]
   if n_elements(w_range) eq 0 then w_range = [-0.2,0.2]
   if n_elements(zeta_range) eq 0 then zeta_range = [-10,10]

   ;; Open the files

   message, /INFORM, string(FORMAT='(%"Opening input file %s")', file_in)
   onc0 = obj_new('MGHncReadFile', file_in)

   message, /INFORM, string(FORMAT='(%"Creating output file %s")', file_out)
   onc1 = obj_new('MGHncFile', file_out, /CREATE, /CLOBBER)

   ;; Copy dimensions and global attributes & update history

   message, /INFORM, 'Copying dimensions & global attributes'
   onc1->DimCopy, onc0
   onc1->AttCopy, onc0, /GLOBAL

   history = onc1->AttGet('history', /GLOBAL)
   fmt = '(%"%s: Packed with MGH_ROMS_NCPACK%s%s")'
   history = string(FORMAT=fmt, $
                    mgh_dt_string(mgh_dt_now()), string(10B), history)
   onc1->AttAdd, /GLOBAL, 'history', temporary(history)

   ;; Work through variables copying definitions and attributes (suitably
   ;; modified for packed variables)

   message, /INFORM, 'Copying variable definitions and attributes'

   var = onc0->VarNames(COUNT=n_var)

   pack_short = bytarr(n_var)

   for v=0,n_var-1 do begin

      onc0->VarInfo, var[v], ALL=info

      case 1B of
;         mgh_roms_ncpack_strmatch(var[v], ['AKt','AKv']): begin
;            pack_short[v] = 1B
;            data_range = ak_range
;         end
         mgh_roms_ncpack_strmatch(var[v], ['bustr','bvstr']): begin
            pack_short[v] = 1B
            data_range = bstr_range
         end
         mgh_roms_ncpack_strmatch(var[v], ['Hbbl']): begin
            pack_short[v] = 1B
            data_range = hbbl_range
         end
         mgh_roms_ncpack_strmatch(var[v], ['Hsbl']): begin
            pack_short[v] = 1B
            data_range = hsbl_range
         end
;         mgh_roms_ncpack_strmatch(var[v], ['omega']): begin
;            pack_short[v] = 1B
;            data_range = omega_range
;         end
         mgh_roms_ncpack_strmatch(var[v], ['salt','salt_*tide*']): begin
            pack_short[v] = 1B
            data_range = salt_range
         end
         mgh_roms_ncpack_strmatch(var[v], ['temp','temp_*tide*']): begin
            pack_short[v] = 1B
            data_range = temp_range
         end
         mgh_roms_ncpack_strmatch(var[v], ['u','v','u_*tide*','v_*tide*']): begin
            pack_short[v] = 1B
            data_range = uv_range
         end
         mgh_roms_ncpack_strmatch(var[v], $
              ['ubar','vbar','ubar_*tide*','vbar_*tide*']): begin
            pack_short[v] = 1B
            data_range = uvbar_range
         end
;         mgh_roms_ncpack_strmatch(var[v], ['w']): begin
;            pack_short[v] = 1B
;            data_range = w_range
;         end
         mgh_roms_ncpack_strmatch(var[v], ['zeta','zeta_*tide*']): begin
            pack_short[v] = 1B
            data_range = zeta_range
         end
         else:
      endcase

      ;; Create variable

      nctype = pack_short[v] ? 'SHORT' : info.datatype

      onc1->VarAdd, var[v], info.dim_names, NCTYPE=nctype

      ;; Copy attributes and add ones relevant to packing. The AttAdd
      ;; method overwrites any existing attributes with the same name.

      onc1->VarCopy, onc0, var[v], /ATTRIBUTES

      if keyword_set(pack_short[v]) then begin
         scale = mgh_norm_coord([-32766S,32767S], data_range)
         onc1->AttAdd, var[v], 'add_offset', scale[0]
         onc1->AttAdd, var[v], 'scale_factor', scale[1]
         onc1->AttAdd, var[v], 'valid_range', [-32766S,32767S]
         onc1->AttAdd, var[v], '_FillValue', -32767S
      endif

   endfor

   mgh_undefine, scale

   ;; Copy variable data.

   message, /INFORM, 'Copying variable data'

   for v=0,n_var-1 do begin

      onc0->VarInfo, var[v], ALL=info

      if keyword_set(pack_short[v]) then begin
         scale = [onc1->AttGet(var[v], 'add_offset'), $
                  onc1->AttGet(var[v], 'scale_factor')]
         fill = onc1->AttGet(var[v], '_FillValue')
      endif

      ;; If the variable has 2 or more dimensions, the last of which
      ;; is unlimited dimension, we will copy one record at a time.

      use_rec = 0B
      if info.n_dims gt 1 then $
           use_rec = onc0->DimInfo(info.dim_names[info.n_dims-1], /IS_UNLIMITED)

      n_rec = use_rec ? info.dimensions[info.n_dims-1] : 1

      for r=0,n_rec-1 do begin

         case 1B of
            use_rec: begin
               count = [info.dimensions[0:info.n_dims-2],1]
               offset = [replicate(0, info.n_dims-1),r]
            end
            info.n_dims gt 0: begin
               count = info.dimensions
               offset = replicate(0, info.n_dims)
            end
            else: begin
               count = 1
               offset = 0
            endelse
         endcase

         if keyword_set(pack_short[v]) then begin
            data = onc0->VarGet(var[v], AUTOSCALE=1, $
                                OFFSET=offset, COUNT=count)
            idata = mgh_reproduce(fill, data)
            l_finite = where(finite(data), n_finite)
            if n_finite gt 0 then begin
               dd = data[l_finite]
               ii = round((dd-scale[0])/scale[1])
               if max(ii) gt 32767 || min(ii) lt -32766 then begin
                  fmt = '(%"Clipping out-of-bounds values for variable %s ' + $
                        '(range(%f to %f)")'
                  message, /INFORM, $
                           string(FORMAT=fmt, var[v], mgh_minmax(dd))
                  ii = (ii > (-32766)) < 32767
               endif
               idata[l_finite] = temporary(ii)
            endif
            onc1->VarPut, var[v], temporary(idata), $
                          OFFSET=offset, COUNT=count
            mgh_undefine, data
         endif else begin
            data = onc0->VarGet(var[v], AUTOSCALE=0, $
                                OFFSET=offset, COUNT=count)
            onc1->VarPut, var[v], temporary(data), $
                          OFFSET=offset, COUNT=count
         endelse

      endfor

   endfor

   ;; Clean up

   obj_destroy, [onc0,onc1]

end

