;+
; NAME:
;   MGH_ROMS_BATHSUDS
;
; PURPOSE:
;   Smooth bathymetry on a ROMS grid
;
; CALLING SEQUENCE:
;   Result = MGH_ROMS_BATHSUDS(hs, ar)
;
; POSITIONAL PARAMETERS:
;   hs (input, 2D numeric array)
;     Depths at rho points.
;
; RETURN VALUE:
;   The function returns a smoothed array with the same shape as the original.
;
;   I think it uses a Shapiro filter.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2000-06:
;     Written.
;   Mark Hadfield, 2003-07:
;     Changed provision for controlling weighting of smoother (in a
;     backwards-incompatible way): the r-value above which weighting
;     is applied is now specified directly via the THRESHOLD keyword.
;   Mark Hadfield, 2009-10:
;     Removed call to widget_event(/NOWAIT).
;-

function mgh_roms_bathsuds_pad, h

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE
   compile_opt HIDDEN

   ;; Given a 2D array, return a copy with an extra row/column around
   ;; the edge, padded with adjacent interior values.

   dims = size(h, /DIMENSIONS)

   ;; Create a copy of the input data,

   result = fltarr(dims[0]+2, dims[1]+2)
   result[1:dims[0],1:dims[1]] = h

   result[0,1:dims[1]] = result[1,1:dims[1]]
   result[dims[0]+1,1:dims[1]] = result[dims[0],1:dims[1]]
   result[*,0] = result[*,1]
   result[*,dims[1]+1] = result[*,dims[1]]

   return, result

end


function mgh_roms_bathsuds, h, $
     KEEP_SHALLOW=keep_shallow, N_PASSES=n_passes, THRESHOLD=threshold

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   ;; Check inputs

   if size(h, /N_DIMENSIONS) ne 2 then $
        message, 'The depth array must be two-dimensional.'

   if n_elements(n_passes) eq 0 then n_passes = 10

   if n_elements(threshold) eq 0 then threshold = 0.15

   ;; Dimensions etc

   dims = size(h, /DIMENSIONS)

   n_shallow = n_elements(keep_shallow)

   ;; Define kernel

   kernel = fltarr(3,3)

   a = 0.3

   kernel[0,0] = a^2
   kernel[1,0] = a
   kernel[2,0] = a^2
   kernel[0,1] = a
   kernel[1,1] = 1
   kernel[2,1] = a
   kernel[0,2] = a^2
   kernel[1,2] = a
   kernel[2,2] = a^2

   kernel = kernel/total(kernel)

   ;; Save locations of shallow areas

   if n_shallow gt 0 then begin

      mask_shallow = bytarr([dims,n_shallow])

      for k=0,n_shallow-1 do $
            mask_shallow[*,*,k] = h lt keep_shallow[k]

   endif

   ;; Copy input to output

   hnew = h

   ;; Apply the filter n_passes times

   for pass=0,n_passes-1 do begin

      ;; Keep shallow areas shallow

      if (pass gt 0) and (pass le n_passes-3) then begin

         for k=0,n_shallow-1 do begin
            hsh = keep_shallow[k]
            ww = where(mask_shallow[*,*,k] and (hnew gt hsh), count)
            if count gt 0 then hnew[ww] = hsh
         endfor

      endif

      ;; The delta array will hold the change calculated
      ;; by applying the filter.

      delta = replicate(!values.f_nan, dims[0:1])

      ;; Pad the depth array

      hpad = mgh_roms_bathsuds_pad(hnew)

      ;; Work through points one at a time
      ;; Indices i & j point to locations in the original data

      for i=0,dims[0]-1 do for j=0,dims[1]-1 do begin

         ;; Extract a 3x3 local subset from the padded 2D arrays

         d = hpad[i:i+2,j:j+2]-hpad[i+1,j+1]

         ;; At each point calculate delta

         delta[i,j] = total(d * kernel)

      endfor

      mgh_undefine, hpad

      ;; Weight delta

      rvw = 0.5 * (1 + tanh((mgh_roms_r_value(hnew)-threshold)/0.03))

      hnew = hnew + rvw * delta

   endfor

   return, hnew

end
