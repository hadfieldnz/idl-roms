;+
; NAME:
;   MGH_ROMS_R_VALUE
;
; PURPOSE:
;   Calculate grid stiffness in the sense of Beckmann and Haidvogel (1993)
;
; POSITIONAL PARAMETERS:
;   hs (input, numeric array, 1-D or 2-D)
;     Depth at rho points.
;
; RETURN VALUE:
;   The function returns a floating/double array with the same size and
;   shape as the input, containing the stiffness of the grid.
;
; PROCEDURE:
;   Calculate values on a staggered grid, then project back onto original grid
;   by taking maximum neighbouring value,
;
;###########################################################################
; Copyright (c) 2015 NIWA:
;   http://www.niwa.co.nz/
; Licensed under the MIT open source license:
;   http://www.opensource.org/licenses/mit-license.php
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2000-??:
;     Written.
;   Mark Hadfield, 2003-06:
;     Add code to handle 1-D arrays, and simplify code.
;   Mark Hadfield, 2011-06:
;     Added MASK keyword.
;-
function mgh_roms_r_value, hs, MASK=mask

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   n_dim = size(hs, /N_DIMENSIONS)

   dim = size(hs, /DIMENSIONS)

   case n_dim of

      1: begin

         ;; Calculate r-value at intermediate points and pad result with 0s at each end

         rr = [0.,abs(1.0*(hs[1:dim[0]-1]-hs[0:dim[0]-2])/(hs[1:dim[0]-1]+hs[0:dim[0]-2])),0.]

         ;; Return, maximum neighbouring each point.

         result = rr[0:dim[0]-1] > rr[1:dim[0]]

      end

      2: begin

         ;; Could probably simplify this on the same lines as 1-D code.

         ;; Create array to hold 4 sets of one-sided r values

         r = fltarr(dim[0], dim[1], 4)

         ;; Create a copy of the input data, with an extra row/column around
         ;; the edge, padded with adjacent interior values.

         hsm = fltarr(dim[0]+2, dim[1]+2)
         hsm[1:dim[0],1:dim[1]] = hs
         hsm[0,1:dim[1]] = hsm[1,1:dim[1]]
         hsm[dim[0]+1,1:dim[1]] = hsm[dim[0],1:dim[1]]
         hsm[*,0] = hsm[*,1]
         hsm[*,dim[1]+1] = hsm[*,dim[1]]

         for i=0,3 do begin

            case i of
               0: tmp = hsm[1:dim[0],0:dim[1]-1] ;;; Shifted -1 in x direction
               1: tmp = hsm[1:dim[0],2:dim[1]+1] ;;; Shifted +1 in x direction
               2: tmp = hsm[0:dim[0]-1,1:dim[1]] ;;; Shifted -1 in y direction
               3: tmp = hsm[2:dim[0]+1,1:dim[1]] ;;; Shifted +1 in y direction
            endcase

            r[*,*,i] = abs((hs-tmp)/(hs+tmp))

         endfor

         result = r[*,*,0] > r[*,*,1] > r[*,*,2] > r[*,*,3]

      end

   endcase

   if n_elements(mask) gt 0 then begin

      if ~array_equal(size(mask, /DIMENSIONS), dim) then $
          message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgdimsize', 'mask'

      l_mask = where(mask eq 0, n_mask)
      if n_mask gt 0 then result[l_mask] = !values.f_nan

   endif

   return, result

end
