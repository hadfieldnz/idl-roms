;+
; NAME:
;   MGH_ROMS_CURL
;
; PURPOSE:
;   Calculate curl of a ROMS vector field
;
; POSITIONAL PARAMETERS:
;   u (input, numeric 2-D array)
;     Vector component in xi direction, defined at u points
;
;   v (input, numeric 2-D array)
;     Vector component in eta direction, defined at v points
;
;   pm (input, numeric 2-D array)
;     Grid metric in xi direction, defined at rho points.
;
;   pn (input, numeric 2-D array)
;     Grid metric in eta direction, defined at rho points.
;
; KEYWORD PARAMETERS:
;   mask (input, numeric 2-D array)
;     Sea mask (0=land, 1=sea) defined at rho points.
;
; RETURN VALUE:
;   The function returns a double-precision 2D array containing curl values,
;   defined at psi points.
;
; PROCEDURE:
;   Straightforward implementation of Wilkin, Mainbridge & Hedstrom
;   (1995) equation 10.
;
;###########################################################################
; Copyright (c) 2003-2015 NIWA:
;   http://www.niwa.co.nz/
; Licensed under the MIT open source license:
;   http://www.opensource.org/licenses/mit-license.php
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2003-07:
;     Written.
;   Mark Hadfield, 2003-08:
;     Improved formulation: the original one did not account correctly
;     for variations in pm & pn.
;-
function mgh_roms_curl, u, v, pm, pn

  compile_opt DEFINT32
  compile_opt STRICTARR
  compile_opt STRICTARRSUBS
  compile_opt LOGICAL_PREDICATE

  if size(u, /N_ELEMENTS) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'u'

  if size(v, /N_ELEMENTS) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'v'

  if size(pm, /N_ELEMENTS) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'pm'

  if size(pn, /N_ELEMENTS) eq 0 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'pn'

  if size(u, /N_DIMENSIONS) ne 2 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'u'

  if size(v, /N_DIMENSIONS) ne 2 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'v'

  if size(pm, /N_DIMENSIONS) ne 2 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pm'

  if size(pn, /N_DIMENSIONS) ne 2 then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'pn'

  dim = size(u, /DIMENSIONS) + [1,0]

  d0 = dim[0]  &  d1 = dim[1]

  if ~ array_equal(size(v, /DIMENSIONS), dim-[0,1]) then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'u', 'v'

  if ~ array_equal(size(pm, /DIMENSIONS), dim) then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'u & v', 'pm'

  if ~ array_equal(size(pn, /DIMENSIONS), dim) then $
    message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'u & v', 'pn'

  ;; Here we go!

  pm_psi = mgh_stagger(pm, DELTA=[-1,-1])
  pn_psi = mgh_stagger(pn, DELTA=[-1,-1])

  pm_u = mgh_stagger(pm, DELTA=[-1,0])

  pn_v = mgh_stagger(pn, DELTA=[0,-1])

  result = ((v[1:d0-1,*]/pn_v[1:d0-1,*]-v[0:d0-2,*]/pn_v[0:d0-2,*]) - $
            (u[*,1:d1-1]/pm_u[*,1:d1-1]-u[*,0:d1-2]/pm_u[*,0:d1-2])) $
           * pm_psi*pn_psi

  ;; Set values along exterior boundary to zero, as they involve velocities
  ;; "outside" the domain. Alternative is to use 1-sided derivatives.

  result[0,*] = 0
  result[d0-2,*] = 0
  result[*,0] = 0
  result[*,d1-2] = 0

  return, result

end
