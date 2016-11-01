; svn $Id$
function mgh_roms_slope, hs, pm, pn

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   dims = size(hs, /DIMENSIONS)

   ;; Calculate xi-slope on U grid and eta-slope on V grid

   su = - (hs[1:dims[0]-1,*]-hs[0:dims[0]-2,*])*mgh_stagger(pm, DELTA=[-1,0])

   sv = - (hs[*,1:dims[1]-1]-hs[*,0:dims[1]-2])*mgh_stagger(pn, DELTA=[0,-1])

   ;; Interpolate back onto RHO grid and return components, bundled into
   ;; a complex number.

   return, complex(mgh_stagger(temporary(su), DELTA=[1,0]), $
                   mgh_stagger(temporary(sv), DELTA=[0,1]))

end
