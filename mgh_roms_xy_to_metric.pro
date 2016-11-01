;+
; NAME:
;   MGH_ROMS_XY_TO_METRIC
;
; PURPOSE:
;   Given x & y, or lon & lat, data, calculate coordinate metric data.
;
; CALLING SEQUENCE:
;   MGH_ROMS_XY_TO_METRIC, x, y,
;        ANGLE=angle, LONLAT=lonlat, PM=pm, PN=pn, DNDX=dndx, DMDE=dmde
;
; POSITIONAL PARAMETERS:
;   x & y (input, numeric arrays)
;     X & Y positions. The size and shape required depend on the
;     setting of the GRID keyword.
;
; KEYWORD PARAMETERS:
;   ANGLE (output, 2D array)
;     Angle of local xi direction relative to x/east axis (radians anticlockwise)
;
;   GRID (input, switch)
;     Set this keyword if x & y arrays are 1D arrays describing a rectangular
;     grid. Otherwise x & y must be 2D arrays which match in size and shape.
;     It appears that calculations for the GRID=0 are incomplete.
;
;   LONLAT (input, integer)
;     Controls whether x & y are to be interpreted as Cartesian
;     coordinates in m (LONLAT=0), longitude and latitude in degrees
;     on a spherical earth (LONLAT=1) or longitude and latitude in
;     degrees on an elliptical earth (LONLAT=2)
;
;   PM (output, 2D array)
;     Curvilinear coordinate metric in xi (m^-1)
;
;   DNDX (output, 2D array)
;     Xi derivative of inverse metric factor 1/pn (m)
;
;   DMDE (output, 2D array matching x & y)
;     Eta derivative of inverse metric factor 1/pm (m)
;
; EXPLANATION:
;   The curvilinear coordinate metrics are described in the Gridpak
;   manual at:
;
;     ftp://ahab.rutgers.edu/pub/gridpak/grid_manual.ps.gz
;
;   See sphere.F and geodesic_dist.F in the gridpak package for Fortran
;   code to the same thing.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2002-12:
;     Written.
;   Mark Hadfield, 2003-09:
;     Added angle code.
;-
pro mgh_roms_xy_to_metric, x, y, $
     ANGLE=angle, DNDX=dndx, DMDE=dmde, GRID=grid, LONLAT=lonlat, PM=pm, PN=pn

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(x) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'x'

   if n_elements(y) eq 0 then $
        message, BLOCK='mgh_mblk_motley', NAME='mgh_m_undefvar', 'y'

   if n_elements(lonlat) eq 0 then lonlat = 1

   if keyword_set(grid) then begin

      if size(x, /N_DIMENSIONS) ne 1 then $
         message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'x'
      if size(y, /N_DIMENSIONS) ne 1 then $
         message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'y'

      d0 = size(x, /N_ELEMENTS)
      d1 = size(y, /N_ELEMENTS)

      dim = [d0,d1]

   endif else begin

      if size(x, /N_DIMENSIONS) ne 2 then $
         message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'x'
      if size(y, /N_DIMENSIONS) ne 2 then $
         message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrgnumdim', 'y'

      dim = size(x, /DIMENSIONS)

      d0 = dim[0]
      d1 = dim[1]

      if ~ array_equal(dim, size(x, /DIMENSIONS)) then $
         message, BLOCK='mgh_mblk_motley', NAME='mgh_m_mismatcharr', 'x', 'y'

   endelse

   ;; We want to handle precision issues in a flexible way. Establish
   ;; if any of the inputs are in double precision and make appropriate
   ;; definitions

   ftype = size(0.E, /TYPE)
   dtype = size(0.D, /TYPE)

   double = size(x, /TYPE) eq dtype || size(y, /TYPE) eq dtype

   if keyword_set(double) then begin
      one = 1.0D
      im = dcomplex(0, 1)
      pi = !dpi
      rtype = dtype
      radius = 6371.D3
   endif else begin
      one = 1.0
      im = complex(0, 1)
      pi = !pi
      rtype = ftype
      radius = 6371.E3
   endelse

   ;; Calculate pm, pn & angle, using distance functions that depend on the
   ;; lonlat setting. There are three lonlat options:
   ;;   0: Cartesian grid;
   ;;   1: Lon/lat grid, spherical earth;
   ;;   2: Lon/lat grid, ellipsoidal earth.

   case lonlat of

      0: begin

         message, 'Calculating grid spacing & angle from x & y data', /INFORM

         if keyword_set(grid) then begin
            pm = mgh_inflate(dim, mgh_stagger(one/mgh_diff(x), DELTA=1), 1)
            pn = mgh_inflate(dim, mgh_stagger(one/mgh_diff(y), DELTA=1), 2)
            angle = make_array(dim, TYPE=rtype)
            message, 'This option not completely implemented!'
         endif else begin
            ;; Express length & orientation of grid segments in the
            ;; xi and eta directions as complex vectors.
            cx = mgh_stagger(mgh_diff(x+im*y, 1), DELTA=[1,0])
            cy = mgh_stagger(mgh_diff(x+im*y, 2), DELTA=[0,1])
            pm = one/abs(cx)
            pn = one/abs(cy)
            angle = atan((cx-im*cy)/2, /PHASE)
            mgh_undefine, cx, cy
         endelse

      end

      1: begin

         message, 'Calculating grid spacing & angle from lon & lat assuming a spherical Earth', /INFORM

         if keyword_set(grid) then begin
            pm = one/(radius*mgh_stagger(mgh_diff(x)*pi/180, DELTA=1)#cos(pi*y/180))
            pn = one/mgh_inflate(dim, radius*mgh_stagger(mgh_diff(y), DELTA=1)*pi/180, 2)
            angle = make_array(dim, SIZE=rtype)
            message, 'This option not completely implemented!'
         endif else begin
            ;; Range & azimuth along line segments in the xi
            ;; direction:
            rx = make_array(dim-[1,0], TYPE=rtype)
            ax = make_array(dim-[1,0], TYPE=rtype)
            for j=0,d1-1 do begin
               for i=0,d0-2 do begin
                  ra = map_2points(x[i,j], y[i,j], x[i+1,j], y[i+1,j])
                  rx[i,j] = ra[0]*radius*pi/180
                  ax[i,j] = ra[1]
               endfor
            endfor
            ;; Range & azimuth along line segments in the eta
            ;; direction:
            ry = make_array(dim-[0,1], TYPE=rtype)
            ay = make_array(dim-[0,1], TYPE=rtype)
            for j=0,d1-2 do begin
               for i=0,d0-1 do begin
                  ra = map_2points(x[i,j], y[i,j], x[i,j+1], y[i,j+1])
                  ry[i,j] = ra[0]*radius*pi/180
                  ay[i,j] = ra[1]
               endfor
            endfor
            ;; Quantities required ar pm, pn and angle
            pm = mgh_stagger(one/rx, DELTA=[1,0])
            pn = mgh_stagger(one/ry, DELTA=[0,1])
            cx = mgh_stagger(exp(im*(90-ax)*pi/180), DELTA=[1,0])
            cy = mgh_stagger(exp(-im*ay*pi/180), DELTA=[0,1])
            angle = atan((cx+cy)/2, /PHASE)
            mgh_undefine, rx, ry, ax, ay, cx, cy
         endelse

      end

      2: begin

         message, 'Calculating grid spacing & angle from lon & lat assuming an ellipsoidal Earth', /INFORM

         if keyword_set(grid) then begin
            ;; Range along line segments in the xi direction:
            rx = make_array(dim-[1,0], TYPE=rtype)
            for j=0,d1-1 do begin
               for i=0,d0-2 do begin
                  eqldaz, x[i], y[j], x[i+1], y[j], void, dd
                  rx[i,j] = dd
               endfor
            endfor
            ;; Range along line segments in the eta direction:
            ry = make_array(dim-[0,1], TYPE=rtype)
            for j=0,d1-2 do begin
               for i=0,d0-1 do begin
                  eqldaz, x[i], y[j], x[i], y[j+1], void, dd
                  ry[i,j] = dd
               endfor
            endfor
            pm = mgh_stagger(one/rx, DELTA=[1,0])
            pn = mgh_stagger(one/ry, DELTA=[0,1])
            angle = make_array(dim, SIZE=rtype)
            message, 'This option not completely implemented!'
         endif else begin
            ;; Range & azimuth along line segments in the xi
            ;; direction:
            rx = make_array(dim-[1,0], TYPE=rtype)
            ax = make_array(dim-[1,0], TYPE=rtype)
            for j=0,d1-1 do begin
               for i=0,d0-2 do begin
                  eqldaz, x[i,j], y[i,j], x[i+1,j], y[i+1,j], void, dd, az
                  rx[i,j] = dd
                  ax[i,j] = az
               endfor
            endfor
            ;; Range & azimuth along line segments in the eta
            ;; direction:
            ry = make_array(dim-[0,1], TYPE=rtype)
            ay = make_array(dim-[0,1], TYPE=rtype)
            for j=0,d1-2 do begin
               for i=0,d0-1 do begin
                  eqldaz, x[i,j], y[i,j], x[i,j+1], y[i,j+1], void, dd, az
                  ry[i,j] = dd
                  ay[i,j] = az
               endfor
            endfor
            pm = mgh_stagger(one/rx, DELTA=[1,0])
            pn = mgh_stagger(one/ry, DELTA=[0,1])
            cx = mgh_stagger(exp(im*(90-ax)*pi/180), DELTA=[1,0])
            cy = mgh_stagger(exp(-im*ay*pi/180), DELTA=[0,1])
            angle = atan((cx+cy)/2, /PHASE)
         endelse

      end

   endcase

   ;; Calculate derivatives and interpolate back to the main grid.

   if arg_present(dndx) then $
        dndx = mgh_stagger(mgh_diff(one/pn, 1), DELTA=[1,0])

   if arg_present(dmde) then $
        dmde = mgh_stagger(mgh_diff(one/pm, 2), DELTA=[0,1])

end
