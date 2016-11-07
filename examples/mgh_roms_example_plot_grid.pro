;+
; NAME:
;   MGH_ROMS_EXAMPLE_PLOT_GRID
;
; PURPOSE:
;   Grid plot examples
;
;###########################################################################
; Copyright (c) 2016 NIWA:
;   http://www.niwa.co.nz/
; Licensed under the MIT open source license:
;   http://www.opensource.org/licenses/mit-license.php
;###########################################################################
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2016-11:
;     Written.
;-
pro mgh_roms_example_plot_grid, option

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(option) eq 0 then option = 0

   root = file_dirname(routine_filepath('mgh_roms_example_plot_grid'))

   case option of

      0: begin

         name = 'WC13 depth'

         file = filepath('wc13_grd.nc', ROOT=root, SUBDIR=['data','WC13','Data'])

         mgh_new, 'mgh_roms_plot_grid', file, $
            GRAPH_PROPERTIES={name: name}

      end

      1: begin

         name = 'WC13 Coriolis parameter'

         file = filepath('wc13_grd.nc', ROOT=root, SUBDIR=['data','WC13','Data'])

         mgh_new, 'mgh_roms_plot_grid', file, $
            VARIABLE='f', DATA_MULTIPLIER=1E4, $
            GRAPH_PROPERTIES={name: name}, $
            COLORBAR_PROPERTIES={title: 'f (10^-4 s^-1)'}

      end

      2: begin

         name = 'WC13 inverse xi spacing'

         file = filepath('wc13_grd.nc', ROOT=root, SUBDIR=['data','WC13','Data'])

         mgh_new, 'mgh_roms_plot_grid', file, $
            VARIABLE='pm', DATA_MULTIPLIER=1E3, $
            GRAPH_PROPERTIES={name: name}, $
            COLORBAR_PROPERTIES={title: 'pm (km^-1)'}

      end

      3: begin

         name = 'WC13 inverse eta spacing'

         file = filepath('wc13_grd.nc', ROOT=root, SUBDIR=['data','WC13','Data'])

         mgh_new, 'mgh_roms_plot_grid', file, $
            VARIABLE='pn', DATA_MULTIPLIER=1E3, $
            GRAPH_PROPERTIES={name: name}, $
            COLORBAR_PROPERTIES={title: 'pn (km^-1)'}

      end

   endcase


end
