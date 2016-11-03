;+
; NAME:
;   MGH_ROMS_EXAMPLE_PLOT_GRID
;
; PURPOSE:
;   Lake Jersey grid plot example
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

         name = 'INLET_TEST grid'

         file = filepath('inlet_test_grid.nc', ROOT=root, SUBDIR=['data','inlet_test','Data'])

         mgh_new, 'mgh_roms_plot_grid', file, $
            GRAPH_PROPERTIES={name: name}

      end

   endcase


end
