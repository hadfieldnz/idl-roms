;+
; NAME:
;   MGH_ROMS_FILE
;
; PURPOSE:
;   Create a ROMS file-sequence object for a specified subdirectory of
;   the SAN area. By default the function searches for ROMS history
;   files and returns them wrapped in a history-file object.
;
; CALLING SEQUENCE:
;   object = MGH_ROMS_FILE(subdirectory)
;   object = MGH_ROMS_FILE(subdirectory, TYPE=type)
;
; POSITIONAL ARGUMENTS:
;   subdirectory (input, string array)
;     Subdirectory in which to search for ROMS output files.
;
; KEYWORD ARGUMENTS:
;   FILE_RANGE (input, 2-element integer)
;     Specify that the file-sequence object should wrap a subset of
;     the files found.
;
;   MTIME (output, long or long64 scalar)
;     If this argument is present, return the latest last-modified
;     time stamp for the files in the sequence.
;
;   NAME (input, scalar string)
;     Name for the file-sequence object. Default is created from
;     subdirectory and type.
;
;   PATTERN (input, scalar string)
;     Search pattern. The default depends on the value of TYPE.
;
;   TYPE (input, scalar string)
;     Specify the file type
;
; RETURN VALUE:
;   The function returns a single object reference.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2002-05:
;     Written.
;   Mark Hadfield, 2002-07:
;     File type now controlled by TYPE keyword, not by a separate
;     switch-type keyword for each type (FLOAT, STATION etc).
;   Mark Hadfield, 2011-08:
;     - Added GUNZIP and MTIME keywords.
;   Mark Hadfield, 2011-08:
;     - GUNZIP keyword renamed GZIP; default is now 2, to take advantage of
;       new functionality in MGH_SAN_FILE.
;   Mark Hadfield, 2012-01:
;     - GZIP keyword deleted to match MGH_SAN_FILE.
;   Mark Hadfield, 2014-07:
;     - Re-indented source.
;   Mark Hadfield, 2015-06:
;     - Made the object name more detailed to avoid ambiguity in temporary file names.
;   Mark Hadfield, 2016-02:
;     Not a modification, just a comment: recently I have been using the hashcode of the
;     list of files to construct temporary file names, so the modification of 2015-06
;     is redundant (but not harmful as far as I can see.)
;   Mark Hadfield, 2016-08:
;     - Changed the default value of PATTERN for each type to something
;       more restrictive (eg. roms_avg_[0-9][0-9][0-9][0-9].nc) to avoid
;       unintended inclusions.
;   Mark Hadfield, 2019-04:
;     - Added a "regrid" type associated with MGHromsRegrid objects.
;-
function mgh_roms_file, subdirectory, $
     FILE_CLASS=file_class, FILE_RANGE=file_range, $
     MIRROR=mirror, MTIME=mtime, NAME=name, PATTERN=pattern, TYPE=type, $
     VOLUME=volume

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   on_error, 2

   if n_elements(type) eq 0 then type = 'history'

   case strlowcase(type) of

      'average': begin
         if n_elements(pattern) eq 0 then pattern = 'roms_avg_[0-9][0-9][0-9][0-9].nc'
         if n_elements(file_class) eq 0 then file_class = 'MGHromsHistory'
      end

      'climate': begin
         if n_elements(pattern) eq 0 then pattern = 'roms_clm_[0-9][0-9][0-9][0-9].nc'
         if n_elements(file_class) eq 0 then file_class = 'MGHromsHistory'
      end

      'forcing': begin
         if n_elements(pattern) eq 0 then pattern = 'roms_frc_[0-9][0-9][0-9][0-9].nc'
         if n_elements(file_class) eq 0 then file_class = 'MGHromsHistory'
      end

      'float': begin
         if n_elements(pattern) eq 0 then pattern = 'roms_flt_[0-9][0-9][0-9][0-9].nc'
         if n_elements(file_class) eq 0 then file_class = 'MGHromsFloat'
      end

      'grid': begin
         if n_elements(pattern) eq 0 then pattern = 'roms_grd.nc'
         if n_elements(file_class) eq 0 then file_class = 'MGHromsHistory'
      end

      'history': begin
         if n_elements(pattern) eq 0 then pattern = 'roms_his_[0-9][0-9][0-9][0-9].nc'
         if n_elements(file_class) eq 0 then file_class = 'MGHromsHistory'
      end

      'restart': begin
         if n_elements(pattern) eq 0 then pattern = 'roms_rst_[0-9][0-9][0-9][0-9].nc'
         if n_elements(file_class) eq 0 then file_class = 'MGHromsHistory'
      end

      'station': begin
         if n_elements(pattern) eq 0 then pattern = 'roms_sta_[0-9][0-9][0-9][0-9].nc'
         if n_elements(file_class) eq 0 then file_class = 'MGHromsStation'
      end

      'boundary': begin
         if n_elements(pattern) eq 0 then pattern = 'roms_bry*.nc'
         if n_elements(file_class) eq 0 then file_class = 'MGHncSequence'
      end

      'regrid': begin
         if n_elements(pattern) eq 0 then pattern = 'roms_frc_[0-9][0-9][0-9][0-9].nc'
         if n_elements(file_class) eq 0 then file_class = 'MGHromsRegrid'
      end

      'miscellaneous': begin
         ;; No default pattern
         if n_elements(file_class) eq 0 then file_class = 'MGHncSequence'
      end

   endcase

   match = mgh_san_search(pattern, COUNT=n_file, SUBDIRECTORY=subdirectory, VOLUME=volume)

   if n_file eq 0 then begin
      dir = mgh_san_path(SUBDIRECTORY=subdirectory, VOLUME=volume)
      message, 'No match for pattern '+pattern+' in directory '+dir
   endif

   file = mgh_san_file(match, MTIME=mtime, SUBDIRECTORY=subdirectory, VOLUME=volume)

   mgh_resolve_indices, n_file, file_range, 1, indices

   result = obj_new(file_class, file)

   ;; Create object name from subdirectory, type, pattern & file range. Note
   ;; that the subdirectory argument may be an array or a "/"-delimited scalar string

   if n_elements(name) eq 0 then begin
      name = mgh_str_subst(strjoin(subdirectory,'-'), '/', '-')+'_type:'+strlowcase(type)
      if n_elements(pattern) gt 0 then $
         name = name + '_pattern:' + mgh_str_vanilla(pattern)
      if n_elements(file_range) gt 0 then $
         name = name + '_range:' + strjoin(strtrim(file_range, 2), '-')
   endif

   result->SetProperty, NAME=name

   return, result

end
