;+
; NAME:
;   MGH_SAN_SEARCH
;
; PURPOSE:
;   This function searches for files matching a specified pattern in a specified
;   directory on a remote data store (SAN). It searches for the unsuffixed pattern
;   and for patterns with the added suffixes '.bz2' and '.gz'. The result is always
;   unsuffixed and duplicates are removed and the result sorted.
;
; CATEGORY:
;   Miscellaneous.
;
; CALLING SEQUENCE:
;   result = mgh_san_search(pattern, SUBDIRECTORY=subdirectory, COUNT=count, $
;                           MIRROR=mirror, VOLUME=volume)
;
; POSITIONAL PARAMETERS:
;   pattern (input, string, optional)
;     A string, which may include wild-card characters, specifying the
;     path to be matched.
;
; KEYWORD PARAMETERS:
;   COUNT (output, integer scalar)
;     This keyword returns the number of files found.
;
;   MIRROR (input, switch)
;     If set, search on the mirror, if not search in the master SAN area.
;     Default is 0.
;
;   SUBDIRECTORY (input, string)
;     The SAN subdirectories
;
;   VOLUME (input, scalar string)
;     The SAN volume index.
;
;   Unrecognised keywords are passed to FILE_SEARCH.
;
; RETURN VALUE:
;   The function returns a string array containing the names of the files
;   matching the pattern. The base part of each file name (i.e. the part
;   formed with the SUBDIRECTORY, MIRROR and VOLUME keywords) is stripped off.
;
; DEPENDENCIES:
;   Requires MGH_SAN_PATH to form path names. This requires a system variable,
;   !MGH_SAN, to be defined, as described in the documentation.
;
; MODIFICATION HISTORY:
;   Mark Hadfield, 2008-09:
;     Written.
;   Mark Hadfield, 2011-08:
;     - Added GUNZIP keyword.
;   Mark Hadfield, 2011-08:
;     - GUNZIP keyword renamed GZIP and option 2 added, provided
;       automagical handling of remote gzipped files.
;   Mark Hadfield, 2012-10:
;     - Changed the method for forming the result when FULL_PATH is not
;       defined. The new method handles the case where SUBDIRECTORY
;       contains '..' entries. However it may break in other situations.
;   Mark Hadfield, 2012-11:
;     - The previous month's fix did break in certain situations, eg.
;       in handling CCMP files. Now trying an alternative.
;   Mark Hadfield, 2015-07:
;     - Removed support for the FULL_PATH keyword, as it hasn't been
;       used for some time.
;     - Removed references to the GZIP keyword, as it was rendered
;       unnecessary and non-functional during the 2011 modifications.
;     - Removed keyword inheritance by FILE_SEARCH. Though there may be
;       good reasons for passing extra keywords to FILE_SEARCH, I haven't
;       run into any yet, and when I do I will treat them on a case-by-case
;       basis.
;     - Implemented support for access to remote volumes by SSH.
;-
function mgh_san_search, name, $
     COUNT=count, MIRROR=mirror, $
     SUBDIRECTORY=subdirectory, VOLUME=volume

   compile_opt DEFINT32
   compile_opt STRICTARR
   compile_opt STRICTARRSUBS
   compile_opt LOGICAL_PREDICATE

   if n_elements(name) eq 0 then name = '*'

   if n_elements(volume) eq 0 then volume = '/hpcf/working/hadfield'

   if n_elements(mirror) eq 0 then mirror = !false

   if ~ isa(volume, 'STRING') then $
      message, BLOCK='mgh_mblk_motley', NAME='mgh_m_wrongtype', 'volume'

   if keyword_set(mirror) then begin
      access = 'mounted'
   endif else begin
      vol = (!mgh_san)[volume]
      access = mgh_struct_has_tag(vol, 'access') ? vol.access : 'mounted'
   endelse

   pattern = mgh_san_path(name+['','.gz','.bz2'], SUBDIRECTORY=subdirectory, MIRROR=mirror, VOLUME=volume)

   case access of
      'mounted': $
         match = file_search(pattern, COUNT=count)
      'ssh': $
         match = mgh_ssh_file_search(pattern, COUNT=count)
   endcase

   if count eq 0 then return, ''

   for i=0,count-1 do begin
      case !true of
         strmatch(match[i], '*.bz2'): $
            match[i] = strmid(match[i], 0, strlen(match[i])-4)
         strmatch(match[i], '*.gz'): $
            match[i] = strmid(match[i], 0, strlen(match[i])-3)
         else:
      endcase
   endfor

   match = match[uniq(match, sort(match))]

   count = n_elements(match)

   ;; Strip off the base part of the path plus the separator. The
   ;; file_search function is used here to normalise the base path.

   base = mgh_san_path('', SUBDIRECTORY=subdirectory, MIRROR=mirror, VOLUME=volume)
   case access of
      'mounted': $
         return, strmid(match, strlen(file_search(base))+1)
      'ssh': $
         return, strmid(match, strlen(mgh_ssh_file_search(base)))
   endcase

end
