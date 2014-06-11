;+
; NAME:
;    ADD
;
; PURPOSE:
;    Adds an element to (possibly non-existent) array.
;
; CATEGORY:
;    Misc
;
; CALLING SEQUENCE:
;    PUSH, Array, Value
;
; INPUTS:
;    Array:   An array to which the values are to be added. The array
;             is created if it doesn't already exist.
;
;    Value:   A value to be added to the end of Array. Must be the
;             same type as Array.
;
pro add, array, value

case size(array, /n_dimen) of
   0: array = value
   1: array = array + value
   2: array = array + value
   else: message, 'Array cannot have more than 1 dimension.'
endcase

end
