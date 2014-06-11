PRO merge_vector, v1, v2

; merge two vectors (1-D) into one
; v1 will be replaced with the new merged vector
; 



n1 = n_elements(v1)
n2 = n_elements(v2)

v3 = make_array(n1+n2, type=size(v1,/type))

v3[0:n1-1] = v1
v3[n1:n1+n2-1] = v2
v1 = v3


END