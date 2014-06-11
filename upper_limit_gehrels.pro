FUNCTION upper_limit_gehrels,cl,n
; compute approximate expression of the upper limit at the cl 
; confidence level for n counts using formula (10) of 
; Gehrels 1986, ApJ, 303, 336
;
; Nicolas Grosso, August 2003.
;
;--------------------------------------------------------------------------

y_a=GAUSS_CVF(cl)
S=ABS(y_a)
ul=n+S*SQRT(n+1.)+(S^2+2.)/3.
RETURN,ul
END


