FUNCTION lower_limit_gehrels, cl, n
; compute approximate expression of the upper limit at the cl 
; confidence level for n counts using formula (12) of 
; Gehrels 1986, ApJ, 303, 336
;
; Nicolas Grosso, August 2003
;
;--------------------------------------------------------------------------

y_a=GAUSS_CVF(cl)
S=ABS(y_a)
ll=n*(1.-1./(9.*n) - S/(3*sqrt(n)))^3
RETURN, ll
END