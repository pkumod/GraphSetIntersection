import sys

p = 0.25
w = 32
b = 8

p_tn = (2 ** w - 2 ** (w - b)) / (2 ** w - 1.0) * (1.0 - p)
p_fp = (2 ** (w - b) - 1.0) / (2 ** w - 1.0) * (1.0 - p)

p_nohit = ((2.0 ** w - 2.0 ** (w - b)) / (2.0 ** w - 1.0) * (1.0 - p)) ** 4
p_onehit = 4.0 * ((2.0 ** w - 2.0 ** (w - b)) / (2.0 ** w - 1.0) * (1.0 - p)) ** 3 \
    * (p + (2.0 ** (w - b) - 1.0) / (2.0 ** w - 1.0) * (1.0 - p))

p_mulhit = 1.0 - p_nohit - p_onehit

p_mulmatch = 1 - (p_nohit + p_onehit) ** 4
p_nomatch = p_nohit ** 4

e_bmiss = p_nomatch * 4 + (1.0 - p_nomatch) * 12
e_qfilter = p_nomatch * 2 + p_mulmatch * (2 + 2 * 1.1 + 2 + 2) \
 + (1.0 - p_nomatch - p_mulmatch) * (2 + 2 + 2)

print 'p_tn =', p_tn
print 'p_fp =', p_fp
print 'p_nohit =', p_nohit
print 'p_onehit =', p_onehit
print 'p_mulhit =', p_mulhit
print 'p_mulmatch =', p_mulmatch
print 'p_nomatch =', p_nomatch
print 'e_bmiss =', e_bmiss
print 'e_qfilter =', e_qfilter