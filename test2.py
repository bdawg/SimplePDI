def testme(a, answer=answer):
    answer = a*2

a = None
try:
    a
except:
    print 'Nope'
else:
    print 'Yep'