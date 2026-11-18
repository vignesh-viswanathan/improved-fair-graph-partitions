"""Exact checker for certificate.json. Run: python check.py [certificate.json]."""
import json, sys
if not __debug__: raise RuntimeError('Run without -O: assertions are proof checks.')
from pathlib import Path
from itertools import combinations_with_replacement, product
from math import comb
import sympy as s

A, E = s.Rational(38, 25), s.Rational(1, 240)
a, b, c, t, u = s.symbols('a b c t u')
X, S, T = (1,1,1,1,0,0,0,0), (1,1,0,0,1,1,0,0), (1,0,1,0,1,0,1,0)
ROOTS = {
    'B': ((a,b), ((s.Rational(1,2)-E,1),(0,s.Rational(1,2)+E)), ((False,False),(True,False)), ('B',)),
    'large': ((t,), ((0,1),), ((True,False),), ('large',)),
    'high': ((t,), ((s.Rational(1144,1575),1),), ((False,False),), ('high',)),
    'exception': ((c,u), ((0,8*E),(0,1)), ((True,False),(False,False)), ('exception',)),
    'low': ((c,t,u), ((0,s.Rational(1,2)),(0,s.Rational(1144,1575)),(0,1)), ((True,False),(True,False),(False,False)), ('small','tiny'))}

def rows(kind):
    pairs = list(combinations_with_replacement((1,2,3,4,5,6,7) if kind=='B' else (0,1,2,3,5,6,7), 2))
    def edge(i,j):
        return s.Matrix([[int(p==tuple(sorted((i,j)))) for p in pairs]])
    def gain(i,Q,R):
        return s.Matrix([[sum(Q[w]-A*int(R[v]==R[w]) for v,w in ((j,k),(k,j)) if v==i) for j,k in pairs]])
    def degree(i):
        return s.Matrix([[sum(2*int(X[v]==X[w])-1 for v,w in ((j,k),(k,j)) if v==i) for j,k in pairs]])
    def cut(q):
        return s.Matrix([[(2*int(X[i]==X[j])-1)*(q[i]+q[j]-2*q[i]*q[j]) for i,j in pairs]])
    def local(i,j,z):
        return sizes[j]*degree(i)+sizes[i]*degree(j)+2*z*edge(i,j)
    overlap = 4 if kind=='B' else 0
    J = gain(overlap,S,X)+gain(overlap,T,S)
    if kind=='B':
        sizes = [0,a/2,0,0,b/2,(1-a-b)/2,0,0]
        rr = [gain(i,S,X) for i in (1,4,5)] + [local(1,4,E),local(1,5,E),cut([int(x!=y) for x,y in zip(X,S)])]
    elif kind in ('large','high'):
        if kind=='large':
            ll = degree(5)+t*degree(0)+s.Rational(2,3)*t*edge(0,5)
            cc = cut([t,0,0,0,0,1,0,0])+s.Rational(4,3)*t*t*(1-t)*edge(0,0)
        else:
            ll = degree(0)+t*degree(5)+s.Rational(2,5)*edge(0,5)
            cc = cut([1,0,0,0,0,t,0,0])+t*(1-t)*edge(5,5)/2
        rr = [gain(0,S,X),gain(5,S,X),ll,cc]
    else:
        bb = 3*c/4 if kind=='exception' else t*c
        r = bb*u
        sizes = [bb/2,(1-c-bb)/2,(c-r)/2,r/2,0,c/2,(1-c-bb+r)/2,(bb-r)/2]
        z = c/8 if kind=='exception' else E if kind=='small' else bb/2
        q = [1,(c-bb)*(s.Rational(1,2)-c)/(1-c-bb),(c-bb)*(s.Rational(1,2)+c)/(c-r),0,0,1,0,0]
        cc = cut(q)+z*(s.Rational(17,8)*q[1]*(1-q[1])*edge(1,1)/sizes[1]+edge(2,2)/sizes[2])+z*z*edge(1,2)/(2*sizes[1]*sizes[2])
        rr = [gain(i,S,X) for i in (0,1,5)]+[gain(i,T,S) for i in (0,2,6)]
        rr += [local(i,j,z) for i,j in ((0,5),(0,6),(1,5),(1,6))]+[cc]
    return s.Matrix.vstack(*rr), J

def polynomial(terms,vs):
    assert all(len(ex)==len(vs) and all(type(k) is int and k>=0 for k in ex) for _,ex in terms)
    return sum(s.Integer(co)*s.prod(v**k for v,k in zip(vs,ex)) for co,ex in terms)

def value(rep,polys):
    q, ff = rep
    assert type(q) is str and all(type(i) is int and 0<=i<len(polys) and type(e) is int for i,e in ff)
    return s.Rational(q)*s.prod(polys[i]**e for i,e in ff)

def bernstein(p,vs,box):
    n = s.Poly(p,*vs).degree_list()
    pp = s.Poly(p.subs({v:lo+(hi-lo)*v for v,(lo,hi) in zip(vs,box)}, simultaneous=True),*vs)
    return {i:sum(co*s.prod(s.Rational(comb(ii,k),comb(nn,k)) for ii,k,nn in zip(i,ex,n))
        for ex,co in pp.terms() if all(k<=ii for k,ii in zip(ex,i))) for i in product(*(range(nn+1) for nn in n))}

def halve(coeffs,axis):
    n = max(i[axis] for i in coeffs)
    left, right = {}, {}
    for fixed in (i for i in coeffs if i[axis]==0):
        keys = [fixed[:axis]+(j,)+fixed[axis+1:] for j in range(n+1)]
        row = [coeffs[k] for k in keys]
        left[keys[0]], right[keys[n]] = row[0], row[-1]
        for j in range(1,n+1):
            row = [(x+y)/s.Integer(2) for x,y in zip(row,row[1:])]
            left[keys[j]], right[keys[n-j]] = row[0], row[-1]
    return left, right

def cover(node,coeffs,opens,duals):
    if len(node)==3:
        axis, left, right = node
        assert type(axis) is int and 0<=axis<len(opens)
        halves = [halve(p,axis) for p in coeffs]
        ol, orr = list(opens), list(opens)
        ol[axis], orr[axis] = (opens[axis][0],False), (False,opens[axis][1])
        return cover(left,[p[0] for p in halves],ol,duals)+cover(right,[p[1] for p in halves],orr,duals)
    assert len(node)==2 and type(node[0]) is int and 0<=node[0]<len(duals)
    reps = duals[node[0]]['multipliers']+duals[node[0]]['residuals']
    signs = dict(node[1])
    needed = {i for q,ff in reps if s.Rational(q) for i,e in ff}
    denominators = {i for q,ff in reps if s.Rational(q) for i,e in ff if e<0}
    for i in needed:
        assert signs[i] in (-1,1)
        cc = {k:signs[i]*v for k,v in coeffs[i].items()}
        assert all(v>=0 for v in cc.values())
        if i in denominators:
            n = tuple(max(k[j] for k in cc) for j in range(len(opens)))
            faces = [[i for i,o in ((0,op[0]),(nn,op[1])) if not o] or [None] for nn,op in zip(n,opens)]
            assert all(any(v>0 for k,v in cc.items() if all(f is None or f==j for f,j in zip(face,k))) for face in product(*faces))
    assert all(not s.Rational(q) or s.Rational(q)*s.prod(signs[i]**(e%2) for i,e in ff)>=0 for q,ff in reps)
    return 1

def verify(path):
    data = json.loads(Path(path).read_text())
    assert set(data)==set(ROOTS)
    total = 0
    for name,(vs,box,opens,kinds) in ROOTS.items():
        part = data[name]
        polys = [polynomial(p,vs) for p in part['polynomials']]
        for d in part['duals']:
            assert d['model'] in kinds and d['target'] in ('overlap','sum')
            M,J = rows(d['model'])
            if d['target']=='sum':
                assert name in ('exception','low')
                J = sum((M[i,:] for i in range(6)),s.zeros(1,M.cols))
            lam,mu = [s.Matrix([[value(r,polys) for r in d[k]]]) for k in ('multipliers','residuals')]
            assert lam.shape==(1,M.rows) and mu.shape==J.shape
            assert all(s.cancel(v)==0 for v in J+lam*M+mu)
        count = cover(part['cover'],[bernstein(p,vs,box) for p in polys],opens,part['duals'])
        total += count
        print(f'{name}: identities and signs verified on {count} box(es)',flush=True)
    print(f'PASS: exact certificate for alpha=38/25, augmented N>={1/E}; {total} boxes.')

if __name__=='__main__':
    verify(sys.argv[1] if len(sys.argv)>1 else Path(__file__).with_name('certificate.json'))
