"""
Algebraic feedback shift register sequences,
N-adic sequences, Boolean monotone and bent functions. 
Also, we include numerous examples of the 
sage.crypto.boolean_function module.


AUTHORS:
   David Joyner, wdjoyner@gmail.com
   Charles Celerier, charles.celerier@gmail.com
   copyright 2011-2012

LICENSE:
   Modified BSD

last modified 2012-04-10

#####

Sage implements in the module sage.crypto.boolean_function
several classes and methods for boolean functions.

We recall some background from [1].

A function $f$ on $GF(2)^n$ can be uniquely represented by
a polynomial on $GF(2)$ whose degree in each variable in each term is 
at most $1$. Namely,

\[
f (x_1 ,\dots , x_n ) =
\sum_{(a_1,\dots ,a_n)\in GF(2)^n}
g(a_1 ,\dots , a_n )x^{a_1}\cdot \dots \cdot x^{a_n},
\]
where $g$ is also a function on $GF(2)^n$. This polynomial 
representation of $f$ is called
the {\it algebraic normal form} (ANF) of the function.
The {\it algebraic degree} of $f$, denoted by $\deg(f )$,
is defined as the number of variables in the longest term 
of the ANF of $f$.


A Boolean function $f$ is said to be {\it balanced} if 
${\rm wt}(f ) = {\rm wt}(f + 1) = 2^{n-1}$. (Here $1$ denotes the 
constant function $1$ on $GF(2)^n$.
A {\it subfunction} of a Boolean function $f$ is a function $f$ 
obtained by substituting some constants for some variables in $f$.

The Boolean function $f$ is called 
{\it correlation immune of order $m$} if 
$[\rm wt}(f') = {\rm wt}(f )/2^m$ for any its subfunctions 
$f'$ of $n - m$ variables. 



If $\deg(f ) \leq 1$ then $f$ is called an {\it affine function}. 
If $f$ is an affine function and $f(0) = 0$ then $f$ is called a 
{\it linear function}.

For two Boolean functions $f_1$ and $f_2$ on $GF(2)^n$, 
we define the {\it distance} between $f_1$ and $f_2$ by 

\[
d(f_1 , f_2 ) = | \{x \in GF(2)^n \ |\ f_1 (x) \not= f_2 (x)\} |.
\]
(Here $|S|$ denotes the cardinality of a set $S$.)
The {\it weight} of a Boolean function is the distance between it
and the constant function $0$, denoted ${\rm wt}(f)$.
It is easy to see that 

\[
d(f_1,f_2)={\rm wt}(f_1+f_2).
\]
The minimum distance between $f$ and the set of all affine functions 
is called the {\it nonlinearity} of $f$ and denoted by 
$[\rm nl}(f )$.



The {\it Walsh transform} of a Boolean function $f$ is an integer-
valued function over $GF(2)^n$ that can be defined as

\[
W_f(u) =
\sum_{x in GF(2)^n}
(-1)^{f(x)+ \langle u,x\rangle}.
\]
The Walsh coefficients satisfy {\it Parseval’s equation}

\[
\sum_{x in GF(2)^n}
W_f(u)^2 = 2^{2n}.
\]


For any Boolean function $f$ on $GF(2)^n$, we have

\[
wt(f ) = 2^{n-1} - {\frac{1}{2}W_f(0),
\]
and
\[
{\rm nl}(f ) = 2^{n-1} - {\frac{1}{2}\max_{u \in GF(2)^n} |W_f (u)|.
\]


For each $u \in GF(2)^n$, the {\it autocorrelation coefficient} of the function 
$f$ at the vector $u$ is defined 

\[
\Delta_f(u) =
\sum_{x in GF(2)^n}
(-1)^{f(x)+ f(x+u)}.
\]

The {\it absolute indicator} of $f$ is defined as

\[
\Delta_f = \max_{u \in GF(2)^n-\{0\}} |\Delta_f (u)|.
\]


If f is any Boolean function having the property
that it is supported on vectors having even weight then 
we say f is an {\it even} Boolean function.

Conjecture: If f is a monotone even Boolean function then
there is an x in GF(2)^n such that W(f)(x) = 0.

For example, this has been checked for all the 168 monotone Boolean functions 
of 4 variables.


EXAMPLES:
    sage: from sage.crypto.boolean_function import *
    sage: P.<x0,x1,x2,x3> = BooleanPolynomialRing()
    sage: b = x0*x1 + x2*x3
    sage: f = BooleanFunction(b)
    sage: [b(x[0],x[1],x[2],x[3]) for x in GF(2)^4]
    [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0]
    sage: f.truth_table()
    (False, False, False, True, False, False, False, True, False, False, 
     False, True, True, True, True, False)
    sage: WT = f.walsh_hadamard_transform(); WT
    (-4, -4, -4, 4, -4, -4, -4, 4, -4, -4, -4, 4, 4, 4, 4, -4)
    sage: f.absolute_walsh_spectrum()
    {4: 16}
    sage: f.nonlinearity()
    6
    sage: 2^(4-1) - (1/2)*max([abs(x) for x in WT])
    6
    sage: f.autocorrelation()
    (16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    sage: f.absolute_autocorrelation()
    {16: 1, 0: 15}
    sage: f.absolut_indicator() # this is a mis-spelling in Sage
    0
    sage: f.is_bent()
    True
    sage: f.is_balanced()
    False
    sage: f.is_symmetric()
    False
    sage: f.sum_of_square_indicator()
    256
    sage: f.correlation_immunity()
    0


REFERENCES:
  [1] Yuriy Tarannikov, Peter Korolev, and Anton Botev, 
      {\it Autocorrelation coefficients and correlation immunity of Boolean 
      functions}, in {\bf Advances in Cryptology—ASIACRYPT 2001}. Available:
      \newline
      \url{http://www.iacr.org/cryptodb/data/paper.php?pubkey=515}

"""

######## utility functions #########################

def GFn_to_integer(v):
    """
    Returns the integer value of the vector v interpreted as a binary number

    INPUT:
        v - vector in GF(2)^n

    EXAMPLES:
        sage: V = GF(2)^2
        sage: [GFn_to_integer(v) for v in V]

    """
    S = 0
    for i in range(len(v)):
        S+=2^i*ZZ(v[i])
    return S

def Hamming_weight(v):
    """
    Returns the Hamming weight of the vector v

    INPUT:
        v - vector in GF(2)^n

    EXAMPLES:
        sage: V=GF(2)^5
        sage: a=V([0,1,1,0,1]); b=V([1,1,1,0,0]); c=V([0,0,1,1,0])
        sage: Hammin
        HammingCode     Hamming_weight  
        sage: Hamming_weight(a)
        3
        sage: Hamming_weight(b)
        3
        sage: Hammi
        HammingCode     Hamming_weight  
        sage: Hamming_weight(c)
        2

    """
    return len(v.support())


###################################################


def connection_element(Qs, p):
    """
    Input consists of 
        Qs - r+1 elements from T ("taps")
        p - prime

    This function returns the connection element

    EXAMPLES:
        sage: p = 5; Qs = [2,2,1,2]
        sage: connection_element(Qs, p)
        283

    """
    r = len(Qs)-1
    c = -Qs[0]+sum([Qs[i]*p**i for i in range(1,r+1)])
    return c


def generating_quotient(As, Qs, m, p):
    """
    Input consists of 
        S, T - complete sets of reps mod p
        As - r elements from S
        Qs - r+1 elements from T
        m - integer
        p - prime
        N - positive integer > r

    This function returns a rational whose p-adic expansion 
    agrees with As.

    EXAMPLES:
        sage: p = 5; Qs = [2,2,1,2]; As = [1,2,2]; m = 4
        sage: generating_quotient(As, Qs, m, p)
        -487/283
        sage: alpha = Qp(5)(-487/283)
        sage: alpha.padded_list(20)
        [1, 2, 2, 1, 0, 1, 0, 1, 2, 0, 0, 0, 3, 0, 2, 0, 2, 4, 0, 0]

    This gives the same output as afsr does.

    """
    r = len(Qs)-1
    DS = sum([sum([Qs[i]*As[n-i]*p**n for i in range(1,n+1)]) for n in range(1,r)])
    u = -m*p**r - Qs[0]*sum([As[i]*p**i for i in range(r)]) + DS
    return u/connection_element(Qs, p)


def afsr(S, T, As, Qs, m, p, N):
    """
    Input consists of 
        S, T - complete sets of reps mod p
        As - r elements from S
        Qs - r+1 elements from T
        m - integer
        p - prime
        N - positive integer > r

    This function returns a list of length N generated
    by the AFSR associated to (S, T, As, Qs, m, p), as in [KX]
  
    EXAMPLES:
        sage: p = 5; N = 20; S = range(p); T = S
        sage: As = [1,2,2]; Qs = [1,2,1,2]; m = 4
        sage: afsr(S, T, As, Qs, m, p, N)
        [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        sage: As = [1,2,2]; Qs = [2,2,1,2]; m = 4
        sage: afsr(S, T, As, Qs, m, p, N)
        [1, 2, 2, 1, 0, 1, 0, 1, 2, 0, 0, 0, 3, 0, 2, 0, 2, 4, 0, 0]
        sage: Qp(5)(-487/283)
         1 + 2*5 + 2*5^2 + 5^3 + 5^5 + 5^7 + 2*5^8 + 3*5^12 + 2*5^14 + 2*5^16 + 4*5^17 + O(5^20)
        sage: alpha.padded_list(20)
        [1, 2, 2, 1, 0, 1, 0, 1, 2, 0, 0, 0, 3, 0, 2, 0, 2, 4, 0, 0]
        sage: multiplicative_order(mod(5,283)) # the period
        282


    REFERENCE:
        [KX] Klapper and Xu, "Algebraic feedback shift registers",
             Theoretical Com. Sci., vol 226, 1999,
             http://www.cs.uky.edu/~klapper/papers.html

    """
    r = len(As)
    # test if S and T will work
    if len(S) <> p:
        raise NameError("arg S is not a complete set of reps mod ", p)
    if len(T) <> p:
        raise NameError("arg T is not a complete set of reps mod ", p)
    Zp = IntegerModRing(p)
    S0 = [Zp(s) for s in S]
    S0.sort()
    if S0 <> [Zp(x) for x in range(p)]:
        raise NameError("arg S is not a complete set of reps mod ", p)
    T0 = [Zp(t) for t in T]
    T0.sort()
    if T0 <> [Zp(x) for x in range(p)]:
        raise NameError("arg T is not a complete set of reps mod ", p)
    # test if As are in S
    for x in As:
        if not(x in S):
            raise NameError("arg As is not a set of fills for ", S)
    # test if Qs are in T
    for x in Qs:
        if not(x in T):
            raise NameError("arg Qs is not a set of taps for ", T)
    if len(Qs) <> r+1:
        raise NameError("arg Qs is not a set of taps for ", T)
    # force m to be an initial memory value
    m = ZZ(m)
    # force p to be a prime
    p = next_prime(p-1)
    # force N to be an integer
    N = ZZ(N)
    conn_elmt = -Qs[0]+sum([Qs[i]*p**i for i in range(1,r+1)])
    newAs = As
    new_m = m
    a = S[0]
    if N <= r:
        return As
    for j in range(N-r):
        tau = sum([Qs[i]*newAs[-i] for i in range(1,r+1)]) + new_m
        for s in S:
            if Zp(tau) == Zp(s*Qs[0]):
                newAs.append(s)
                new_m = ZZ((tau-s*Qs[0])/p)
                #print j, tau, s, new_m, r
    return newAs


def afsr2(u, q, p, N): 
    """
    Input consists of 
        u, q - positive integers (q<p^(20))
        p - positive integer
        N - positive integer > r

    This function returns a triple:
     -  u,q in least denominator form,
     -  a list of length N generated by the AFSR associated 
        to (u, q, p), as in Theorem 10, [KX])
  
    EXAMPLES:
        sage: u = 972; q = 1001; p = 5; N = 20
        sage: adic_seq(u, q, p, N)
        (972, 1001, [2, 4, 3, 1, 1, 4, 0, 4, 0, 3, 1, 3, 4, 1, 3, 2, 4, 2, 3, 3])
        sage: afsr2(u, q, p, N)
        (972, 1001, [2, 4, 3, 1, 1, 4, 0, 4, 0, 3, 1, 3, 4, 1, 3, 2, 4, 2, 3, 3])

    REFERENCE:
        [KX] Klapper and Xu, "Algebraic feedback shift registers",
             Theoretical Com. Sci., vol 226, 1999,
             http://www.cs.uky.edu/~klapper/papers.html

    """
    S = range(p)
    U = range(q)
    G = gcd(u,q)
    u = ZZ(u/G)
    q = ZZ(q/G)
    gamma = 1
    for i in U:
        if i*p%q==1:
            gamma = i
    rho = 1
    for i in S:
        if i*q%p==1:
            rho = i
    aseq=[];useq=[];j=0;
    useq+=[(u*gamma**k)%q-q for k in range(N)]
    while (u%(p**(j+1))==0 and  u>=0):
        aseq+=[0]
        j+=1
    aseq+=[rho*((u*gamma**j)%q)%p]
    aseq+=[rho*useq[k]%p for k in range(j+1,N)]
    return u,q,aseq

#######################################

def index(x,M):
    """
    Inputs a non-zero integer and outputs the 
    integer part of its log.

    """
    return int(log(abs(x))/log(M))

def height(x,y,M):
    """
    Inputs a pair of non-zero integers and outputs the 
    max of the integer part of their log base M.

    EXAMPLE:
        sage: height(32,33,2)
        5
    """
    return int(max(log(abs(x))/log(M),log(abs(y))/log(M)))

def interpolation_set(M):
    """
    Returns the interpolation set [1] associated to the 
    modulus M.

    Needed in rational_synthesis_xu.

    REFERENCE:
        [KX] Klapper and Xu, "Register Synthesis for Algebraic Feedback 
             Shift Registers Based on Non-Primes," (section 3)
             http://www.cs.uky.edu/~klapper/papers.html

    """
    if M==2:
        L = 5
    if M==3:
        L = 15
    if M>3:
        L = int(M**3/2)
    return range(-L,L+1)

def findst(h1,h2,i,M):
    """
    Solves h1*s+h2*t = 0 mod M^5.

    Needed in rational_synthesis_xu.

    """
    P = interpolation_set(M)
    for s in P:
        for t in P:
            if (s*h1+t*h2)%M**(i+5) == 0:
                return (s,t)

def maxpower(u,q,M):
    """
    Returns max power of M dividing u,q

    Needed in rational_synthesis_xu.

    EXAMPLES:
        sage: maxpower(1000, 200, 10)
        2
        sage: maxpower(1200,9,9)
        0


    """
    r = gcd(u,q)
    m = int(log(r,M))
    for t in range(2*m):
        if r%M**t != 0:
            return t-1
    return 0

def ifsinP(i,h,r,a,M,P):
    """"
    Returns s if s is in P and is non-zero.
    Returns 0 otherwise.

    Needed in rational_synthesis_xu.

    """
    B = 5
    ans = 0
    Pn0 = [x for x in P if x!=0]
    for s in Pn0:
        if s*(h[i]-a*r[i])%M**(i+B)==0:
            return s
    return 0

def adic_seq(u,q,x,N):
    """
    Returns N elements of the x-adic sequence for u/q,
    and u,q in least denominator form.

    Raises a ValueError if the denominator is divisible
    by the modulus x.

    EXAMPLES:
        sage: adic_seq(972,1001,10,15)
        (972, 1001, [2, 7, 9, 8, 2, 0, 1, 7, 9, 8, 2, 0, 1, 7, 9])
        sage: adic_seq(972,1001,10,20)
        (972, 1001, [2, 7, 9, 8, 2, 0, 1, 7, 9, 8, 2, 0, 1, 7, 9, 8, 2, 0, 1, 7])
    
    """
    b = []
    b_seq = []
    G = gcd(u,q)
    u = ZZ(u/G)
    q = ZZ(q/G)
    if gcd(q,x) != 1: #check if the denominator is relatively prime to x
        raise ValueError, "After reduction of u and q, x and q are not relatively prime" 
    x = ZZ(x)
    N = int(N)
    for i in range(N):
        if i == 0:
            S = range(x)
        else:
            S = range(int(b_seq[i-1]),int(x^(i+1)),int(x^i))
        for j in S: 
            if (u-q*j)%x^(i+1) == 0:
                if i == 0:
                    b.append(ZZ(j))
                else:
                    b.append(ZZ((j-b_seq[i-1])/x^(i)))
                b_seq.append(ZZ(j))
    return u,q,b


def rational_synthesis_xu(As, M, verbose = False):
    """
    This outputs the rational which "best fits" the
    decimal expansion in a.

    INPUT:
         As - a list of k integers in [0,M-1]
         M - a modulus (such as 10)

    EXAMPLES:
        sage: adic_seq(972,1001,10,20)
        (972, 1001, [2, 7, 9, 8, 2, 0, 1, 7, 9, 8, 2, 0, 1, 7, 9, 8, 2, 0, 1, 7])
        sage: rational_synthesis_xu([2, 7, 9, 8, 2, 0, 1, 7, 9, 8, 2], 10)
        (72849565132, 384693842281)
        sage: rational_synthesis_xu([2, 7, 9, 8, 2, 0, 1, 7, 9, 8, 2, 0, 1, 7, 9], 10)
        (10846647308, 62610590489)
        sage: rational_synthesis_xu([2, 7, 9, 8, 2, 0, 1, 7, 9, 8, 2, 0, 1, 7, 9, 8, 2, 0, 1, 7], 10)
        (972, 1001)
        sage: rational_synthesis_xu([0,1,0,0,0,0,0],10)
        (10, 1)
        sage: rational_synthesis_xu([0,0,1,0,0,0,0],10)
        (100, 1)
        sage: rational_synthesis_xu([1,0,0,0],10)
        (1, 1)
        sage: rational_synthesis_xu([1,0,0],10)
        (-9, 91)
        sage: rational_synthesis_xu([2, 7, 8, 2, 1, 7, 8, 2, 1, 7, 8, 2, 1],10)
        (72, 101)
        sage: rational_synthesis_xu([2, 7, 8, 2, 1, 7, 8, 2, 1],10)
        (149711144, 196122577)

    QUESTION: Do you need 3 times the period to arrive with the
    correct answer?

    """
    k = len(As)
    B = 3
    P = interpolation_set(M)
    a = 1+M*sum([As[j]*M**j for j in range(k)])
    m = 0
    h = [0 for i in range(k+1)]
    r = [0 for i in range(k+1)]
    h[0] = 0; r[0] = 1
    h[1] = 1+M*sum([As[j]*M**j for j in range(B)])
    r[1] = 1+M**B
    for i in range(1, k):
        if verbose:
            print "start", numerator((-1+(h[i]/r[i]))/M), denominator((-1+(h[i]/r[i]))/M)
        if (h[i]-a*r[i])%(M**(i+1)) <> 0:
            s = ifsinP(i,h,r,a,M,P)
            t = 0
            if s<>0:
                h[i+1] = s*h[i]
                r[i+1] = s*r[i]
                if verbose:
                    print "type 1 update", i, numerator((-1+(h[i]/r[i]))/M), denominator((-1+(h[i]/r[i]))/M)
            else:
                s,t = findst(h[i]-r[i]*a,M**(i-m)*(h[m]-r[m]*a),i,M)
                h[i+1] = s*h[i]+t*M**(i-m)*h[m]
                r[i+1] = s*r[i]+t*M**(i-m)*r[m]
                if verbose:
                    print "type 2 update", i, numerator((-1+(h[i]/r[i]))/M), denominator((-1+(h[i]/r[i]))/M)
            cond1 = bool(height(h[i+1],r[i+1],M)>height(h[i],r[i],M))
            cond2 = bool(height(h[i],r[i],M) <= i-m+height(h[m],r[m],M))
            if (cond1 and cond2 and t<> 0):
                m = i
                if verbose:
                    print "turn updating", i, numerator((-1+(h[i]/r[i]))/M), denominator((-1+(h[i]/r[i]))/M)
        else:
              h[i+1] = h[i]
              r[i+1] = r[i]
    u = numerator((-1+(h[k]/r[k]))/M)
    q = denominator((-1+(h[k]/r[k]))/M)
    t = maxpower(u,q,M)
    return u/M**t,q/M**t

###################################################



def walsh_transform_1d(s,k):
    """
    This computes the Walsh(-Fourier) transform of the 
    sequence s of length n.
    
    INPUT:
        s - a sequence of n elements in GF(2), represented as a list
        k - an integer in the interval [0,n-1].

    EXAMPLES:
        sage: F = GF(2)
        sage: V = F^4
        sage: s = [F.random_element() for i in range(100)]
        sage: max([walsh_transform_1d(s,i) for i in range(100)])
        10
        sage: min([walsh_transform_1d(s,i) for i in range(100)])
        -14
        sage: s = [F.random_element() for i in range(10000)]
        sage: max([walsh_transform_1d(s,i) for i in range(100)])
        124
        sage: min([walsh_transform_1d(s,i) for i in range(100)])
        -16

    """
    n = len(s)
    terms = [(-1)^(s[i]+i*k) for i in range(n)]
    return sum(terms)

def walsh_transform(f, a):
    """
    This computes the Walsh(-Fourier) transform of the 
    Boolean function f on $GF(2)^n$ at $a \in GF(2)^n$.
    
    INPUT:
        s - a sequence of n elements in GF(2), represented as a list
        a - an element of GF(2)^n.

    EXAMPLES:
        sage: F = GF(2)
        sage: V = F^4
        sage: f = bent_function_standard
        sage: a = V([1,0,1,1])
        sage: walsh_transform(f, a)
        -4
        sage: a = V([0,0,0,0])
        sage: walsh_transform(f, a)
        4
        sage: WT = [walsh_transform(f, a) for a in V]; WT
        [4, 4, 4, -4, 4, 4, 4, -4, 4, 4, 4, -4, -4, -4, -4, 4]
        sage: 2^(4-1) - (1/2)*max([abs(x) for x in WT])
        6

    WARNING: This seems to differ by a sign from the method
    walsh_hadamard_transform in sage.crypto.boolean_function.

    """
    F = a[0].parent()
    n = len(a)
    terms = [(-1)^(f(x)+x.dot_product(a)) for x in F**n]
    return sum(terms)

###################################################
####### bent Boolean functions
###################################################

def bent_function_standard(x):
    """
    This function returns the value of a bent function at x.
    The bent function is constructed using the "simplest" 
    construction [1].

    INPUT:
        n - an even integer >= 4
        x - an element of GF(2)^n

    EXAMPLES:
        sage: V = GF(2)^4
        sage: s = [bent_function_standard(x) for x in V]
        sage: s
        [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0]


    REFERENCES:
        [1] http://en.wikipedia.org/wiki/Bent_function

    """
    n = x.degree()
    if n<1 or n%2 <> 0:
        raise ValueError("%s must be a positive even integer"%n)
    m = int(n/2)
    terms = [x[2*i]*x[2*i+1] for i in range(m)]
    return sum(terms)
    
def bent_function_MM(x, M = 1, g = lambda x: GF(2)(0) ):
    """
    This function returns the value of a bent function at x.
    The bent function is constructed using the Maiorana-McFarland
    construction [1].

    INPUT:
        x - an element of GF(2)^n
        M - an invertible n/2 x n/2 matrix over GF(2)
        g - a boolean function on GF(2)^(n/2)

    EXAMPLES:
        sage: V = GF(2)^4
        sage: s = [bent_function_MM(x) for x in V]
        sage: s
        [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0]
        sage: b = lambda x: bent_function_MM(x)
        sage: [walsh_transform(b, a) for a in V]
        [4, 4, 4, 4, 4, -4, 4, -4, 4, 4, -4, -4, 4, -4, -4, 4]
        sage: V = GF(2)^8
        sage: M = GL(4, GF(2)).random_element(); M
        [0 0 1 1]
        [1 1 0 1]
        [1 0 1 0]
        [1 0 0 0]
        sage: g = lambda x: bent_function_standard(x) 
        sage: s = [bent_function_MM(x, M, g) for x in V]
        sage: len(s)
        256
        sage: b = lambda x: bent_function_MM(x, M, g)
        sage: W = [walsh_transform(b, a) for a in V]
        #list of \pm 16 of length 256

    REFERENCES:
        [1] C. Carlet, "Boolean Functions for Cryptography and 
            Error Correcting Codes," preprint (available on www).

    """
    n = x.degree()
    if n<1 or n%2 <> 0:
        raise ValueError("%s must be a positive even integer"%n)
    m = int(n/2)
    if M==1:
        M = identity_matrix(GF(2), m)
    V = GF(2)^m
    v1 = V([x[i] for i in range(m)])
    v2 = V([x[i] for i in range(m,n)])
    return v1*M*v2+g(v2)

def bent_function_rothaus(x, g = lambda x: GF(2)(0)):
    """
    This function returns the value of a bent function at x.
    The bent function is constructed using the Rothaus 
    construction [1]. (This construction was generalized by
    the Maiorana-McFarland construction.

    INPUT:
        x - an element of GF(2)^n
        g - a boolean function on GF(2)^(n/2)

    EXAMPLES:
        sage: V = GF(2)^4
        sage: s = [bent_function_rothaus(x) for x in V]
        sage: s
        [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0]
        sage: b = lambda x: bent_function_rothaus(x)
        sage: [walsh_transform(b, a) for a in V]
        [4, 4, 4, 4, 4, -4, 4, -4, 4, 4, -4, -4, 4, -4, -4, 4]

    REFERENCES:
        [1] O. S. Rothaus: On "Bent" Functions. 
            J. Comb. Theory, Ser. A 20 (1976) 300-305

    """
    n = x.degree()
    if n<1 or n%2 <> 0:
        raise ValueError("%s must be a positive even integer"%n)
    m = int(n/2)
    V = GF(2)^m
    v1 = V([x[i] for i in range(m)])
    v2 = V([x[i] for i in range(m,n)])
    M = identity_matrix(GF(2), m)
    return v1*M*v2+g(v1)

###################################################

def boolean_fcn_synthesis(f, n, p):
    """
    Returns rational number whose p-adic expansion matches the values of 
    f:GF(p)^n -> GF(p).

    INPUT:

    EXAMPLES:
        sage: f = lambda x: x[0]*x[1]+x[2]*x[3]
        sage: boolean_fcn_synthesis(f, 4, 2)
        -30856/65535
        sage: boolean_fcn_synthesis(f, 4, 3)
        -158298070546506595701961644853448561652/221713244121518884974124815309574946401
        sage: V = GF(2)^3
        sage: V.list()
        [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
        sage: f = lambda x: x[0]               
        sage: boolean_fcn_synthesis(f, 3, 2)
        -2/3
        sage: f = lambda x: x[2]
        sage: boolean_fcn_synthesis(f, 3, 2)
        -16/17
        sage: adic_seq(-16, 17, 2, 15)
        (-16, 17, [0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1])
    """
    F = GF(p)
    V = GF(p)**n
    return (1-p**(p**n))**(-1)*sum([int(f(V[j]))*p**j for j in range(p**n)])


###################################################
####### monotone Boolean functions
###################################################

def is_monotone(f, n = 0):
    """
    Returns True if x<y => f(x)<f(y), where f is a Boolean
    function of n variables. The function f can be entered either
    by assigning values, using lambda notation and specifying
    the number of variables n explicitl, or as a BooleanFunction
    instance, where n is determined from the function.

    EXAMPLES:
        sage: V = GF(2)^4
        sage: Vlist = V.list()
        sage: flist = [0,0,0,1,0,1,0,1,0,0,0,1,0,1,1,1]
        sage: f = lambda x: GF(2)(flist[Vlist.index(x)])
        sage: WT = [walsh_transform(f, a) for a in V]; WT
        [2, 10, 6, -2, 6, -2, 2, -6, 2, 2, -2, -2, -2, -2, 2, 2]
        sage: is_monotone(f, 4)
        True
        sage: V = GF(2)^8 
        sage: Vlist = (GF(2)^3).list()    
        sage: F = [v for v in V if is_monotone(lambda x: GF(2)(list(v)[Vlist.index(x)]), 3)==True]
        sage: len(F)
        20
        sage: [(v, H*vector(QQ, [(-1)^x for x in list(v)])) for v in F]
        [((0, 0, 0, 0, 0, 0, 0, 0), (8, 0, 0, 0, 0, 0, 0, 0)), 
         ((0, 0, 0, 0, 0, 0, 0, 1), (6, 2, 2, -2, 2, -2, -2, 2)),
         ((0, 0, 0, 1, 0, 0, 0, 1), (4, 4, 4, -4, 0, 0, 0, 0)),
         ((0, 0, 0, 0, 0, 1, 0, 1), (4, 4, 0, 0, 4, -4, 0, 0)),
         ((0, 0, 0, 1, 0, 1, 0, 1), (2, 6, 2, -2, 2, -2, 2, -2)),
         ((0, 1, 0, 1, 0, 1, 0, 1), (0, 8, 0, 0, 0, 0, 0, 0)),
         ((0, 0, 0, 0, 0, 0, 1, 1), (4, 0, 4, 0, 4, 0, -4, 0)),
         ((0, 0, 0, 1, 0, 0, 1, 1), (2, 2, 6, -2, 2, 2, -2, -2)),
         ((0, 0, 1, 1, 0, 0, 1, 1), (0, 0, 8, 0, 0, 0, 0, 0)),
         ((0, 0, 0, 0, 0, 1, 1, 1), (2, 2, 2, 2, 6, -2, -2, -2)),
         ((0, 0, 0, 1, 0, 1, 1, 1), (0, 4, 4, 0, 4, 0, 0, -4)),
         ((0, 1, 0, 1, 0, 1, 1, 1), (-2, 6, 2, 2, 2, 2, -2, -2)),
         ((0, 0, 1, 1, 0, 1, 1, 1), (-2, 2, 6, 2, 2, -2, 2, -2)),
         ((0, 1, 1, 1, 0, 1, 1, 1), (-4, 4, 4, 4, 0, 0, 0, 0)),
         ((0, 0, 0, 0, 1, 1, 1, 1), (0, 0, 0, 0, 8, 0, 0, 0)),
         ((0, 0, 0, 1, 1, 1, 1, 1), (-2, 2, 2, -2, 6, 2, 2, -2)),
         ((0, 1, 0, 1, 1, 1, 1, 1), (-4, 4, 0, 0, 4, 4, 0, 0)),
         ((0, 0, 1, 1, 1, 1, 1, 1), (-4, 0, 4, 0, 4, 0, 4, 0)),
         ((0, 1, 1, 1, 1, 1, 1, 1), (-6, 2, 2, 2, 2, 2, 2, 2)),
         ((1, 1, 1, 1, 1, 1, 1, 1), (-8, 0, 0, 0, 0, 0, 0, 0))]
        sage: B = BooleanPolynomialRing(5,'x')
        sage: x0,x1,x2,x3,x4 = B.gens()
        sage: f = BooleanFunction(x0 + x1)
        sage: f
        Boolean function with 5 variables
        sage: is_monotone(f)
        False
        sage: g=BooleanFunction(x0)
        sage: is_monotone(g)
        True
        sage: h=BooleanFunction(x0+x1+x0*x1+x2+x3+x2*x3)
        sage: is_monotone(h)
        False
        sage: h=BooleanFunction(x0+x1+x2+x0*x1+x0*x2+x1*x2+x0*x1*x2)
        sage: is_monotone(h)
        True
        sage: f=BooleanFunction(x0*x2+x1+x0*x2*x1)
        sage: is_monotone(f)
        True

    """
    b2i = GFn_to_integer
    try:
        V = GF(2)**n
        for x in V:
            sx = x.nonzero_positions()
            for y in V:
                sy = y.nonzero_positions()
                if Set(sx).issubset(Set(sy)):
                    if not(f(x)<=f(y)):
                        return False, x, y
        return True
    except:
            n = f.nvariables()
            V = GF(2)^n
            Vlist = V.list()
            t = True
            i = 0
            while (t and i<2^n):
                if f(i):
                    J = [b2i(v) for v in V if Set(v.support()).issuperset(Set((Vlist[i]).support()))]
                    for j in J:
                        if f(j)==0:
                            t = False
                i+=1
            return t

def boolean_monomial(u):
    """
    Returns the monomial x^u.

    INPUT:
        u - element of GF(2)^n

    OUTPUT:
        f - Boolean monomial

    EXAMPLES:
        sage: F=GF(2)
        sage: V=F^4
        sage: u=V([1,0,1,1])
        sage: boolean
        boolean_fcn_synthesis  boolean_monomial       
        sage: boolean_monomial(u)
        Defining x0, x1, x2, x3
        x0*x2*x3
        sage: g0=boolean_monomial(u)
        Defining x0, x1, x2, x3
        sage: g0
        x0*x2*x3
    """
    from sage.crypto.boolean_function import BooleanPolynomial
    n = len(u)
    F = GF(2)
    V = F^n
    B = BooleanPolynomialRing(n,'x')
    B.inject_variables()
    f = B(0); f=f+1;
    for i in range(n):
      if u[i]==1:
        f=f*(B.gens())[i]
    return f

def monotone_from_support(L):
    """
    Returns a monotone Boolean function f whose set of least
    vectors in in the list L of vectors in GF(2)^n

    INPUT:
        L - list of elements of GF(2)^n

    OUTPUT:
        f - monotone Boolean function on n variables

    EXAMPLES:
        sage: V = GF(2)^4  
        sage: f = monotone_from_support([V([1,0,0,0]),V([0,0,0,1])])              
        sage: f.truth_table()                                 
        (False, True, False, True, False, True, False, True, True, True, 
         True, True, True, True, True, True)

    """
    from sage.crypto.boolean_function import BooleanPolynomial,BooleanFunction
    F = GF(2)
    n = L[0].degree()
    V = F^n
    B = BooleanPolynomialRing(n,'x')
    B.inject_variables()
    f = B(0); f = f+1;
    for u in L:
      f=f*(boolean_monomial(u)+1)
    return BooleanFunction(f+1)

def print_truth_table(f):

    """
    INPUT:
        f - Boolean function 

    OUTPUT:
        Prints the truth table of f.

    This provides a way to find how a Boolen function can be written
    linear combination of the atomic monotone functions in a nice way.

    EXAMPLES:
        sage: V = GF(2)^4  
        sage: f = monotone_from_support([V([1,0,0,0]),V([0,0,0,1])])              
        sage: f.truth_table()                                 
        (False, True, False, True, False, True, False, True, True, True, 
         True, True, True, True, True, True)
        sage: print_truth_table(f)
        (0, 0, 0, 0) 0 False
        (1, 0, 0, 0) 1 True
        (0, 1, 0, 0) 0 False
        (1, 1, 0, 0) 0 True
        (0, 0, 1, 0) 0 False
        (1, 0, 1, 0) 0 True
        (0, 1, 1, 0) 0 False
        (1, 1, 1, 0) 0 True
        (0, 0, 0, 1) 1 True
        (1, 0, 0, 1) 1 True
        (0, 1, 0, 1) 0 True
        (1, 1, 0, 1) 0 True
        (0, 0, 1, 1) 0 True
        (1, 0, 1, 1) 0 True
        (0, 1, 1, 1) 0 True
        (1, 1, 1, 1) 0 True
    """ 
    b2i = GFn_to_integer
    n = f.nvariables()
    V = GF(2)^(n)
    W = GF(2)^(2^n)
    A = matrix([W(monotone_from_support([x]).truth_table()) for x in V])
    fvec = W(f.truth_table())
    c_v = fvec*A^(-1)
    for x in V:
        print str(x) + " " + str(c_v[b2i(x)]) + " " + str(f(b2i(x)))


def supp_compare(u,v):
    """
    INPUT:
        u,v - vectors in GF(2)^n

    OUTPUT:
        Returns True is u<=v and False otherwise.

    This function tests whether the support of u is a subset of the support
    of v.
    
    EXAMPLES:

    """ 
    supp_u = Set(u.nonzero_positions())
    supp_v = Set(v.nonzero_positions())

    return supp_u.issubset(supp_v)

###################################################
####### Cayley graph for Boolean functions
###################################################

"""
Some routines for computing with the Cayley graph associated with a
Boolean function.


Let f be a Boolean function on GF(2)^n. The Cayley graph 
of f is defined to be the graph

             Gamma_f = (GF(2)^n, E_f ), 

whose vertex set is GF(2)^n and the set of edges is defined by 

            E_f ={(u,v) in GF(2)^n | f(u+v)=1}.

The adjacency matrix A_f is the matrix whose entries are 

        A_{i,j} = f(b(i) + b(j)), 

where b(k) is the binary representation of the integer k.
Note Gamma_f is a regular graph of degree wt(f), where wt denotes
the Hamming weight of f when regarded as a vector of values
(of length 2^n).

Recall that, given a graph Gamma and its adjacency matrix A, the spectrum 
Spec(Gamma) is the multi-set of eigenvalues of A. It turns out that the
spectrum of Gamma_f is equal to the Walsh-Hadamard transform of
f when regarded as a vector of (integer) 0,1-values (of length 2^n).
(This nice fact seems to have first appeared in [2], [3].)

A graph is regular of degree r (or r-regular) if every vertex 
has degree r (number of edges incident to it).	We say that
an r-regular graph Gamma is a strongly regular graph (srg) 
with parameters (v, r, d, e)  (for nonnegative integers e, d) 
provided, for all vertices u, v the number of vertices adjacent 
to both u, v is equal to
              
          e, if u, v are adjacent,
          d, if u, v are nonadjacent.

It turns out tht f is bent iff Gamma_f is strongly regular and
e = d (see [3] and [4]).


Theorem (Stanica) 
Let f be a Boolean function, and let G_f be 
the associated Cayley graph with g being the multiplicity of its 
lowest eigenvalue lambda_v. Assume that the eigenvalues of G_f are ordered as 
lambda_1 \leq lambda_2\leq \dots \leq lambda_v.	Then, 
min{ g + 1, 1 - lambda_v/lambda_2} \leq chi(Gf ) \leq |Omega_f | 
(provided lambda_2 != 0).


The following Sage computations illustrate these and other theorems 
in [1], [2], [3], [4].

EXAMPLES:
    sage: V = GF(2)^4
    sage: f = lambda x: x[0]*x[1]+x[2]*x[3]
    sage: Gamma = boolean_cayley_graph(f, 4)
    sage: Gamma.spectrum()
    [6, 2, 2, 2, 2, 2, 2, -2, -2, -2, -2, -2, -2, -2, -2, -2]
    sage: [walsh_transform(f, a) for a in V]
    [4, 4, 4, -4, 4, 4, 4, -4, 4, 4, 4, -4, -4, -4, -4, 4]
    sage: Omega_f = [v for v in V if f(v)==1] 
    sage: len(Omega_f)
    6
    sage: Gamma.is_bipartite()
    False
    sage: Gamma.is_hamiltonian()
    True
    sage: Gamma.is_planar()     
    False
    sage: Gamma.is_regular()
    True
    sage: Gamma.is_eulerian()
    True
    sage: Gamma.is_connected()
    True
    sage: Gamma.is_triangle_free()
    False
    sage: Gamma.diameter()
    2
    sage: Gamma.degree_sequence()
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]
    sage: H = matrix(QQ, 16, 16, [(-1)^(Vlist[x[0]]).dot_product(Vlist[x[1]]) for x in C])
    sage: flist = vector(QQ, [int(f(v)) for v in V])
    sage: H*flist  
    (6, -2, -2, 2, -2, -2, -2, 2, -2, -2, -2, 2, 2, 2, 2, -2)
    sage: A = matrix(QQ, 16, 16, [f(Vlist[x[0]]+Vlist[x[1]]) for x in C])
    sage: A.eigenvalues()
    [6, 2, 2, 2, 2, 2, 2, -2, -2, -2, -2, -2, -2, -2, -2, -2]



##### another example


    sage: V = GF(2)^3
    sage: f = lambda x: x[0]*x[1]+x[2]
    sage: Omega_f = [v for v in V if f(v)==1] 
    sage: len(Omega_f)
    4
    sage: Gamma = boolean_cayley_graph(f, 4)
    sage: # not strongly regular
    sage: Gamma.spectrum()
    [4, 2, 0, 0, 0, -2, -2, -2]
    sage: 
    sage: Gamma.is_bipartite()
    False
    sage: Gamma.is_hamiltonian()
    True
    sage: Gamma.is_planar()     
    False
    sage: Gamma.is_regular()
    True
    sage: Gamma.is_eulerian()
    True
    sage: Gamma.is_connected()
    True
    sage: Gamma.is_triangle_free()
    False
    sage: Gamma.diameter()
    2
    sage: Gamma.degree_sequence()
    [4, 4, 4, 4, 4, 4, 4, 4]
    sage: H = matrix(QQ, 8, 8, [(-1)^(Vlist[x[0]]).dot_product(Vlist[x[1]]) for x in C])
    sage: A = matrix(QQ, 8, 8, [f(Vlist[x[0]]+Vlist[x[1]]) for x in C])
    sage: A.eigenvalues()
    [4, 2, 0, 0, 0, -2, -2, -2]
    sage: V = GF(2)^3
    sage: C = CartesianProduct(range(8), range(8))
    sage: f = lambda x: x[0]*x[1]+x[2]
    sage: # not bent
    sage: Omega_f = [f(x) for x in Vlist]
    sage: Omega_f
    [0, 0, 0, 1, 1, 1, 1, 0]
    sage: Gamma = boolean_cayley_graph(f, 3)
    sage: Gamma = Graph(E)
    sage: Gamma.vertices()
    [0, 1, 2, 3, 4, 5, 6, 7]
    sage: Gamma.spectrum()
    [4, 2, 0, 0, 0, -2, -2, -2]
    sage: Gamma.coloring()
    [[5, 2], [6, 7], [0, 1], [3, 4]]


REFERENCES:
 [1] Pantelimon Stanica, Graph eigenvalues and Walsh spectrum 
     of Boolean functions, INTEGERS 7(2) (2007), #A32.

 [2] Anna Bernasconi, Mathematical techniques for the analysis 
     of Boolean functions, Ph. D. dissertation TD-2/98, 
     Universit di Pisa-Udine, 1998.

 [3] Anna Bernasconi and Bruno Codenotti, Spectral Analysis 
     of Boolean Functions as a Graph Eigenvalue Problem,
     IEEE TRANSACTIONS ON COMPUTERS, VOL. 48, NO. 3, MARCH 1999.

 [4] A. Bernasconi, B. Codenotti, J.M. VanderKam. A Characterization 
     of Bent Functions in terms of Strongly Regular Graphs, 
     IEEE Transactions on Computers, 50:9 (2001), 984-985. 

 [5] Michel Mitton, Minimal polynomial of Cayley graph
     adjacency matrix for Boolean functions, preprint, 2007.

 [6] ------, On the Walsh-Fourier analysis of Boolean
     functions, preprint, 2006.

"""

def boolean_cayley_graph(f, n):
    """

    EXAMPLES:
        sage: V = GF(2)^4
        sage: Vlist = V.list()
        sage: flist = [0,0,0,1,0,1,0,1,0,0,0,1,0,1,1,1]
        sage: f = lambda x: GF(2)(flist[Vlist.index(x)])
        sage: boolean_cayley_graph(f, 4)
        Graph on 16 vertices

    """
    V = GF(2)**n
    C = CartesianProduct(range(2**n), range(2**n))
    Vlist = V.list()                  
    E = [(x[0],x[1]) for x in C if f(Vlist[x[0]]+Vlist[x[1]])==1]
    E = Set([Set(s) for s in E])
    E = [tuple(s) for s in E] 
    Gamma = Graph(E)
    return Gamma
