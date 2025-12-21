#set page(height: auto)

#let inv(n) = math.frac([1], n);

= Multiply

$
  (a+b i)(c + d i)=\
  a c + a d i + b c i + b d i^2=\
  a c + (a d + b c)i - b d=\
  a c - b d + (a d + b c)i
$

= Divide
$
  frac(a + b i, c + d i)=\
  frac((a + b i)(c - d i), (c + d i)(c + d i))=\
  frac(a c + b d + (b c - a d)i, c^2 + d^2)=\
  frac(a c + b d, c^2 + d^2) + frac(b c - a d, c^2 + d^2)i=\
$

= Square root
If $b$ is positive:
$
  sqrt(z)=sqrt(frac(a+sqrt(a^2+b^2), 2))+sqrt(frac(-a+sqrt(a^2+b^2), 2))i
$

If $b$ is negative:
$
  sqrt(z)=sqrt(frac(a+sqrt(a^2+b^2), 2))-sqrt(frac(-a+sqrt(a^2+b^2), 2))i
$

= Inverse
$
  inv(z) =\
  frac(a - b i, |z|^2)=\
  frac(a, |z|^2) - frac(b, |z|^2) i
$

= Natural Logarithm
$
  ln(z) =\
  ln(|z|e^(i theta))=\
  ln|z|+ln(e^(i theta))=\
  ln|z|+i theta
$

= Exponentiation

$
  z_1^(z_2)=\
  z_1^(c + d i)=\
  z_1^(c) z_1^(d i)=\
  z_1^(c) e^(ln(z_1)d i)=\
$

= Trig

$
  sin(z) & =sin(a)cosh(b) + cos(a)sinh(b)i \
  cos(z) & =cos(a)cosh(b) - sin(a)sinh(b)i \
  tan(z) & =frac(sin(z), cos(z)) \
  cot(z) & =frac(cos(z), sin(z)) \
  sec(z) & =inv(cos(z)) \
  csc(z) & =inv(sin(z)) \
$

= Inverse Trig

$
    arcsin(z) & =-i ln(sqrt(1-z^2) + i z) \
    arccos(z) & =i ln(1/i sqrt(1-z^2) + z) \
    arctan(z) & =arcsin(frac(z, sqrt(1+z^2))) \
  "arccot"(z) & =arctan(inv(z)) \
  "arcsec"(z) & =arccos(inv(z)) \
  "arccsc"(z) & =arcsin(inv(z)) \
$

= Hyperbolic Trig

$
  sinh(z) & =sinh(a)cos(b)+cosh(a)sin(b)i \
  cosh(z) & =cosh(a)cos(b)+sinh(a)sin(b)i \
  tanh(z) & =frac(sinh(z), cosh(z)) \
  coth(z) & =frac(cosh(z), sinh(z)) \
  sech(z) & =inv(cos(z)) \
  csch(z) & =inv(sin(z)) \
$

= Inverse Hyperbolic Trig

$
  "arcsinh"(z) & =ln(sqrt(z^2+1)+z) \
  "arccosh"(z) & =ln(sqrt(z^2-1)+z) \
  "arctanh"(z) & =1/2ln(frac(1+z, 1-z)) \
  "arccoth"(z) & ="arctanh"(inv(z))=1/2ln(frac(1-z, 1+z)) \
  "arcsech"(z) & ="arccosh"(inv(z)) \
  "arccsch"(z) & ="arcsinh"(inv(z)) \
$
