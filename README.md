## The structure of shell and basis function for gtf in fch file

### `Shell type`:
- -4: 9g
- -3: 7f
- -2: 5d
- -1: sp
- 0: 1s
- 1: 3p
- 2: 6d
- 3: 10f
- 4: 15g

```
Shell types                                I   N=          13
           0          -1          -1           2           0           0
           1           0           0           1           0           0
           1
```

### `Number of primitives per shell`:

The contraction number of each shell, the number of elements of primitives per shell equals to the number of shell types.

```
Number of primitives per shell             I   N=          13
           6           3           1           1           3           1
           1           3           1           1           3           1
           1
```

### `Primitive exponents` and `Contraction coefficients`:

They are the same of every compenents for the same shell type, i.e. the same for the X Y Z components of the shell type `1`, then its element number is equal to the sum of each contraction number of the corresponding shell.

```
Primitive exponents                        R   N=          26
  3.04752488E+03  4.57369518E+02  1.03948685E+02  2.92101553E+01  9.28666296E+00
  3.16392696E+00  7.86827235E+00  1.88128854E+00  5.44249258E-01  1.68714478E-01
  8.00000000E-01  1.87311370E+01  2.82539436E+00  6.40121692E-01  1.61277759E-01
  1.10000000E+00  1.87311370E+01  2.82539436E+00  6.40121692E-01  1.61277759E-01
  1.10000000E+00  1.87311370E+01  2.82539436E+00  6.40121692E-01  1.61277759E-01
  1.10000000E+00
Contraction coefficients                   R   N=          26
  1.83473713E-03  1.40373228E-02  6.88426223E-02  2.32184443E-01  4.67941348E-01
  3.62311985E-01 -1.19332420E-01 -1.60854152E-01  1.14345644E+00  1.00000000E+00
  1.00000000E+00  3.34946043E-02  2.34726953E-01  8.13757326E-01  1.00000000E+00
  1.00000000E+00  3.34946043E-02  2.34726953E-01  8.13757326E-01  1.00000000E+00
  1.00000000E+00  3.34946043E-02  2.34726953E-01  8.13757326E-01  1.00000000E+00
  1.00000000E+00
```
For SP shell(`-1` in Shell types), we only take the coefficients of the last three in `P(S=P) Contraction coefficients`, i.e. the P type Shell instead using the coefficients in `Contraction coefficients`.

```fortran
if ( shl_type(ishl) == -1 .and. ibsshl /= 1 ) then
    ! sp type shell, and the last three GTFs are p-functions
    basis%gtfcoeff(ictr) = SP_ctr_cff(ictr+iprim-1) 
else
    ! sp type shell, the first is s-function
    basis%gtfcoeff(ictr) = ctr_cff(ictr+iprim-1)
```


```
P(S=P) Contraction coefficients            R   N=          26
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  6.89990666E-02  3.16423961E-01  7.44308291E-01  1.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00  0.00000000E+00
  0.00000000E+00
```

### Orbit coefficients:

We only concentrate on the number of the contracted basis shell, i.e. the sum of the components of corresponding shell types is the number of basis function used in the linear combination of orbits.

```
Alpha Orbital Energies                     R   N=          30
 -1.01884160E+01 -6.84051688E-01 -4.17833764E-01 -4.17833764E-01 -2.24551475E-01
  1.19389709E-01  1.85516107E-01  1.85516107E-01  5.35951550E-01  5.35951550E-01
  5.37185565E-01  6.81886841E-01  8.97771005E-01  8.97771005E-01  9.52649005E-01
  1.43519157E+00  1.43519157E+00  1.92622445E+00  1.92622445E+00  1.97602171E+00
  2.00907442E+00  2.28717507E+00  2.57565312E+00  2.57565312E+00  2.74315224E+00
  2.74315224E+00  3.28192103E+00  3.48818760E+00  3.48818760E+00  4.27195390E+00
```

```
Alpha MO coefficients                      R   N=         900
```

Which means the number of basis function is 30 = 1 + 4 + 4 + 6 + 1 + 1 + 3 + 1 + 1 + 3 + 1 + 1 + 3.