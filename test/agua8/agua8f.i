&INT0 
 TITLE = 'NO3',
 NATOM = 24, 
 NOPT=2,
 NORM = T, 
 ANG = T, 
 EX = T, 
 COUL = F, 
 SCF1 = T, 
 PROP = F, 
 FIELD = F, 
 SOL = F, 
 RESP1 = F/
  8        37.230000       71.550003       45.750000
  1       36.520000       71.430000       45.050003
  1       37.239998       70.760002       46.350002
  8        42.930000       68.610001       48.889999
  1       42.000000       68.540001       48.520004
  1       43.040001       69.500000       49.330002
  8        40.590000       70.209999       47.849998
  1       39.809998       69.899994       47.320000
  1       41.170002       70.790001       47.280003
  8        39.639999       66.770004       46.100002
  1       39.779999       65.930000       46.620003
  1       39.090000       67.410004       46.650002
  8        40.779999       70.959999       50.529999
  1       40.570000       70.520004       49.660000
  1       41.750000       71.220001       50.560001
  8        37.610001       68.959999       50.060001
  1       38.139999       68.470001       50.739998
  1       37.830002       69.940002       50.110001
  8        42.280003       71.540001       46.339996
  1       42.940002       72.279999       46.430000
  1       42.719997       70.750000       45.919998
  8        38.070000       68.440002       47.459999
  1       38.020000       68.430000       48.460003
  1       37.239998       68.029999       47.080002
&geo/
gaussian
 8  15   6
 6 2 1 4 1 1
 0 0 0 1 1 2
   5222.9022000   -0.001936
    782.5399400   -0.014851
    177.2674300   -0.073319
     49.5166880   -0.245116
     15.6664400   -0.480285
      5.1793599   -0.335943
     10.6014410    0.078806
      0.9423170   -0.567695
      0.2774746    1.000000
     33.4241260    0.017560
      7.6221714    0.107630
      2.2382093    0.323526
      0.6867300    0.483223
      0.1938135    1.000000
      0.8000000    1.000000
 8 10 10
 1 1 1 1 1 1 1 1 1 1
 0 0 0 0 0 0 0 1 1 2
    628.6475400    1.000000
    143.9976180    1.000000
     40.0859040    1.000000
     11.9849376    1.000000
      1.4560475    1.000000
      4.7140760    1.000000
      0.4059979    1.000000
      4.7140760    1.000000
      0.4059979    1.000000
      1.0000000    1.000000
gaussian
 1   6   3
 4 1 1
 0 0 0
     50.9991780    0.0096604761
      7.4832181    0.073728860
      1.7774676    0.29585808
      0.5193295    0.71590532
      0.1541100    1.000000
      0.7500000    1.000000
 1 4 4
 1 1 1 1  1 1 1
 0 0 0 0  1 1 1
     45.0000000    1.000000
      7.5000000    1.000000
      0.3000000    1.000000
      1.5000000    1.000000
endbasis
&SCFINP
 OPEN = F, 
 NMAX=300
 NCO = 40,
 NUNP = 0, 
 ATRHO = F, 
 VCINP = F, 
 DIRECT = T, 
 EXTR = F, 
 SHFT = F, 
 SHI =  1., 
 IDAMP = 0, 
 GOLD =  5., 
 TOLD =  1.E-06, 
 WRITE = F, 
 MEMO = T/
&EXCH 
 IEXCH=9
 INTEG = T, 
 DENS = T, 
 IGRID = 2, 
 IGRID2 = 1/


