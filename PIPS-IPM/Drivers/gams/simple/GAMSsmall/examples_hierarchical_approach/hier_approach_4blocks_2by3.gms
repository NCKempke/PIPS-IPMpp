* ./gamsexample.sh -NP=4 -BLOCKS=5 -GAMSFILE=./hier_approach_4blocks_2by3

Set i rows    / i1*i15 /
    j columns / j1*j19 /;

parameter g(j) obj coefficients / j1 1, j2 1, j3 1, j4 1, j5 1, j6 1, j7 1, j8 1, j9 1, j10 1, j11 1, j12 2, j13 1, j14 1 , j15 1, j16 1, j17 1, j18 1, j19 1 /
          bA(i) right hand side       / i1 6, i2 3, i3 4, i4 1, i5 4, i6 1, i7 4, i8 1, i9 4, i10 1, i11 10, i12 6, i13 5, i14 8, i15 17 /
*          clow(i) c left hand side   / i1 9, i2 3, i3 5, i4 -0.5, i5 5, i6 -0.5, i7 5, i8 -0.5, i9 13.5, i10 3, i11 11, i12 18 /
          cupp(i) c right hand side   / i1 13, i2 7, i3 9, i4 2.5, i5 9, i6 2.5, i7 9, i8 2.5, i9 8.5, i10 2.3, i11 20.5, i12 17, i13 12.1, i14 16, i15 34.1 /


* in this example the singletonColumnPresolver should substitute j1 (free) and put it into the objective
Table A(i,j)
    j1    j2    j3    j4    j5    j6    j7    j8    j9   j10   j11   j12   j13   j14   j15   j16   j17   j18   j19
i1   1     1     1           1     1     1
i2         1           1    -1     1     1
i3   1    -1                 1                 1     1     1 
i4         1                -1                -1     1     1 
i5   1    -1                       1                             1     1     1
i6         1                      -1                            -1     1     1
i7   1    -1                 1                                                     1     1     1
i8         1                -1                                                    -1     1     1
i9   1    -1                       1                                                                 1     1     1
i10        1                      -1                                                                -1     1     1
i11  1     1     1           1     1                 1     2     1     1       
i12 -1                 1    -1     1                                   1     2     1     2 
i13 -1    -1     1     1                                                                 1     2     1     1
i14  1    -1                 1    -1           1     1           1     1           1     1           1     1
i15  1     1                 1     1     1           2     1           2     1           2     1           2     1
;                         
*    1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
* expected values for x full determined by Ax=b
* objective = 14
Table C(i,j)
    j1    j2    j3    j4    j5    j6    j7    j8    j9   j10   j11   j12   j13   j14   j15   j16   j17   j18   j19
i1   2     2     2           2     2     2
i2         2           2    -2     2     2
i3   2    -2                 2                 2     2     2 
i4         2                -2                -2     2     2 
i5   2    -2                       2                             2     2     2
i6         2                      -2                            -2     2     2
i7   2    -2                 2                                                     2     2     2
i8         2                -2                                                    -2     2     2
i9   2    -2                       2                                                                 2     2     2
i10        2                      -2                                                                -2     2     2
i11  2     2     2           2     2                 2     4     2     2       
i12 -2                 2    -2     2                                   2     4     2     4
i13 -2    -2     2     2                                                                 2     4     2     2
i14  2    -2                 2    -2           2     2           2     2           2     2           2     2
i15  2     2                 2     2     2           4     2           4     2           4     2           4     2;

Variables          x(j) * / j4.lo -5, j4.up 5, j7.lo 1, j10.lo 1, j13.lo 1, j16.lo 1, j19.lo 1 /
Variable           z      objective variable
Equations          e(i)   equality equations
*                   ge(i)  greater than inequality equations
                   le(i)  less than inequality equations
                   defobj objective function;

defobj.. z =e= sum(j, g(j)*x(j));
e(i)..   sum(j, A(i,j)*x(j)) =e= bA(i);
* ge(i)..  sum(j, C(i,j)*x(j)) =g= clow(i);
le(i)..  sum(j, C(i,j)*x(j)) =l= cupp(i);

Model m /all/ ;


$ifthen %METHOD%==PIPS
*annotations for variables:
  x.stage('j8') = 2;
  x.stage('j9') = 2;
  x.stage('j10') = 2;
  x.stage('j11') = 3;
  x.stage('j12') = 3;
  x.stage('j13') = 3;
  x.stage('j14') = 4;
  x.stage('j15') = 4;
  x.stage('j16') = 4;
  x.stage('j17') = 5;
  x.stage('j18') = 5;
  x.stage('j19') = 5;
*annotations for equations:
  e.stage('i1') = 1;
  e.stage('i2') = 1;
  e.stage('i3') = 2;
  e.stage('i4') = 2;
  e.stage('i5') = 3;
  e.stage('i6') = 3;
  e.stage('i7') = 4;
  e.stage('i8') = 4;
  e.stage('i9') = 5;
  e.stage('i10') = 5;
  e.stage('i11') = 6;
  e.stage('i12') = 6;
  e.stage('i13') = 6;
  e.stage('i14') = 6;
  e.stage('i15') = 6;
*  ge.stage('i1') = 1;
*  ge.stage('i2') = 1;
*  ge.stage('i3') = 2;
*  ge.stage('i4') = 2;
*  ge.stage('i5') = 3;
*  ge.stage('i6') = 3;
*  ge.stage('i7') = 4;
*  ge.stage('i8') = 4;
*  ge.stage('i9') = 5;
*  ge.stage('i10') = 5;
*  ge.stage('i11') = 5;
*  ge.stage('i12') = 5;
  le.stage('i1') = 1;
  le.stage('i2') = 1;
  le.stage('i3') = 2;
  le.stage('i4') = 2;
  le.stage('i5') = 3;
  le.stage('i6') = 3;
  le.stage('i7') = 4;
  le.stage('i8') = 4;
  le.stage('i9') = 5;
  le.stage('i10') = 5;
  le.stage('i11') = 6;
  le.stage('i12') = 6;
  le.stage('i13') = 6;
  le.stage('i14') = 6;
  le.stage('i15') = 6;
  defobj.stage  = 6;


* For creation of gdx files:
$ echo jacobian hier_approach_4blocks_2by3.gdx > convertd.opt
  option lp=convertd;
  m.optfile = 1;
  solve m use lp min z;
$else
  option lp=cplex;
$ onecho > cplex.opt
*  lpmethod 4
  sopreind 0
$ offecho
  m.optfile = 1;
  solve m use lp min z;
$endif
