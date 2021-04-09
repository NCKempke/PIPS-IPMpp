* ./gamsexample.sh -NP=3 -BLOCKS=4 -GAMSFILE=./hier_approach_2blocks_2by2

Set i rows    / i1*i12 /
    j columns / j1*j12 /;

parameter g(j) obj coefficients / j1 1, j2 1, j3 1, j4 1, j5 1, j6 1, j7 1, j8 1, j9 1, j10 1 /
          bA(i) right hand side  / i1 5, i2 2, i3 3, i4 0, i5 3, i6 0, i7 3, i8 0, i9 7, i10 2, i11 6, i12 10 /
*          clow(i) c left hand side   / i1 9, i2 3, i3 5, i4 -0.5, i5 5, i6 -05, i7 5, i8 -0.5, i9 13.5, i10 3, i11 11, i12 18 /
          cupp(i) c right hand side   / i1 11, i2 5, i3 7, i4 0.5, i5 7, i6 0.5, i7 7, i8 0.5, i9 14.5, i10 5, i11 13, i12 22 /


* in this example the singletonColumnPresolver should substitute j1 (free) and put it into the objective
Table A(i,j)
    j1    j2    j3    j4    j5    j6    j7    j8    j9   j10   j11   j12
i1   1     1     1           1     1
i2         1           1    -1     1
i3   1    -1                 1           1     1
i4         1                -1          -1     1
i5   1    -1                       1                 1     1
i6         1                      -1                -1     1
i7   1    -1                 1                                   1     1
i8         1                -1                                  -1     1                  
i9   1     1     1           1     1           1     1
i10 -1                 1    -1     1                       1     1
i11  1    -1                 1    -1     1     1     1     1     1     1
i12  1     1                 1     1           2           2           2
; 
*    1     1     1     1     1     1     1     1     1     1     1     1
* expected values for x full determined by Ax=b
* objective = 10 
Table C(i,j)
    j1    j2    j3    j4    j5    j6    j7    j8    j9   j10   j11   j12
i1   2     2     2           2     2
i2         2           2    -2     2
i3   2    -2                 2           2     2
i4         2                -2          -2     2
i5   2    -2                       2                 2     2
i6         2                      -2                -2     2
i7   2    -2                 2                                   2     2
i8         2                -2                                  -2     2                  
i9   2     2     2           2     2           2     2
i10 -2                 2    -2     2                       2     2
i11  2    -2                 2    -2     2     2     2     2     2     2
i12  2     2                 2     2           4           4           4
;

Variables          x(j) * / j4.lo -5, j4.up 5 /
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
  x.stage('j7') = 2;
  x.stage('j8') = 2;
  x.stage('j9') = 3;
  x.stage('j10') = 3;
  x.stage('j11') = 4;
  x.stage('j12') = 4;
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
  e.stage('i11') = 5;
  e.stage('i12') = 5;
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
  le.stage('i11') = 5;
  le.stage('i12') = 5;
  defobj.stage  = 5;


* For creation of gdx files:
$ echo jacobian hier_approach_2blocks_2by2.gdx > convertd.opt
  option lp=convertd;
  m.optfile = 1;
  solve m use lp min z;
$else
  option lp=cplex;
$ onecho > cplex.opt
  lpmethod 4
  solutiontype 2
$ offecho
  m.optfile = 1;
  solve m use lp min z;
$endif
