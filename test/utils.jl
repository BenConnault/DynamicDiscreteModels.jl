@test_approx_eq fatdiagonal((.7,.2),2) [.7 .3;.3 .7]
@test_approx_eq fatdiagonal((.7,.2),3) [.7 .15 .15;.15 .7 .15; .15 .15 .7 ]
@test_approx_eq fatdiagonal((.7,.2),4) [.7 .1 .1 .1;.1 .7 .1 .1; .1 .1 .7 .1; .1 .1 .1 .7]
@test_approx_eq fatdiagonal((.7,.2),5) [.7 .1 .1 .05 .05;.1 .7 .1 .05 .05; .05  .1 .7 .1 .05;.05 .05 .1 .7 .1; .05 .05 .1 .1 .7]

