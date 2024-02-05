

In [213]: subplot(311)
Out[213]: <matplotlib.axes.Subplot instance at 0x88bae18>

In [214]: ylabel("R (um)")
Out[214]: <matplotlib.text.Text instance at 0x88c7b90>

In [215]: pcolor(zar*1e2,Rar*1e6,chiar2d)
Out[215]: <matplotlib.collections.PolyCollection instance at 0x8894b48>

In [216]: subplot(312)
Out[216]: <matplotlib.axes.Subplot instance at 0x8894dd0>

In [217]: plot(zar*1e2,T/T.max(),zar*1e2,P/P.max(),zar*1e2,M)
Out[217]:
    [<matplotlib.lines.Line2D instance at 0x8896d40>,
      <matplotlib.lines.Line2D instance at 0x8896d88>,
      <matplotlib.lines.Line2D instance at 0x8896ab8>]

    In [218]: ylabel("T/T0, P/P0, M")
    Out[218]: <matplotlib.text.Text instance at 0x8898c68>

    In [219]: legend(("T/T0", "P/P0", "M"),loc='best')
    Out[219]: <matplotlib.legend.Legend instance at 0x8896ef0>

    In [220]: subplot(313)
    Out[220]: <matplotlib.axes.Subplot instance at 0x8896f38>

    In [221]: plot(zar*1e2,log10(totalions/totalions.max()))
    Out[221]: [<matplotlib.lines.Line2D instance at 0x8897f38>]

    In [222]: xlabel("z (cm)")
    Out[222]: <matplotlib.text.Text instance at 0x8899368>

    In [223]: ylabel("log10(qty)")
    Out[223]: <matplotlib.text.Text instance at 0x8899f38>

    In [224]: plot(zar*1e2,log10(totalions/totalions.max()))
    Out[224]: [<matplotlib.lines.Line2D instance at 0x889d3f8>]

    In [225]: legend(("", "normed total ion flux"),loc='best')
    Out[225]: <matplotlib.legend.Legend instance at 0x889d440>

    In [226]: show()

    


subplot(311)
title("Ion concentration versus z and R")
xlabel("z (cm)")
ylabel("R (um)")
pcolor(zar,Rar*1e6,chiar2d)

subplot(312)
plot(zar,T/T.max(),zar,P/P.max(),zar,M)


subplot(313)
plot(zar,totalions)

        
################################################################################

newchi = ones([65,10000])

In [263]: for i in 65:
    .....:
        .....:
            .....:
                .....:       print i
                .....:
                    .....:
                        ---------------------------------------------------------------------------
                        TypeError                                 Traceback (most recent call last)

                        /home/pwilliam/projects/capillary/<ipython console> in <module>()

                        TypeError: 'int' object is not iterable

                        In [264]:

                            In [265]:

                                In [266]: for i in range(65):
                                    .....:     for j in range(10000):
                                        .....:         z = i + j
                                        .....:

                                            In [267]: z
                                            Out[267]: 10063

                                            In [268]:

                                                In [269]:

                                                    In [270]: for j in range(10000):
                                                        .....:     newchi[32,j] = chiar2d[0,j]
                                                        .....:

                                                            In [271]: shape(chiar2d)
                                                            Out[271]: (32, 10000)

                                                            In [272]: chiar[31,0]
                                                            ---------------------------------------------------------------------------
                                                            IndexError                                Traceback (most recent call last)

                                                            /home/pwilliam/projects/capillary/<ipython console> in <module>()

                                                            IndexError: invalid index

                                                            In [273]: chiar2d[31,0]
                                                            Out[273]: 1.0

                                                            In [274]: chiar2d[32,0]
                                                            ---------------------------------------------------------------------------
                                                            IndexError                                Traceback (most recent call last)

                                                            /home/pwilliam/projects/capillary/<ipython console> in <module>()

                                                            IndexError: index (32) out of range (0<=index<=32) in dimension 0

                                                            In [275]:

                                                                In [276]:

                                                                    In [277]: for i in range(32):
                                                                        .....:     for j in range(10000):
                                                                            .....:         z = i + j
                                                                            .....:

                                                                                In [278]: for i in range(1,32):
                                                                                    .....:     for j in range(10000):
                                                                                        .....:         newchi[32 + i,j] = chiar2d[i,j]
                                                                                        .....:         newchi[32 - i,j] = chiar2d[i,j]
                                                                                        .....:

                                                                                            In [279]: pcolor(newchi)
                                                                                            Out[279]: <matplotlib.collections.PolyCollection instance at 0x10bd0fc8>

                                                                                            In [280]: show()

                                                                                            In [281]: newchi[0,54]
                                                                                            Out[281]: 1.0

                                                                                            In [282]: for j in range(10000):
                                                                                                .....:     newchi[0,j] = 0.0
                                                                                                .....:     newchi[64,j] = 0.0
                                                                                                .....:

                                                                                                    In [283]: Rar
                                                                                                    Out[283]:
                                                                                                        array([  0.00000000e+00,   9.52620968e-06,   1.90524194e-05,
                                                                                                                          2.85786290e-05,   3.81048387e-05,   4.76310484e-05,
                                                                                                                          5.71572581e-05,   6.66834677e-05,   7.62096774e-05,
                                                                                                                          8.57358871e-05,   9.52620968e-05,   1.04788306e-04,
                                                                                                                          1.14314516e-04,   1.23840726e-04,   1.33366935e-04,
                                                                                                                          1.42893145e-04,   1.52419355e-04,   1.61945565e-04,
                                                                                                                          1.71471774e-04,   1.80997984e-04,   1.90524194e-04,
                                                                                                                          2.00050403e-04,   2.09576613e-04,   2.19102823e-04,
                                                                                                                          2.28629032e-04,   2.38155242e-04,   2.47681452e-04,
                                                                                                                          2.57207661e-04,   2.66733871e-04,   2.76260081e-04,
                                                                                                                          2.85786290e-04,   2.95312500e-04])

                                                                                                        In [284]: shape(Rar)
                                                                                                        Out[284]: (32,)

                                                                                                        In [285]: newrar = zeros(65)

                                                                                                        In [286]: newrar[5]
                                                                                                        Out[286]: 0.0

                                                                                                        In [287]: for i in range(size(Rar)):
                                                                                                            .....:     i += 1
                                                                                                            .....:

                                                                                                                In [288]: for i in range(1,size(Rar)):
                                                                                                                    .....:     newrar[32 + i] = Rar[i]
                                                                                                                    .....:     newrar[32 - i] = -Rar[i]
                                                                                                                    .....:

                                                                                                                        In [289]: newrar[32] = 0.0

                                                                                                                        In [290]: newrar[64] = 300.e-6

                                                                                                                        In [291]: newrar[0] = -300.e-6


                                                                                                                        
