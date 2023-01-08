<script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS_CHTML" ></script>

<script type="text/x-mathjax-config">
 MathJax.Hub.Config({
 tex2jax: {
 inlineMath: [["$","$"]]
 }
 });
</script>

`make_latex_text.py` の出力結果のプレビュー確認用

cell 1, 2, 3
$$
A^T =
\begin{bmatrix}
    A _{1, 1} & \cdots & A _{11, 1} \\
    \vdots & & \vdots \\
    A _{1, 5} & \cdots & A _{11, 5} \\
    A _{1, 6} & \cdots & A _{11, 6} \\
    A _{1, 7} & \cdots & A _{11, 7} \\
    \vdots & & \vdots \\
    A _{1, 11} & \cdots & A _{11, 11} \\
    A _{1, 12} & \cdots & A _{11, 12} \\
    \vdots & & \vdots \\
    A _{1, 16} & \cdots & A _{11, 16} \\
    A _{1, 17} & \cdots & A _{11, 17} \\
    A _{1, 18} & \cdots & A _{11, 18} \\
    \vdots & & \vdots \\
    A _{1, 22} & \cdots & A _{11, 22} \\
    A _{1, 23} & \cdots & A _{11, 23} \\
    \vdots & & \vdots \\
    A _{1, 27} & \cdots & A _{11, 27} \\
    A _{1, 28} & \cdots & A _{11, 28} \\
    A _{1, 29} & \cdots & A _{11, 29} \\
    \vdots & & \vdots \\
    A _{1, 33} & \cdots & A _{11, 33} \\
\end{bmatrix} =
\begin{bmatrix}
    a _{1, 1} & \cdots & a _{11, 1} \\
    \vdots & & \vdots \\
    a _{1, 5} & \cdots & a _{11, 5} \\
    a _{1, 6} - a _{1, 12} & \cdots & a _{11, 6} - a _{11, 12} \\
    a _{1, 7} & \cdots & a _{11, 7} \\
    \vdots & & \vdots \\
    a _{1, 11} & \cdots & a _{11, 11} \\
    a _{1, 13} & \cdots & a _{11, 13} \\
    \vdots & & \vdots \\
    a _{1, 17} & \cdots & a _{11, 17} \\
    a _{1, 18} - a _{1, 24} & \cdots & a _{11, 18} - a _{11, 24} \\
    a _{1, 19} & \cdots & a _{11, 19} \\
    \vdots & & \vdots \\
    a _{1, 23} & \cdots & a _{11, 23} \\
    a _{1, 25} & \cdots & a _{11, 25} \\
    \vdots & & \vdots \\
    a _{1, 29} & \cdots & a _{11, 29} \\
    a _{1, 30} - a _{1, 36} & \cdots & a _{11, 30} - a _{11, 36} \\
    a _{1, 31} & \cdots & a _{11, 31} \\
    \vdots & & \vdots \\
    a _{1, 35} & \cdots & a _{11, 35} \\
\end{bmatrix} =
\begin{bmatrix}
    \vec{a _{1}} \, ^T \\
    : \\
    \vec{a _{5}} \, ^T \\
    \vec{a _{6}} \, ^T - \vec{a _{12}} \, ^T \\
    \vec{a _{7}} \, ^T \\
    : \\
    \vec{a _{11}} \, ^T \\
    \vec{a _{13}} \, ^T \\
    : \\
    \vec{a _{17}} \, ^T \\
    \vec{a _{18}} \, ^T - \vec{a _{24}} \, ^T \\
    \vec{a _{19}} \, ^T \\
    : \\
    \vec{a _{23}} \, ^T \\
    \vec{a _{25}} \, ^T \\
    : \\
    \vec{a _{29}} \, ^T \\
    \vec{a _{30}} \, ^T - \vec{a _{36}} \, ^T \\
    \vec{a _{31}} \, ^T \\
    : \\
    \vec{a _{35}} \, ^T \\
\end{bmatrix}
$$

cell 4, 5
$$
\vec{x} =
\begin{bmatrix}
    x _{1} \\
    \vdots \\
    x _{5} \\
    x _{6} \\
    x _{7} \\
    \vdots \\
    x _{11} \\
    x _{12} \\
    \vdots \\
    x _{16} \\
    x _{17} \\
    x _{18} \\
    \vdots \\
    x _{22} \\
    x _{23} \\
    \vdots \\
    x _{27} \\
    x _{28} \\
    x _{29} \\
    \vdots \\
    x _{33} \\
\end{bmatrix} =
\begin{bmatrix}
    y _{1} \\
    \vdots \\
    y _{5} \\
    (y _{6} + y _{12}) / 2 \\
    y _{7} \\
    \vdots \\
    y _{11} \\
    y _{13} \\
    \vdots \\
    y _{17} \\
    (y _{18} + y _{24}) / 2 \\
    y _{19} \\
    \vdots \\
    y _{23} \\
    y _{25} \\
    \vdots \\
    y _{29} \\
    (y _{30} + y _{36}) / 2 \\
    y _{31} \\
    \vdots \\
    y _{35} \\
\end{bmatrix}
$$

cell 7
$$
\begin{bmatrix}
    x _{1} \\
    x _{2} \\
    x _{3} \\
    x _{4} \\
    x _{5} \\
    x _{6} \\
    x _{7} \\
    x _{8} \\
    x _{9} \\
    x _{10} \\
    x _{11} \\
    \vdots \\
\end{bmatrix} =
\begin{bmatrix}
    y _{1} \\
    y _{2} \\
    y _{3} \\
    y _{4} \\
    y _{5} \\
    (y _{6} + y _{12}) / 2 \\
    y _{7} \\
    y _{8} \\
    y _{9} \\
    y _{10} \\
    y _{11} \\
    \vdots \\
\end{bmatrix}
$$

cell 6

| $x _{1}$ | $\cdots$ | $x _{5}$ | $x _{6}$ | $x _{7}$ | $\cdots$ | $x _{11}$ | $x _{12}$ | $\cdots$ | $x _{16}$ | $x _{17}$ | $x _{18}$ | $\cdots$ | $x _{22}$ | $x _{23}$ | $\cdots$ | $x _{27}$ | $x _{28}$ | $x _{29}$ | $\cdots$ | $x _{33}$ |
| - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - |
| $y _{1}$ | $\cdots$ | $y _{5}$ | $(y _{6} + y _{12}) / 2$ | $y _{7}$ | $\cdots$ | $y _{11}$ | $y _{13}$ | $\cdots$ | $y _{17}$ | $(y _{18} + y _{24}) / 2$ | $y _{19}$ | $\cdots$ | $y _{23}$ | $y _{25}$ | $\cdots$ | $y _{29}$ | $(y _{30} + y _{36}) / 2$ | $y _{31}$ | $\cdots$ | $y _{35}$ |

| $x _{1}$ | $x _{2}$ | $x _{3}$ | $x _{4}$ | $x _{5}$ | $x _{6}$ | $x _{7}$ | $x _{8}$ | $x _{9}$ | $x _{10}$ | $x _{11}$ |
| - | - | - | - | - | - | - | - | - | - | - |
| $y _{1}$ | $y _{2}$ | $y _{3}$ | $y _{4}$ | $y _{5}$ | $y _{6} + y _{12}$ | $y _{7}$ | $y _{8}$ | $y _{9}$ | $y _{10}$ | $y _{11}$ |

cell 8
$$
\begin{bmatrix}
    \vec{\vphantom{d} a _{1}} \\
    \vec{\vphantom{d} a _{2}} \\
    \vec{\vphantom{d} a _{3}} \\
    \vec{\vphantom{d} a _{4}} \\
    \vec{\vphantom{d} a _{5}} \\
    \vec{\vphantom{d} a _{6}} \\
    \vec{\vphantom{d} a _{7}} \\
    \vec{\vphantom{d} a _{8}} \\
    \vec{\vphantom{d} a _{9}} \\
    \vec{\vphantom{d} a _{10}} \\
    \vec{\vphantom{d} a _{11}} \\
    \vdots \\
\end{bmatrix} =
\begin{bmatrix}
    \mathcal{Z}(\vec{d _{1}}) \\
    \mathcal{Z}(\vec{d _{2}}) \\
    \mathcal{Z}(\vec{d _{3}}) \\
    \mathcal{Z}(\vec{d _{4}}) \\
    \mathcal{Z}(\vec{d _{5}}) \\
    \mathcal{Z}(\vec{d _{6}} - \vec{d _{12}}) \\
    \mathcal{Z}(\vec{d _{7}}) \\
    \mathcal{Z}(\vec{d _{8}}) \\
    \mathcal{Z}(\vec{d _{9}}) \\
    \mathcal{Z}(\vec{d _{10}}) \\
    \mathcal{Z}(\vec{d _{11}}) \\
    \vdots \\
\end{bmatrix}
$$

cell 9

| $\vec{a _{1}}$ | $\cdots$ | $\vec{a _{5}}$ | $\vec{a _{6}}$ | $\vec{a _{7}}$ | $\cdots$ | $\vec{a _{11}}$ | $\vec{a _{12}}$ | $\cdots$ | $\vec{a _{16}}$ | $\vec{a _{17}}$ | $\vec{a _{18}}$ | $\cdots$ | $\vec{a _{22}}$ | $\vec{a _{23}}$ | $\cdots$ | $\vec{a _{27}}$ | $\vec{a _{28}}$ | $\vec{a _{29}}$ | $\cdots$ | $\vec{a _{33}}$ |
| - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - |
| $\mathcal{Z}( \vec{d _{1}})$ | $\cdots$ | $\mathcal{Z}( \vec{d _{5}})$ | $\mathcal{Z}( \vec{d _{6}} - \vec{d _{12}})$ | $\mathcal{Z}( \vec{d _{7}})$ | $\cdots$ | $\mathcal{Z}( \vec{d _{11}})$ | $\mathcal{Z}( \vec{d _{13}})$ | $\cdots$ | $\mathcal{Z}( \vec{d _{17}})$ | $\mathcal{Z}( \vec{d _{18}} - \vec{d _{24}})$ | $\mathcal{Z}( \vec{d _{19}})$ | $\cdots$ | $\mathcal{Z}( \vec{d _{23}})$ | $\mathcal{Z}( \vec{d _{25}})$ | $\cdots$ | $\mathcal{Z}( \vec{d _{29}})$ | $\mathcal{Z}( \vec{d _{30}} - \vec{d _{36}})$ | $\mathcal{Z}( \vec{d _{31}})$ | $\cdots$ | $\mathcal{Z}( \vec{d _{35}})$ |

cell 10

||||
| - | - | - |
| $x _{1}$ = $y _{1}$ | $x _{12}$ = $y _{13}$ | $x _{23}$ = $y _{25}$ |
| $x _{2}$ = $y _{2}$ | $x _{13}$ = $y _{14}$ | $x _{24}$ = $y _{26}$ |
| $x _{3}$ = $y _{3}$ | $x _{14}$ = $y _{15}$ | $x _{25}$ = $y _{27}$ |
| $x _{4}$ = $y _{4}$ | $x _{15}$ = $y _{16}$ | $x _{26}$ = $y _{28}$ |
| $x _{5}$ = $y _{5}$ | $x _{16}$ = $y _{17}$ | $x _{27}$ = $y _{29}$ |
| $x _{6}$ = $(y _{6} + y _{12}) / 2$ | $x _{17}$ = $(y _{18} + y _{24}) / 2$ | $x _{28}$ = $(y _{30} + y _{36}) / 2$ |
| $x _{7}$ = $y _{7}$ | $x _{18}$ = $y _{19}$ | $x _{29}$ = $y _{31}$ |
| $x _{8}$ = $y _{8}$ | $x _{19}$ = $y _{20}$ | $x _{30}$ = $y _{32}$ |
| $x _{9}$ = $y _{9}$ | $x _{20}$ = $y _{21}$ | $x _{31}$ = $y _{33}$ |
| $x _{10}$ = $y _{10}$ | $x _{21}$ = $y _{22}$ | $x _{32}$ = $y _{34}$ |
| $x _{11}$ = $y _{11}$ | $x _{22}$ = $y _{23}$ | $x _{33}$ = $y _{35}$ |

cell 11
||||
| - | - | - |
| $\vec{a _{1}} = \mathcal{Z}(\vec{d _{1}})$ | $\vec{a _{12}} = \mathcal{Z}(\vec{d _{13}})$ | $\vec{a _{23}} = \mathcal{Z}(\vec{d _{25}})$ |
| $\vec{a _{2}} = \mathcal{Z}(\vec{d _{2}})$ | $\vec{a _{13}} = \mathcal{Z}(\vec{d _{14}})$ | $\vec{a _{24}} = \mathcal{Z}(\vec{d _{26}})$ |
| $\vec{a _{3}} = \mathcal{Z}(\vec{d _{3}})$ | $\vec{a _{14}} = \mathcal{Z}(\vec{d _{15}})$ | $\vec{a _{25}} = \mathcal{Z}(\vec{d _{27}})$ |
| $\vec{a _{4}} = \mathcal{Z}(\vec{d _{4}})$ | $\vec{a _{15}} = \mathcal{Z}(\vec{d _{16}})$ | $\vec{a _{26}} = \mathcal{Z}(\vec{d _{28}})$ |
| $\vec{a _{5}} = \mathcal{Z}(\vec{d _{5}})$ | $\vec{a _{16}} = \mathcal{Z}(\vec{d _{17}})$ | $\vec{a _{27}} = \mathcal{Z}(\vec{d _{29}})$ |
| $\vec{a _{6}} = \mathcal{Z}(\vec{d _{6}} - \vec{d _{12}})$ | $\vec{a _{17}} = \mathcal{Z}(\vec{d _{18}} - \vec{d _{24}})$ | $\vec{a _{28}} = \mathcal{Z}(\vec{d _{30}} - \vec{d _{36}})$ |
| $\vec{a _{7}} = \mathcal{Z}(\vec{d _{7}})$ | $\vec{a _{18}} = \mathcal{Z}(\vec{d _{19}})$ | $\vec{a _{29}} = \mathcal{Z}(\vec{d _{31}})$ |
| $\vec{a _{8}} = \mathcal{Z}(\vec{d _{8}})$ | $\vec{a _{19}} = \mathcal{Z}(\vec{d _{20}})$ | $\vec{a _{30}} = \mathcal{Z}(\vec{d _{32}})$ |
| $\vec{a _{9}} = \mathcal{Z}(\vec{d _{9}})$ | $\vec{a _{20}} = \mathcal{Z}(\vec{d _{21}})$ | $\vec{a _{31}} = \mathcal{Z}(\vec{d _{33}})$ |
| $\vec{a _{10}} = \mathcal{Z}(\vec{d _{10}})$ | $\vec{a _{21}} = \mathcal{Z}(\vec{d _{22}})$ | $\vec{a _{32}} = \mathcal{Z}(\vec{d _{34}})$ |
| $\vec{a _{11}} = \mathcal{Z}(\vec{d _{11}})$ | $\vec{a _{22}} = \mathcal{Z}(\vec{d _{23}})$ | $\vec{a _{33}} = \mathcal{Z}(\vec{d _{35}})$ |

cel 12
||||
| - | - | - |
| $x _{1} = 1 \dashrightarrow \vec{a _{1}} = \mathcal{Z}(\vec{d _{1}})$ | $x _{12} = 1 \dashrightarrow \vec{a _{12}} = \mathcal{Z}(\vec{d _{13}})$ | $x _{23} = 1 \dashrightarrow \vec{a _{23}} = \mathcal{Z}(\vec{d _{25}})$ |
| $x _{2} = 1 \dashrightarrow \vec{a _{2}} = \mathcal{Z}(\vec{d _{2}})$ | $x _{13} = 1 \dashrightarrow \vec{a _{13}} = \mathcal{Z}(\vec{d _{14}})$ | $x _{24} = 1 \dashrightarrow \vec{a _{24}} = \mathcal{Z}(\vec{d _{26}})$ |
| $x _{3} = 1 \dashrightarrow \vec{a _{3}} = \mathcal{Z}(\vec{d _{3}})$ | $x _{14} = 1 \dashrightarrow \vec{a _{14}} = \mathcal{Z}(\vec{d _{15}})$ | $x _{25} = 1 \dashrightarrow \vec{a _{25}} = \mathcal{Z}(\vec{d _{27}})$ |
| $x _{4} = 1 \dashrightarrow \vec{a _{4}} = \mathcal{Z}(\vec{d _{4}})$ | $x _{15} = 1 \dashrightarrow \vec{a _{15}} = \mathcal{Z}(\vec{d _{16}})$ | $x _{26} = 1 \dashrightarrow \vec{a _{26}} = \mathcal{Z}(\vec{d _{28}})$ |
| $x _{5} = 1 \dashrightarrow \vec{a _{5}} = \mathcal{Z}(\vec{d _{5}})$ | $x _{16} = 1 \dashrightarrow \vec{a _{16}} = \mathcal{Z}(\vec{d _{17}})$ | $x _{27} = 1 \dashrightarrow \vec{a _{27}} = \mathcal{Z}(\vec{d _{29}})$ |
| $x _{6} = 1 \dashrightarrow \vec{a _{6}} = \mathcal{Z}(\vec{d _{6}} - \vec{d _{12}})$ | $x _{17} = 1 \dashrightarrow \vec{a _{17}} = \mathcal{Z}(\vec{d _{18}} - \vec{d _{24}})$ | $x _{28} = 1 \dashrightarrow \vec{a _{28}} = \mathcal{Z}(\vec{d _{30}} - \vec{d _{36}})$ |
| $x _{7} = 1 \dashrightarrow \vec{a _{7}} = \mathcal{Z}(\vec{d _{7}})$ | $x _{18} = 1 \dashrightarrow \vec{a _{18}} = \mathcal{Z}(\vec{d _{19}})$ | $x _{29} = 1 \dashrightarrow \vec{a _{29}} = \mathcal{Z}(\vec{d _{31}})$ |
| $x _{8} = 1 \dashrightarrow \vec{a _{8}} = \mathcal{Z}(\vec{d _{8}})$ | $x _{19} = 1 \dashrightarrow \vec{a _{19}} = \mathcal{Z}(\vec{d _{20}})$ | $x _{30} = 1 \dashrightarrow \vec{a _{30}} = \mathcal{Z}(\vec{d _{32}})$ |
| $x _{9} = 1 \dashrightarrow \vec{a _{9}} = \mathcal{Z}(\vec{d _{9}})$ | $x _{20} = 1 \dashrightarrow \vec{a _{20}} = \mathcal{Z}(\vec{d _{21}})$ | $x _{31} = 1 \dashrightarrow \vec{a _{31}} = \mathcal{Z}(\vec{d _{33}})$ |
| $x _{10} = 1 \dashrightarrow \vec{a _{10}} = \mathcal{Z}(\vec{d _{10}})$ | $x _{21} = 1 \dashrightarrow \vec{a _{21}} = \mathcal{Z}(\vec{d _{22}})$ | $x _{32} = 1 \dashrightarrow \vec{a _{32}} = \mathcal{Z}(\vec{d _{34}})$ |
| $x _{11} = 1 \dashrightarrow \vec{a _{11}} = \mathcal{Z}(\vec{d _{11}})$ | $x _{22} = 1 \dashrightarrow \vec{a _{22}} = \mathcal{Z}(\vec{d _{23}})$ | $x _{33} = 1 \dashrightarrow \vec{a _{33}} = \mathcal{Z}(\vec{d _{35}})$ |

cell 13
||||
| - | - | - |
| $x _{1} \left ( = y _{1}\right ) = 1 \dashrightarrow \vec{a _{1}} = \mathcal{Z}(\vec{d _{1}})$ | $x _{12} \left ( = y _{13}\right ) = 1 \dashrightarrow \vec{a _{12}} = \mathcal{Z}(\vec{d _{13}})$ | $x _{23} \left ( = y _{25}\right ) = 1 \dashrightarrow \vec{a _{23}} = \mathcal{Z}(\vec{d _{25}})$ |
| $x _{2} \left ( = y _{2}\right ) = 1 \dashrightarrow \vec{a _{2}} = \mathcal{Z}(\vec{d _{2}})$ | $x _{13} \left ( = y _{14}\right ) = 1 \dashrightarrow \vec{a _{13}} = \mathcal{Z}(\vec{d _{14}})$ | $x _{24} \left ( = y _{26}\right ) = 1 \dashrightarrow \vec{a _{24}} = \mathcal{Z}(\vec{d _{26}})$ |
| $x _{3} \left ( = y _{3}\right ) = 1 \dashrightarrow \vec{a _{3}} = \mathcal{Z}(\vec{d _{3}})$ | $x _{14} \left ( = y _{15}\right ) = 1 \dashrightarrow \vec{a _{14}} = \mathcal{Z}(\vec{d _{15}})$ | $x _{25} \left ( = y _{27}\right ) = 1 \dashrightarrow \vec{a _{25}} = \mathcal{Z}(\vec{d _{27}})$ |
| $x _{4} \left ( = y _{4}\right ) = 1 \dashrightarrow \vec{a _{4}} = \mathcal{Z}(\vec{d _{4}})$ | $x _{15} \left ( = y _{16}\right ) = 1 \dashrightarrow \vec{a _{15}} = \mathcal{Z}(\vec{d _{16}})$ | $x _{26} \left ( = y _{28}\right ) = 1 \dashrightarrow \vec{a _{26}} = \mathcal{Z}(\vec{d _{28}})$ |
| $x _{5} \left ( = y _{5}\right ) = 1 \dashrightarrow \vec{a _{5}} = \mathcal{Z}(\vec{d _{5}})$ | $x _{16} \left ( = y _{17}\right ) = 1 \dashrightarrow \vec{a _{16}} = \mathcal{Z}(\vec{d _{17}})$ | $x _{27} \left ( = y _{29}\right ) = 1 \dashrightarrow \vec{a _{27}} = \mathcal{Z}(\vec{d _{29}})$ |
| $x _{6} \left ( = (y _{6} + y _{12}) / 2\right ) = 1 \dashrightarrow \vec{a _{6}} = \mathcal{Z}(\vec{d _{6}} - \vec{d _{12}})$ | $x _{17} \left ( = (y _{18} + y _{24}) / 2\right ) = 1 \dashrightarrow \vec{a _{17}} = \mathcal{Z}(\vec{d _{18}} - \vec{d _{24}})$ | $x _{28} \left ( = (y _{30} + y _{36}) / 2\right ) = 1 \dashrightarrow \vec{a _{28}} = \mathcal{Z}(\vec{d _{30}} - \vec{d _{36}})$ |
| $x _{7} \left ( = y _{7}\right ) = 1 \dashrightarrow \vec{a _{7}} = \mathcal{Z}(\vec{d _{7}})$ | $x _{18} \left ( = y _{19}\right ) = 1 \dashrightarrow \vec{a _{18}} = \mathcal{Z}(\vec{d _{19}})$ | $x _{29} \left ( = y _{31}\right ) = 1 \dashrightarrow \vec{a _{29}} = \mathcal{Z}(\vec{d _{31}})$ |
| $x _{8} \left ( = y _{8}\right ) = 1 \dashrightarrow \vec{a _{8}} = \mathcal{Z}(\vec{d _{8}})$ | $x _{19} \left ( = y _{20}\right ) = 1 \dashrightarrow \vec{a _{19}} = \mathcal{Z}(\vec{d _{20}})$ | $x _{30} \left ( = y _{32}\right ) = 1 \dashrightarrow \vec{a _{30}} = \mathcal{Z}(\vec{d _{32}})$ |
| $x _{9} \left ( = y _{9}\right ) = 1 \dashrightarrow \vec{a _{9}} = \mathcal{Z}(\vec{d _{9}})$ | $x _{20} \left ( = y _{21}\right ) = 1 \dashrightarrow \vec{a _{20}} = \mathcal{Z}(\vec{d _{21}})$ | $x _{31} \left ( = y _{33}\right ) = 1 \dashrightarrow \vec{a _{31}} = \mathcal{Z}(\vec{d _{33}})$ |
| $x _{10} \left ( = y _{10}\right ) = 1 \dashrightarrow \vec{a _{10}} = \mathcal{Z}(\vec{d _{10}})$ | $x _{21} \left ( = y _{22}\right ) = 1 \dashrightarrow \vec{a _{21}} = \mathcal{Z}(\vec{d _{22}})$ | $x _{32} \left ( = y _{34}\right ) = 1 \dashrightarrow \vec{a _{32}} = \mathcal{Z}(\vec{d _{34}})$ |
| $x _{11} \left ( = y _{11}\right ) = 1 \dashrightarrow \vec{a _{11}} = \mathcal{Z}(\vec{d _{11}})$ | $x _{22} \left ( = y _{23}\right ) = 1 \dashrightarrow \vec{a _{22}} = \mathcal{Z}(\vec{d _{23}})$ | $x _{33} \left ( = y _{35}\right ) = 1 \dashrightarrow \vec{a _{33}} = \mathcal{Z}(\vec{d _{35}})$ |