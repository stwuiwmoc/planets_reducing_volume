<script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS_CHTML" ></script>

<script type="text/x-mathjax-config">
 MathJax.Hub.Config({
 tex2jax: {
 inlineMath: [["$","$"]]
 }
 });
</script>

- [zernike多項式](#zernike多項式)
- [前提](#前提)
- [モーター駆動量ベクトル](#モーター駆動量ベクトル)
- [最適モーター駆動量ベクトルの導出](#最適モーター駆動量ベクトルの導出)
- [WH駆動量ベクトルの導入](#wh駆動量ベクトルの導入)
- [制約付き](#制約付き)
- [行列置き場](#行列置き場)

# zernike多項式

$$
W(\rho , \theta) = \sum _j a _j Z _j (\rho , \theta)
$$

ここで、 $a _j$ はzernike 係数と呼ばれる。 Zernike多項式の各項 $Z _j (\rho , \theta)$ は

$$
Z _j (\rho , \theta) =
\left\{
    \begin{array}{ll}
        \sqrt{2 (n + 1)} R _n ^m (\rho) \cos{(m \theta )} & (m \neq 0 \: \& \: j = even) \\
        \sqrt{2 (n + 1)} R _n ^m (\rho) \sin{(m \theta )} & (m \neq 0 \: \& \: j = odd) \\
        \sqrt{(n + 1)} R _n ^0 (\rho) & (m = 0)
    \end{array}
\right.
$$

$$
R _n ^m (\rho) =
\sum _{s = 0} ^{(n - m) / 2}
\cfrac
    {(-1) ^s (n - s) !}
    {s ! \left (  \cfrac{n + m}{2} - s \right ) ! \left (  \cfrac{n - m}{2} - s \right ) !}
\rho ^{n - 2s}
$$

$() _{int}$ はカッコ内の整数部分を表す

$$
n = [(2j - 1) ^{1/2} + 0.5] _{int} - 1
$$

$$
m =
\left\{
    \begin{array}{ll}
        2 \{ [2j + 1 - n(n + 1)] / 4 \} _{int} & (n = even) \\
        2 \{ [2(j + 1) - n(n + 1)] / 4 \} _{int} - 1 & (n = odd)
    \end{array}
\right.
$$

# 前提

以下ではスカラーを小文字 $a$ 、ベクトルを矢印記法 $\vec{a}$ 、行列を大文字 $A$ で表記する。ベクトルは全て縦ベクトルであるが、幅をとるので文中では転置して $\vec{a} ^T = (a_1, \cdots , a_n)$ のように表記する。

$$
\vec{a} =
\begin{bmatrix}
    a_1 \\
    \vdots \\
    a_n
\end{bmatrix}
$$

# モーター駆動量ベクトル

鏡面上に均等に配置された $m$ 個の点を考え、各点での鏡面の高さ方向の変形量 $f _1, f _2, \cdots, f _{m}$  [m] を要素とする変形量ベクトル $\vec{f}$ を $\vec{f} ^T \equiv (f _1 , \cdots, f _{m})$ と定義する。

また、WHのM01, M02, ... , M36 のモーターの駆動量 $y _1, y _2, \cdots, y _{36}$ [mm] を要素とするモーター駆動量ベクトル $\vec{y}$ を $\vec{y} ^T \equiv (y _1 , \cdots , y _{36})$ と定義する。

任意のモーター駆動量ベクトル $\vec{y}$ を与えた時の鏡面の変形量ベクトル $\vec{f}$ を計算するためのモーター駆動量変換行列 $D$ （以後、作用行列と呼ぶ）を以下の式を満たすように定義する。

$$
\vec{f} = D \vec{y}
$$
$$
\begin{bmatrix}
    f _1 \\
    \vdots \\
    f _{m}
\end{bmatrix} =
\begin{bmatrix}
    D _{1, 1} & \cdots & D _{1, 36} \\
    \vdots & \ddots & \vdots \\
    D _{m, 1} & \cdots & D _{m, 36}
\end{bmatrix}
\begin{bmatrix}
    y _1 \\
    \vdots \\
    y _{36}
\end{bmatrix}
$$

よって、$D$ の各要素の値が事前に分かっていれば、任意の駆動量 $\vec{y}$ に対して変形量 $\vec{f}$ を得ることができる。

ここで、 $\vec{y} ^T = (1, 0, \cdots, 0)$ のときの変形量 $\vec{d _1}$ を考える。 $\vec{d _1}$ は $D$ を用いて

$$
\vec{d _1} =
\begin{bmatrix}
    D _{1, 1} & \cdots & D _{1, 36} \\
    \vdots & \ddots & \vdots \\
    D _{m, 1} & \cdots & D _{m, 36}
\end{bmatrix}
\begin{bmatrix}
    1 \\
    0 \\
    \vdots \\
    0
\end{bmatrix} =
\begin{bmatrix}
    D _{1, 1} \\
    \vdots \\
    D _{m, 1}
\end{bmatrix}
$$

と計算できるため、 $\vec{y} ^T = (1, 0, \cdots, 0)$ に対する変形量は $\vec{d _1} ^T = (D _{1, 1}, \cdots, D _{m, 1})$ である。同様に、 $\vec{y} ^T = (0, 1, \cdots, 0)$ に対する変形量は $\vec{d _2} ^T = (D _{1, 2}, \cdots, D _{m, 2})$ であり、これを $\vec{d _{36}}$ まで繰り返すと

$$
D =
\begin{bmatrix}
    D _{1, 1} & \cdots & D _{1, 36} \\
    \vdots & \ddots & \vdots \\
    D _{m, 1} & \cdots & D _{m, 36}
\end{bmatrix}
= \left [ \vec{d _1}, \cdots, \vec{d _{36}} \right ]
$$

となるため、 $D$ の各要素を求めるには $\vec{d _1}, \cdots, \vec{d _{36}}$ が分かれば良い。

ここで、 $\vec{d _1}$の物理的意味は、M01のモーターの駆動量 $y _1 = 1$ [mm]、他の駆動量 $y _2, \cdots, y _{36} = 0$ [mm] の時の変形量であり、以下に述べる方法で決定することができる。

まず、WHが駆動していないときの各点の高さ $\vec{h _0} ^T = (h _{1, 0}, \cdots, h _{m, 0})$ [m] を有限要素解析（FEM）で求める。

次に、act1に 1 [Nm] のトルクを与えたときの各点の高さ $\vec{h _1} ^T = (h _{1, 1}, \cdots, h _{m, 1})$ [m] をFEMで求める。

このとき、act1 に 1 [Nm] のトルクを与えたときの各点の変形量は $\vec{h _1} - \vec{h _0}$ [m] である。

M01が 1 [mm] 動いた時の変形量 $\vec{d _1}$ は、係数 $\alpha _1$ を用いて

$$
\vec{d _1} = \alpha _1 \left( \vec{h _1} - \vec{h _0} \right)
$$

ここで $\alpha _1$ は、act1 に 0.05 [Nm] のトルクがかかったときの変形量を、M01が 1 [mm] 動いたときの変形量に換算するための係数（倍率）である。

$D$ を改めて書き下すと

$$
D =
\begin{bmatrix}
    D _{1, 1} & \cdots & D _{1, 36} \\
    \vdots & \ddots & \vdots \\
    D _{m, 1} & \cdots & D _{m, 36}
\end{bmatrix} =
\begin{bmatrix}
    \alpha _1 (h _{1, 1}  - h _{1, 0}) & \cdots & \alpha _{36} (h _{1, 36}  - h _{1, 0})  \\
    \vdots & \ddots & \vdots \\
    \alpha _1 (h _{m, 1}  - h _{m, 0}) & \cdots & \alpha _{36} (h _{m, 36}  - h _{m, 0})  \\
\end{bmatrix}
$$

以上より、係数 $\alpha _1, \cdots, \alpha _{36}$ とFEMの結果があれば、任意のモーター駆動量ベクトル $\vec{y}$ に対して変形量ベクトル $\vec{f}$ を計算できる。

# 最適モーター駆動量ベクトルの導出

$$
\vec{f} = D \vec{y}
$$

上式において、既知の鏡面形状誤差 $\vec{f}$ を最小二乗的に再現（reproduction）するような最適モーター駆動量ベクトル $\vec{y_r}$ を考える。このような $\vec{y _r}$ を求めることができれば、既知の $\vec{f}$ を能動支持機構で最小二乗的に再現した結果である $\vec{f_r}$ を

$$
\vec{f_r} = D \vec{y_r}
$$

と計算できる。よって、能動支持機構で $\vec{f _r}$ を作って鏡面形状誤差 $\vec{f}$ を打ち消すことで、鏡面形状誤差を補正することが可能となる。

ここで、 $\vec{y _r}$ を求めることは

$$
\| D \vec{y} - \vec{f}\|
$$

を最小化するような $\vec{y_r}$ を求めることと同義である。 $D$ の要素である $\vec{d _1}, \vec{d _2}, \cdots, \vec{d _{36}}$ が互いに独立な場合は、このような $\vec{y_r}$ は正規方程式の解としてただ一つに定まることが知られている。

PLANETSのWTとWHでは、M06とM12はどちらも同じ回転軸に沿ったトルクをかける構造となっており、 $\vec{d _6}$ と $\vec{d _{12}}$ はトルクの正負以外は同一である。よって、例えば M06 を1 [mm] 駆動させても、M12を -1 [mm] 駆動させても同じ形状が得られる。figureでは形状が異なるが、これは本研究で使用したFEMが不完全であることに由来しており、piston, tip, tiltを除いたfigureで比較すると正負以外は全く同一な形状であることが確認できる。これはM18とM24、M30とM36でも同様である。

よって、 $D$ の要素 $\vec{d _n}$ が部分的に縮退している（独立ではない）ため、この重複を取り除いた新たな行列を作る必要がある。

# WH駆動量ベクトルの導入

新たにWH駆動量ベクトル $\vec{x} ^T \equiv (x _1, \cdots, x _{33})$ を定義する。これは、 $\vec{y}$ のうち縮退が生じているモーターの駆動量である $y_{6}$ と $y_{12}$ 、  $y_{18}$ と $y_{24}$ 、 $y_{30}$ と $y_{36}$ をそれぞれまとめたものである。以下に $x_1 \sim x_{11}$ と $y_1 \sim y_{12}$ の対応関係を示す。

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

$x_{12} \sim x_{22}$ と $y_{13} \sim y_{24}$ 、 $x_{23} \sim x_{33}$ と $y_{25} \sim y_{36}$ においても同様の対応関係がある。

| $x _{1}$ | $\cdots$ | $x _{5}$ | $x _{6}$ | $x _{7}$ | $\cdots$ | $x _{11}$ | $x _{12}$ | $\cdots$ | $x _{16}$ | $x _{17}$ | $x _{18}$ | $\cdots$ | $x _{22}$ | $x _{23}$ | $\cdots$ | $x _{27}$ | $x _{28}$ | $x _{29}$ | $\cdots$ | $x _{33}$ |
| - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - |
| $y _{1}$ | $\cdots$ | $y _{5}$ | $(y _{6} + y _{12}) / 2$ | $y _{7}$ | $\cdots$ | $y _{11}$ | $y _{13}$ | $\cdots$ | $y _{17}$ | $(y _{18} + y _{24}) / 2$ | $y _{19}$ | $\cdots$ | $y _{23}$ | $y _{25}$ | $\cdots$ | $y _{29}$ | $(y _{30} + y _{36}) / 2$ | $y _{31}$ | $\cdots$ | $y _{35}$ |

CAD図を用いた説明

よって、 $Z_1, Z_2, Z_3$ を無視する場合、 $\vec{a _{6}}$ と $\vec{a _{12}}$ は縮退することが分かる。これは $\vec{a _{18}}$ と $\vec{a _{24}}$ 、 $\vec{a _{30}}$ と $\vec{a _{36}}$ でも同様である。
上で述べた通り、作用行列を構成する各ベクトルが縮退している場合、解が一意に定まらない可能性がある。

任意のWH駆動量ベクトル $\vec{x}$ を与えた時の鏡面の変形のzernike係数ベクトル $\vec{b}$ を計算するための近似作用行列 $A$ を以下の式を満たすように定義する。

$$
\vec{b} = A \vec{x}
$$
$$
\begin{bmatrix}
    b _1 \\
    \vdots \\
    b _{11}
\end{bmatrix} =
\begin{bmatrix}
    A _{1, 1} & \cdots & A _{1, 36} \\
    \vdots & \ddots & \vdots \\
    A _{11, 1} & \cdots & A _{11, 36}
\end{bmatrix}
\begin{bmatrix}
    x _1 \\
    \vdots \\
    x _{36}
\end{bmatrix}
$$

この時 $D, \vec{d}, \vec{x}$ の間での議論と同様に考えると、

$$
A =
\begin{bmatrix}
    A _{1, 1} & \cdots & A _{1, 36} \\
    \vdots & \ddots & \vdots \\
    A _{11, 1} & \cdots & A _{11, 36}
\end{bmatrix}
= \left [ \vec{a _1}, \cdots, \vec{a _{36}} \right ]
$$

である。ここで $\vec{a _1}, \cdots, \vec{a _{36}}$ は $\vec{d _1}, \cdots, \vec{d _{36}}$ をzernike fittingしたものである。よって、任意のWH駆動量ベクトル $\vec{x}$ に対して変形量ベクトルのzernike fitting結果 $\vec{b}$ を計算するには、やはり係数 $\alpha _1, \cdots, \alpha _{36}$ とFEMの結果があればよい。

ここで、既知の $\vec{b}$ を最小二乗的に再現するような $\vec{x_r}$ は

$$
\| A \vec{x} - \vec{b} \|
$$

を最小化する $\vec{x_r}$ である。よって、 $A$ の疑似逆行列 $A^+ = A^T (A A^T)^{-1}$ を用いれば

$$
\vec{x_r} = A^+ \vec{b} \\
$$

となる。ここで求めた $\vec{x_r}$ を用いて

$$
\vec{f_r} = D \vec{x_r}
$$

を計算する。鈴木M論ではこの疑似逆行列を用いる方法で最適WH駆動量ベクトルを求めた。

# 制約付き

実際の能動支持機構では、WH駆動量ベクトル $\vec{x}$ の各要素 $x_n$ の範囲には制約がある。
よって、 $x _{min} \leq x_n \leq x _{max}$ で $\| A \vec{x} - \vec{b} \|$ を最小化する $\vec{x_r}$ を探す必要がある。
この制約付き最小化を解析的に解くことはできないため、本研究では数値解を得るために `scipy.optimize.lsq_linear` を用いた

# 行列置き場

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

| $x _{1}$ | $\cdots$ | $x _{5}$ | $x _{6}$ | $x _{7}$ | $\cdots$ | $x _{11}$ | $x _{12}$ | $\cdots$ | $x _{16}$ | $x _{17}$ | $x _{18}$ | $\cdots$ | $x _{22}$ | $x _{23}$ | $\cdots$ | $x _{27}$ | $x _{28}$ | $x _{29}$ | $\cdots$ | $x _{33}$ |
| - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - |
| $y _{1}$ | $\cdots$ | $y _{5}$ | $(y _{6} + y _{12}) / 2$ | $y _{7}$ | $\cdots$ | $y _{11}$ | $y _{13}$ | $\cdots$ | $y _{17}$ | $(y _{18} + y _{24}) / 2$ | $y _{19}$ | $\cdots$ | $y _{23}$ | $y _{25}$ | $\cdots$ | $y _{29}$ | $(y _{30} + y _{36}) / 2$ | $y _{31}$ | $\cdots$ | $y _{35}$ |

| $x _{1}$ | $x _{2}$ | $x _{3}$ | $x _{4}$ | $x _{5}$ | $x _{6}$ | $x _{7}$ | $x _{8}$ | $x _{9}$ | $x _{10}$ | $x _{11}$ |
| - | - | - | - | - | - | - | - | - | - | - |
| $y _{1}$ | $y _{2}$ | $y _{3}$ | $y _{4}$ | $y _{5}$ | $y _{6} + y _{12}$ | $y _{7}$ | $y _{8}$ | $y _{9}$ | $y _{10}$ | $y _{11}$ |