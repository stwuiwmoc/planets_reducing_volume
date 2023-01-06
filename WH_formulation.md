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
- [作用行列の導入](#作用行列の導入)
- [最適WH駆動量ベクトルの導出](#最適wh駆動量ベクトルの導出)
- [近似作用行列の導入](#近似作用行列の導入)
- [制約付き](#制約付き)

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

# 作用行列の導入

鏡面上に均等に配置されたm個の点を考え、各点での鏡面の高さ方向の変形量 $f _1, f _2, \cdots, f _{m}$  [m] を要素とする変形量ベクトル $\vec{f}$ を $\vec{f} ^T \equiv (f _1 , \cdots, f _{m})$ と定義する。

また、WHのM01, M02, ... , M36 の駆動量 $x _1, x _2, \cdots, x _{36}$ [mm] を要素とするWH駆動量ベクトル $\vec{x}$ を $\vec{x} ^T \equiv (x _1 , \cdots , x _{36})$ と定義する。

任意のWH駆動量ベクトル $\vec{x}$ を与えた時の鏡面の変形量ベクトル $\vec{f}$ を計算するためのWH駆動量変換行列 $D$ （以後、作用行列と呼ぶ）を以下の式を満たすように定義する。

$$
\vec{f} = D \vec{x}
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
    x _1 \\
    \vdots \\
    x _{36}
\end{bmatrix}
$$

よって、$D$ の各要素の値が事前に分かっていれば、任意の駆動量 $\vec{x}$ に対して変形量 $\vec{f}$ を得ることができる。

ここで、 $\vec{x} ^T = (1, 0, \cdots, 0)$ のときの変形量 $\vec{d _1}$ を考える。 $\vec{d _1}$ は $D$ を用いて

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

と計算できるため、 $\vec{x} ^T = (1, 0, \cdots, 0)$ に対する変形量は $\vec{d _1} ^T = (D _{1, 1}, \cdots, D _{m, 1})$ である。同様に、 $\vec{x} ^T = (0, 1, \cdots, 0)$ に対する変形量は $\vec{d _2} ^T = (D _{1, 2}, \cdots, D _{m, 2})$ であり、これを $\vec{d _{36}}$ まで繰り返すと

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

ここで、 $\vec{d _1}$の物理的意味は、M01の駆動量 $x _1 = 1$ [mm]、他の駆動量 $x _2, \cdots, x _{36} = 0$ [mm] の時の変形量であり、以下に述べる方法で決定することができる。

まず、WHが駆動していないときの各点の高さ $\vec{h _0} ^T = (h _{1, 0}, \cdots, h _{m, 0})$ [m] を有限要素解析（FEM）で求める。

次に、act1に0.05 [Nm] のトルクを与えたときの各点の高さ $\vec{h _1} ^T = (h _{1, 1}, \cdots, h _{m, 1})$ [m] をFEMで求める。

このとき、act1 に0.05 [Nm] のトルクを与えたときの各点の変形量は $\vec{h _1} - \vec{h _0}$ [m] である。

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

以上より、係数 $\alpha _1, \cdots, \alpha _{36}$ とFEMの結果があれば、任意のWH駆動量ベクトル $\vec{x}$ に対して変形量ベクトル $\vec{f}$ を計算できる。

# 最適WH駆動量ベクトルの導出

$$
\vec{f} = D \vec{x}
$$

上式において、既知の $\vec{f}$ を最小二乗的に再現（reproduction）するような最適WH駆動量ベクトル $\vec{x_r}$ を考える。すなわち

$$
\| D \vec{x} - \vec{f}\|
$$

を最小化するような $\vec{x_r}$ を考える。このような $\vec{x_r}$ は正規方程式の解としてただ一つに定まることが知られており、$D$ の疑似逆行列 $D^+ = D^T (D D^T)^{-1}$ を用いて

$$
\vec{x_r} = D^+ \vec{f}
$$

である。この $\vec{x_r}$ を用いれば、既知の $\vec{f}$ を最小二乗的に再現した結果である $\vec{f_r}$ は

$$
\vec{f_r} = D \vec{x_r}
$$

となる。

# 近似作用行列の導入

上記の導出法では、

- m 個の点で均等に重みづけして最適化を行うため、形状を正確に反映できない
- $m \times 36$ の行列計算をする必要があり、計算量が多い

という問題がある。このため、より適切に $\vec{x_r}$ を計算するために近似作用行列を導入する。

$\vec{d_n}$ の zernike fittingで得られたzernike係数ベクトルを $\vec{a_n} ^T \equiv (a _{1, n} , \cdots, a _{11, n})$ とする。
同様に、 $\vec{f}$ の zernike fittingで得られたzernike係数ベクトルを $\vec{b} ^T \equiv (b _1, \cdots, b _{11})$ とする。

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
