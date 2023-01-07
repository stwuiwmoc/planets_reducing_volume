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
- [最適モーター駆動量ベクトルの導出](#最適モーター駆動量ベクトルの導出)
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

上式において、既知の $\vec{f}$ を最小二乗的に再現（reproduction）するような最適モーター駆動量ベクトル $\vec{y_r}$ を考える。すなわち

$$
\| D \vec{y} - \vec{f}\|
$$

を最小化するような $\vec{y_r}$ を考える。
このような $\vec{y_r}$ は、 正規方程式の解としてただ一つに定まることが知られている。
$\vec{y_r}$ は $D$ の疑似逆行列 $D^+ = D^T (D D^T)^{-1}$ を用いて

$$
\vec{y_r} = D^+ \vec{f}
$$

ただし、 $D$ の要素である $\vec{d _1}, \vec{d _2}, \cdots, \vec{d _{36}}$ の中に同一な要素が含まれる（縮退している）場合、解が発散し一意に定まらない可能性がある。よって、解 $\vec{y_r}$ が一意に定まるためには $\vec{d _n}$ が互いに独立でなければならない（縮退してはならない）ことに注意が必要である。
この $\vec{y_r}$ を用いれば、既知の $\vec{f}$ を最小二乗的に再現した結果である $\vec{f_r}$ は

$$
\vec{f_r} = D \vec{y_r}
$$

となる。以上の方法を用いることで、主鏡と支持機構を望遠鏡の架台に設置した後で生じた任意の鏡面形状誤差を、能動支持機構を用いて補正することができる。

# 近似作用行列の導入

本研究の能動支持機構は、一般的な使用方法である「支持機構設置後に生じた鏡面形状誤差の補正」だけでなく、「主鏡製造時に除去し切れなかった鏡面形状誤差の補正」にも用いることを想定している。後者の使用方法では、鏡面形状誤差のうちpiston, tip, tilt成分については能動支持機構を望遠鏡架台に設置する際に調整することが可能である。その場合、作用行列 $D$ とその疑似逆行列 $D^+$ を用いる方法は以下に述べる理由から補正効率が悪い。

figureの通り、 $\vec{d _n}$ のRMSのほとんどはpiston, tip, tilt成分である。このため、 $y _r$ の算出においても、これら3成分の補正に駆動量のリソースの多くが用いられている。しかし、後者の使用方法であればpiston, tip, tilt成分は架台への設置時に調整可能なため、これら3成分の調整のために駆動量のリソースを割くのは無駄が多い。

このため、piston, tip, tilt成分を無視した、より効率的な $\vec{y_r}$ を計算するためにzernike多項式へのfittingを用いたzernike近似作用行列を導入する。

$\vec{d_n}$ をZernike多項式 $(j \leqq 11)$ に最小二乗フィットして得られたzernike係数ベクトルを $\vec{a_n} ^T \equiv (a _{1, n} , \cdots, a _{11, n})$ とする。
同様に、 $\vec{f}$ に対するzernike係数ベクトルを $\vec{b} ^T \equiv (b _1, \cdots, b _{11})$ とする。

後者の使用方法では、piston, tip, tilt成分に相当する $Z_1, Z_2, Z_3$ は考慮する必要がない。ここで、 $\vec{a _{6}}$ と $\vec{a _{12}}$ を比較すると、 $Z_4$ 以降では絶対値がほぼ同じで、符号が逆転しているのみであることが分かる。具体的に hogehoge % の違いしかない。このような一致が生じる理由は、WTとWHの構造によって説明できる。

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

$$
A^T =
\begin{bmatrix}
    A _{4, 1} & \cdots & A _{11, 1} \\
    \vdots & & \vdots \\
    A _{4, 5} & \cdots & A _{11, 5} \\
    A _{4, 6} & \cdots & A _{11, 6} \\
    A _{4, 7} & \cdots & A _{11, 7} \\
    \vdots & & \vdots \\
    A _{4, 11} & \cdots & A _{11, 11} \\
    A _{4, 12} & \cdots & A _{11, 12} \\
    \vdots & & \vdots \\
    A _{4, 16} & \cdots & A _{11, 16} \\
    A _{4, 17} & \cdots & A _{11, 17} \\
    A _{4, 18} & \cdots & A _{11, 18} \\
    \vdots & & \vdots \\
    A _{4, 22} & \cdots & A _{11, 22} \\
    A _{4, 23} & \cdots & A _{11, 23} \\
    \vdots & & \vdots \\
    A _{4, 27} & \cdots & A _{11, 27} \\
    A _{4, 28} & \cdots & A _{11, 28} \\
    A _{4, 29} & \cdots & A _{11, 29} \\
    \vdots & & \vdots \\
A _{4, 33} & \cdots & A _{11, 33} \\
\end{bmatrix} =
\begin{bmatrix}
    a _{4, 1} & \cdots & a _{11, 1} \\
    \vdots & & \vdots \\
    a _{4, 5} & \cdots & a _{11, 5} \\
    a _{4, 6} - a _{4, 12} & \cdots & a _{11, 6} - a _{11, 12} \\
    a _{4, 7} & \cdots & a _{11, 7} \\
    \vdots & & \vdots \\
    a _{4, 11} & \cdots & a _{11, 11} \\
    a _{4, 13} & \cdots & a _{11, 13} \\
    \vdots & & \vdots \\
    a _{4, 17} & \cdots & a _{11, 17} \\
    a _{4, 18} - a _{4, 24} & \cdots & a _{11, 18} - a _{11, 24} \\
    a _{4, 19} & \cdots & a _{11, 19} \\
    \vdots & & \vdots \\
    a _{4, 23} & \cdots & a _{11, 23} \\
    a _{4, 25} & \cdots & a _{11, 25} \\
    \vdots & & \vdots \\
    a _{4, 29} & \cdots & a _{11, 29} \\
    a _{4, 30} - a _{4, 36} & \cdots & a _{11, 30} - a _{11, 36} \\
    a _{4, 31} & \cdots & a _{11, 31} \\
    \vdots & & \vdots \\
    a _{4, 35} & \cdots & a _{11, 35} \\
\end{bmatrix} =
\begin{bmatrix}
    \overrightarrow{a _{4 \sim 11, 1}} \, ^T \\
    : \\
    \overrightarrow{a _{4 \sim 11, 5}} \, ^T \\
    \overrightarrow{a _{4 \sim 11, 6}} \, ^T - \overrightarrow{a _{4 \sim 11, 12}} \, ^T \\
    \overrightarrow{a _{4 \sim 11, 7}} \, ^T \\
    : \\
    \overrightarrow{a _{4 \sim 11, 11}} \, ^T \\
    \overrightarrow{a _{4 \sim 11, 13}} \, ^T \\
    : \\
    \overrightarrow{a _{4 \sim 11, 17}} \, ^T \\
    \overrightarrow{a _{4 \sim 11, 18}} \, ^T - \overrightarrow{a _{4 \sim 11, 24}} \, ^T \\
    \overrightarrow{a _{4 \sim 11, 19}} \, ^T \\
    : \\
    \overrightarrow{a _{4 \sim 11, 23}} \, ^T \\
    \overrightarrow{a _{4 \sim 11, 25}} \, ^T \\
    : \\
    \overrightarrow{a _{4 \sim 11, 29}} \, ^T \\
    \overrightarrow{a _{4 \sim 11, 30}} \, ^T - \overrightarrow{a _{4 \sim 11, 36}} \, ^T \\
    \overrightarrow{a _{4 \sim 11, 31}} \, ^T \\
    : \\
    \overrightarrow{a _{4 \sim 11, 35}} \, ^T \\
\end{bmatrix}
$$