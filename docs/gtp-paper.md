**Combinatorial Analysis of Diagonal, Box, and Greater-Than Polynomials as Packing Functions**  
**Summary and Structured Extraction** (Jose Torres-Jimenez et al., 2015)

### 1. Greater-Than Polynomials (GTP)

**Definition**  
GTP is a bijection that maps strictly increasing vectors \( V = (v_1, v_2, \dots, v_m) \) where \( 0 \leq v_1 < v_2 < \dots < v_m \) to non-negative integers \( \alpha \in \mathbb{N} \):

\[
\alpha = \sum_{i=1}^{m} \binom{v_i}{i}
\]

**Direct Method** (Algorithm 1 - DirectGTP)

```pseudocode
α := 0
for i := 1 to m do
    α := α + Binomial(v_i, i)
return α
```

**Binomial(k, r)** (Algorithm 2)

```pseudocode
if k < r then return 0
b := 1
for i := 1 to r do
    b := (b * (k - i + 1)) / i
return b
```

**Inverse Method** (Algorithm 3 - InverseGTP)

```pseudocode
f := m!
for i := m downto 1 do
    a_i := floor( (i! * α)^{1/i} )   // integer i-th root
    v_i := a_i
    Δ := Binomial(v_i, i)
    while NextBinomial(Δ, v_i, i) ≤ α do
        Δ := NextBinomial(Δ, v_i, i)
        v_i := v_i + 1
    α := α - Δ
    f := f / i
return V
```

**NextBinomial** (Algorithm 4) – incremental computation of binomial coefficients.

**Complexity**  
- Direct: \( O(n^2 \cdot m^2) \)  
- Inverse: \( O(n^3 \cdot m) \) (improved version presented in the paper)

**Applications**  
- Generation of all \( t \)-subsets of \( k \) columns (Algorithm 5 - AllSubsets).  
- Random generation of \( t \)-subsets for metaheuristics in Covering Array construction.

---

### 2. Diagonal Polynomials (DP)

**Definition**  
DP maps vectors \( W = (w_1, w_2, \dots, w_m) \) (non-negative integers, no strict order) to \( \mathbb{N} \) by grouping them into diagonals based on the sum of components. For \( m=2 \):

\[
\alpha = w_1 + \binom{1 + w_1 + w_2}{2}
\]

**General Direct Method** (Algorithm 6)

```pseudocode
s := 0, α := 0
for i := 1 to m do
    s := s + w_i
    α := α + Binomial(i-1 + s, i)
return α
```

**Complexity**  
- Direct: \( O(n^2 \cdot m^2) \)  
- Inverse: \( O(n^3 \cdot m^2) \)

**Application with GTP**  
Generation of random solutions for Linear Diophantine Equations with Unit coefficients (LDEU), e.g., \( x_1 + x_2 + \dots + x_r = n \).

---

### 3. Box Polynomials (BP)

**Definition**  
BP uses “box shells” defined by the maximum value \( \beta \) in the vector and the index \( l \) of its first occurrence. Vectors are enumerated according to these classes.

**Direct Method** (Algorithm 7)

The contribution of the class is:

\[
\beta^m + (\beta + 1)^m - \beta^{l-1} (\beta + 1)^{m-l+1}
\]

The specific vector contribution is computed as a mixed-radix number (first \( l-1 \) digits base \( \beta \), remaining digits base \( \beta+1 \)).

**Inverse Method** (Algorithm 8)

1. Compute \( \beta = \lfloor \alpha^{1/m} \rfloor \).  
2. Determine the class index \( l \).  
3. Decode the remaining components using the appropriate bases.

**Complexity** (Most Efficient)

| Polynomial | Direct Method       | Inverse Method      |
|------------|---------------------|---------------------|
| DP         | \( O(n^2 \cdot m^2) \) | \( O(n^3 \cdot m^2) \) |
| GTP        | \( O(n^2 \cdot m^2) \) | \( O(n^3 \cdot m) \)   |
| **BP**     | **\( O(n^2 \cdot m) \)** | **\( O(n^3 \cdot m) \)** |

**Spatial complexity** for all methods: \( O(m) \).

---

### Additional Resources in the Document

- **Table 1**: GTP enumeration example (strictly upper triangular).  
- **Table 2**: DP enumeration example (diagonal structure).  
- **Table 4**: BP enumeration example (box structure).  
- **Table 5**: First 27 vectors of dimension 3 with BP (\( \beta \), \( l \), vector).  
- **Table 6**: Summary of computational complexities.  
- **Figure 1**: Example of a Covering Array CA(9; 2, 4, 3).

**Recommendation**: Box Polynomials (BP) are concluded to be the most efficient packing functions among the three.

