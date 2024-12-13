# How to use MENTAT in *Mathematica* 14.1

MENTAT's design philosophy prioritises understandability and ease-of-use. When writing expressions that you want to be treated according to Fock-state calculation rules, it usually isn't necessary to use any special code or make any API calls. Just write down the expression as you might with pen and paper, and MENTAT will often be able to interpret and simplify the results algebraically without any additional help.

Here are a few tips and pointers:

## Defining states
Defining a Ket state or Bra state is simple. Kets and Bras can be defined explicitly via `Ket[]` or `Bra[]`. The Fock mode numbers are passed to `Ket[]` / `Bra[]` as a list of integers or symbols `{i,j,k,...}`. 
Alternatively, Kets and Bras can be written symbolically via the shortcut `esc` `ket` `esc` or `esc` `bra` `esc` and the mode numbers written directly inside them as a list of integers or symbols separated by commas: $\ket{i,j,k,...}$. `esc` represents Mathematica's Alias Delimiter.

Density matrices are written similarly. To define a density state, simply define the Ket and Bra components and 'join' them with the SmallCircle ($\circ$) operator `esc` `sc` `esc`: $\ket{i,j,k...}\circ\bra{x,y,z...}$.

Interactions with scalars are all done automatically - Kets, Bras, and density matrices can be multiplied or added just like any other function.

## Vector products
Mathematica's command `CenterDot` ($\cdot$) is used for vector multiplication operations within MENTAT. `CenterDot` can be written easily using the shortcut `esc` `.` `esc`. `CenterDot` covers things like vector inner and outer products:

$\alpha\bra{0,1} \cdot \beta\ket{0,1} = \alpha\beta$

$\alpha\ket{0,1} \cdot \beta\bra{0,1} = \alpha\beta\ket{0,1}\circ\bra{0,1}.$

## Operators and conjugate transposition
MENTAT also provides creation and annihilation operator functionality. Annihilation operators are defined as any symbol with a hat and an integer subscript: use `ctrl` `&` and `^` to place the hats, and `ctrl` `_` to place the subscript. Creation operators are the same as annihilation operators, but they have a dagger - add this using `shift` `^` and `esc` `dg` `esc`. Operators interact with states via `CenterDot`:

$\hat{a}_1 \cdot \ket{n} = \sqrt{n}\ket{n-1}$

$\hat{a}_1^\dagger \cdot \ket{n} = \sqrt{n+1}\ket{n+1}.$

Adding the dagger as a superscript is actually a general part of MENTAT syntax, specifying the conjugate transpose. It can be used in other contexts, such as being applied to states:

$(\alpha\ket{1})^\dagger = \alpha^* \bra{1}$.

## Tensor product
The last basic bit of functionality MENTAT offers is the tensor product operation CircleTimes ($\otimes$). The shortcut for `CircleTimes` is `esc` `c*` `esc`. Tensor products are used to create joint quantum states comprising two states from separate Hilbert spaces:

$(\alpha\ket{i}\circ\bra{k})\otimes(\beta\ket{j}\circ\bra{l}) = \alpha \beta \ket{i,j}\circ\bra{k,l}.$

## Further functionality
MENTAT offers also offers higher-order functionality useful for Fock-state calculations, such as partial traces, beamsplitters, etc. Documentation for these functions can be accessed directly via [MENTAT.wl](https://github.com/nicholaszaunders/MENTAT/blob/main/MENTAT.wl).
