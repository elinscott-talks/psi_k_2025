#import "@preview/mannot:0.3.0": *
#import "touying/lib.typ": *
#import "@preview/pinit:0.1.4": *
#import "@preview/xarrow:0.3.0": xarrow
#import "@preview/cetz:0.4.1"
#import "psi-slides-0.6.1.typ": *
#import "@preview/grayness:0.3.0": image-grayscale
#import "@preview/algorithmic:1.0.3"
#import algorithmic: algorithm

// color-scheme can be navy-red, blue-green, or pink-yellow
// #let s = psi-slides.register(aspect-ratio: "16-9", color-scheme: "pink-yellow")
#show: psi-theme.with(aspect-ratio: "16-9",
                      color-scheme: "pink-yellow",
                             config-info(
                                title: [Koopmans functionals],
                                subtitle: [How satisfying piecewise linearity can yield reliable band structures],
                                author: [Edward Linscott],
                                date: datetime(year: 2025, month: 8, day: 27),
                                location: [Psi-k],
                                references: [references.bib],
                             ))

#set footnote.entry(clearance: 0em)
#show bibliography: set text(0.6em)

#let primary = rgb("#dc005a")
#let secondary = rgb("#f0f500")

#let blcite(reference) = {
  text(fill: white, cite(reference))
}

#let delayedmark(start, content, tag: none, color: red, mark: mark, color-before: black, alternatives: none) = {
   let entries = (mark(content, tag: tag, color: color-before),)*(start - 1) + (mark(content, tag: tag, color: color),)
   alternatives(repeat-last: true, ..entries)
}

#let methods-with-marks(self) = {
  let (uncover, only, alternatives) = utils.methods(self)
  let dm = delayedmark.with(alternatives: alternatives)
  (uncover, only, alternatives, dm, dm.with(mark: markhl, color-before: white))
}

#title-slide()

== Predicting spectral properties

Spectral properties are fundamental to understanding materials...

#align(center, 
grid(columns: 3, column-gutter:  1em,
cetz.canvas({
  import cetz.draw: *

  // grid((0,-5), (8,5), stroke: gray + .5pt)

  // Valence
  rect((-1, -1), (1, 1), stroke: none, fill: secondary, alpha: 0.5)
  content((1.75, 1), [$E_F$], align: left)
  circle((0, 0), radius: 0.2, fill: white, stroke: (dash: "dashed", paint: primary))

  // Vacuum
  circle((0, 4), radius: 0.2, fill: primary, stroke: none)
  line((-1, 3.5), (1, 3.5), stroke: (dash: "dashed", paint: primary))
  content((1.75, 3.5), [$E_"vac"$], align: left)

  // Arrow
  arc((0,0), start: -30deg, stop: 30deg, radius: 4, mark: (end: ">", fill: black))
  
  let photon(amplitude: 1, phases: 2, scale: 8, samples: 1000, angle: 0, start-x: 0, start-y: 0, ..args) = {
    line(..(for x in range(0, samples + 1) {
      let x = x / samples
      // A Gaussian envelope with sigma = 1/4 and mean = 1/2 and height = amplitude
      let envelope = amplitude * calc.exp(-calc.pow(((x - 0.5) / (0.25)), 2))

      let phase = (2 * phases * calc.pi) * x

      // Rotate the output by angle
      let xval = x * scale
      let yval = calc.sin(phase) * envelope

      let rotated-x = xval * calc.cos(angle) - yval * calc.sin(angle)
      let rotated-y = xval * calc.sin(angle) + yval * calc.cos(angle)
      ((start-x + rotated-x, start-y + rotated-y),)
    }), ..args)

    let subdivs = 8
    for phase in range(0, phases) {
      let x = phase / phases
      for div in range(1, subdivs + 1) {
        let p = 2 * calc.pi * (div / subdivs)
        let y = calc.sin(p) * amplitude
        let x = x * scale + div / subdivs * scale / phases
      }
    }
  }
  photon(amplitude: 0.8, phases: 9, start-x: -0.25, start-y: 0.25, scale: 3, fill: none, angle: 2.5, mark: (start: ">", fill: black))
}),
  image("figures/arpes.png", height: 45%),
  image("figures/arpes_puppin.png", height: 45%),
))

#blcite(<delaTorre2021>)#blcite(<Puppin2020>)

#pause ... but how can we routinely compute them? #pause
- GW: accurate, expensive, often ill-behaved, diagrammatic
- DFT: plagued by intrinsic errors
#pause

#box[#move(dy: 0.1em, image("figures/lightbulb.png", height: 1em))] Koopmans functionals: overcome limitations of DFT $arrow.r$ a functional that can accurately predict single-particle excitations
#v(-1em)

= The failures of DFT

== Total energy differences vs. eigenvalues

#align(horizon,
grid(align: horizon, columns: 2, column-gutter: 1em,
  [
We all know that DFT underestimates the band gap. But why? #pause

The exact Green's function has poles that correspond to total energy differences

$
  ε_i = cases(E(N) - E_i (N-1) & "if" i in "occ", E_i (N+1) - E(N) & "if" i in "emp")
$

#pause

but DFT does #emph[not]
],[
  #image(width: 20em, "figures/fig_en_curve_gradients_zoom.svg")
]
))

#focus-slide()[Core idea: impose this condition on DFT]

== Imposing generalised piecewise linearity
#align(horizon,
grid(align: horizon, columns: (1fr, auto), column-gutter: 1em,
  [Formally, every orbital $i$ should have an eigenenergy
  $
    epsilon_i^"Koopmans" = ⟨
      phi_i mid(|)hat(H)mid(|)phi_i
    ⟩ = frac(dif E, dif f_i)
  $
  that is
  - independent of $f_i$
  - equal to $Delta E$ of explicit electron addition/removal
],[
  #image(width: 20em, "figures/fig_en_curve_gradients_zoom.svg")
]
))

== Imposing generalised piecewise linearity
#slide(repeat: 3, self => [
  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

#align(horizon,
grid(align: horizon, columns: (1fr, auto), column-gutter: 1em,
[
  $
  E^"KI" &[rho, {rho_i}] =
  E^"DFT" [rho]
  \ & +
  sum_i (
    - delayedmarkhl(#2, integral_0^f_i lr(angle.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i angle.r) dif f, tag: #<remove_nonlin>, color: primary)
  \ &
    + delayedmarkhl(#3, f_i integral_0^1 lr(angle.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i angle.r) dif f, tag: #<restore_linear>, color: primary)
  )
$
// Bakes the total energy differences $E^"DFT" [rho^(f_i arrow.r 1)] - E^"DFT" [rho^(f_i arrow.r 0)]$ into the functional
#uncover("2-")[#annot(<remove_nonlin>, pos: bottom)[#align(center, [removes dependence on $f_i$])]]
#uncover("3-")[#annot(<restore_linear>, pos: bottom)[#align(center, [restores linear dependence on $f_i$])]]

],[
  #image(width: 20em, "figures/fig_en_curve_gradients_zoom.svg")
]))

])
== Comparison with DFT+_U_ (and BLOR)
#slide()[

#set table(
  fill: (x, y) =>
    if calc.rem(y, 2) == 0 { silver } else { silver.lighten(75%) },
    inset: 0.5em,
)
#show table.cell.where(x: 0): set text(style: "italic")

#show table.header: {
  set text(fill: primary, weight: "semibold")
}

#table(align: horizon, columns: (1fr, 2fr, 2fr), stroke: none,
    table.header([], [*DFT+_U_*], [*Koopmans*]),
    table.hline(),
   [seeks to correct...],
   uncover("2-")[erroneous curvature in total energies w.r.t. $N$],
   uncover("4-")[erroneous curvature in total energies w.r.t. $f_i forall i$],
   [in practice...],
   uncover("3-")[corrects curvature in total energies w.r.t. local manifold (BLOR does so more faithfully)],
   uncover("5-")[removes dependence of $epsilon_i$ on $f_i$ and guarantees $epsilon_i = E_i (N plus.minus 1) - E(N)$],
   [correction applied to...],
   [],
   [],
   [orbitals defined by...],
   [],
   [],
   [parametrised by...],
   [],
   [],
)
]
== Electronic screening via parameters
#slide(repeat: 4, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)
  $
    E^"KI" [{rho_i}] = &
    E^"DFT" [rho]
    +
    sum_i (
      - integral_0^f_i lr(angle.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i angle.r) dif f
      + f_i integral_0^1 lr(angle.l phi_i mid(|) hat(h)^"DFT" (f) mid(|) phi_i angle.r) dif f
    )
    #pause
    \ = & E^"DFT" [rho]
    + sum_i {
      - (E^"DFT" [rho] - delayedmark(#3, E^"DFT" [rho^(f_i arrow.r 0)], tag: #<ENm1_hard>, color: primary))
      + f_i (delayedmark(#3, E^"DFT" [rho^(f_i arrow.r 1)], tag: #<ENp1_hard>, color: primary) - delayedmark(#3, E^"DFT" [rho^(f_i arrow.r 0)], tag: #<ENm1b_hard>, color: primary))
    }
    // uncover("5-", 
    // \ arrow.r E^"uKI" [{rho_i}] approx & 
    // E^"DFT" [rho]
    // \ & +
    // sum_i {
    //   - (E^"DFT" [rho] - delayedmark(#3, E^"DFT" [rho - rho_i], tag: #<ENm1>, color: primary))
    //   + f_i (delayedmark(#3, E^"DFT" [rho - rho_i + n_i], tag: #<ENp1>, color: primary) - delayedmark(#3, E^"DFT" [rho - rho_i], tag: #<ENm1b>, color: primary))
    // }
    // )
  $

  #pause
  #annot(<ENm1_hard>, pos: bottom)[#align(center, [cannot evaluate \ directly])]
  #annot(<ENp1_hard>, pos: bottom)[#align(center, [cannot evaluate \ directly])]
  #annot(<ENm1b_hard>, pos: bottom)[#align(center, [cannot evaluate \ directly])]
  #pause

  // Instead use a frozen-orbital picture:
  
  // $
  //  rho^(f_i arrow.r f)(bold(r)) approx rho(bold(r)) + (f - f_i) |phi^N_i (bold(r))|^2
  // $
  // 
  // very easy to evaluate -- but not at all accurate! Correct this _post hoc_ via a screening parameter i.e.
  // 
  // $
  //   E[rho^(f_i arrow.r f)] approx alpha_i E[rho + (f - f_i) |phi^N_i (bold(r))|^2]
  // $
])

#slide[
#align(center + horizon, 
  image("figures/fig_pwl_DFT.svg", height: 100%)
)
]
#slide[
#align(center + horizon, 
  image("figures/fig_pwl_uKI.svg", height: 100%)
)
]
#slide[
#align(center + horizon, 
  image("figures/fig_pwl_alphaKI.svg", height: 100%)
)
]

#slide(repeat: 5, self => [

  #let (uncover, only, alternatives, delayedmark, delayedmarkhl) = methods-with-marks(self)

$
  E^"KI"_bold(alpha) [rho, {rho_i}] approx & 
  E^"DFT" [rho]
  \ & +
  sum_i delayedmark(#3, alpha_i, tag: #<alpha>, color: primary) {
    - (E^"DFT" [rho] - delayedmark(#2, E^"DFT" [rho - rho_i], tag: #<ENm1>, color: primary))
    + f_i (delayedmark(#2, E^"DFT" [rho - rho_i + n_i], tag: #<ENp1>, color: primary) - delayedmark(#2, E^"DFT" [rho - rho_i], tag: #<ENm1b>, color: primary))
  }
$

#pause
#annot(<ENm1>, pos: bottom)[uses frozen orbitals]
#annot(<ENp1>, pos: bottom)[uses frozen orbitals]
#annot(<ENm1b>, pos: bottom)[uses frozen orbitals]
#pause
#annot(<alpha>, pos: bottom)[screening parameter]
#pause
which is easy to evaluate _e.g._
$ H^"KI"_(i j) = angle.l phi_j|hat(h)^"DFT" + alpha_i hat(v)_i^"KI"|phi_i angle.r #h(2cm) hat(v)^"KI"_i = - E_"Hxc" [rho - n_i] + E_"Hxc" [rho] - integral v_"Hxc" (bold(r)', [rho]) n_i d bold(r)' $

#pause
Screening parameters _not_ a fitting parameter!

])

// == Screening
// 
// #grid(columns: (1fr, 2fr), 
// [
// #align(center + horizon, 
//   image("figures/fig_pwl.png", width: 100%)
// )
// ],
// [
// 
// #pause
// Construct $alpha_i$ from explicit $Delta$SCF calculations@Nguyen2018@DeGennaro2022a
// 
// $
//   alpha_i = alpha_i^0 (Delta E_i - lambda_(i i)(0)) / (lambda_(i i)(alpha^0) - lambda_(i i)(0)) "where" lambda_(i i)(alpha) = angle.l phi_i|hat(h)^"DFT" + alpha hat(v)_i^"KI"|phi_i angle.r $
// 
// #pause
// Recast via linear response@Colonna2018:
// 
// $
//   alpha_i = (angle.l n_i mid(|) epsilon^(-1) f_"Hxc" mid(|) n_i angle.r) / (angle.l n_i mid(|) f_"Hxc" mid(|) n_i angle.r)
// $
// 
// which can be efficiently computed via DFPT@Colonna2022
// ],
// )

== Orbital-density dependence
#slide()[
The potential is orbital-density-dependent!
#v(-0.5em)
  $ v^"KI"_(i in"occ") = - E_"Hxc" [rho - n_i] + E_"Hxc" [rho] - integral v_"Hxc" (bold(r)', [rho]) n_i d bold(r)' $

#pause

- loss of unitary invariance@Nguyen2018
#v(-1em)
#align(center,
  grid(columns: (auto, auto), column-gutter: 1em,
  image("figures/fig_nguyen_variational_orbital.png", width: 10em),
  image("figures/fig_nguyen_canonical_orbital.png", width: 10em),
  [two variational orbitals],
  [a canonical orbital],
  )
) #pause
- we can use MLWFs@Marzari2012 #pause
- we know $hat(H)|phi_i angle.r$ but not $hat(H)$ #pause
- a natural generalisation of DFT towards spectral functional theory@Ferretti2014
]

== To summarise...
$
  E^"KI"_bold(alpha) [rho, {rho_i}] =
  E^"DFT" [rho] +
  sum_i alpha_i { &
    - (E^"DFT" [rho] - E^"DFT" [rho - rho_i])
  \ &
    + f_i (E^"DFT" [rho - rho_i + n_i] - E^"DFT" [rho - rho_i])
  }
$

- an orbital-by-orbital correction to DFT
- screening parameters
- orbital-density-dependence
- total energy at integer occupations unchanged!

== Comparison with DFT+_U_ (and BLOR)
#slide()[

#set table(
  fill: (x, y) =>
    if calc.rem(y, 2) == 0 { silver } else { silver.lighten(75%) },
    inset: 0.5em,
)

#show table.cell.where(y: 0): set text(weight: "bold", fill: primary)

#show table.cell.where(x: 0): set text(style: "italic")

#show table.cell: it => {
  set text(size: 0.9em)
  it
}

#table(align: horizon, columns: (1fr, 2fr, 2fr), stroke: none,
    table.header([], [DFT+_U_], [Koopmans]),
    table.hline(),
   [seeks to correct...],
   [erroneous global curvature in total energies w.r.t. $N$],
   [erroneous global curvature in total energies w.r.t. #uncover("6-")[*canonical*] orbital occupancies],
   [in practice...],
   [corrects curvature in total energies w.r.t. local manifold (BLOR does so more faithfully)],
   [removes dependence of $epsilon_i$ on #uncover("7-")[*variational*] orbital occupations and guarantees $epsilon_i = E_i (N plus.minus 1) - E(N)$],
   [correction applied to...],
   uncover("2-")[selected subspaces (e.g. _3d_ orbitals)],
   uncover("4-")[the entire system],
   [orbitals defined by...],
   uncover("3-")[Hubbard projectors (atom-centred, frozen, incomplete)],
   uncover("5-")[variational (localised) orbitals],
   [parametrised by...],
   uncover("8-")[${U_I}$], //, defined w.r.t. charge-neutral excitations if using LR],
   uncover("9-")[${alpha_i}$], // defined w.r.t. charged excitations]

)
  
]

= Results

== Molecular systems

=== Ionisation potentials@Colonna2019
#align(center + horizon,
image("figures/colonna_2019_gw100_ip.jpeg", width: 100%)
)

=== UV photoemission spectra@Nguyen2015
#align(center + horizon,
image("figures/fig_nguyen_prl_spectra_pink.png", width: 100%)
)


== Extended systems
#slide[
=== Prototypical semiconductors and insulators @Nguyen2018

#show table.cell: it => {
  if it.x == 3 or it.x == 4 {
    set text(fill: primary, weight: "semibold")
    it
  } else {
    it
  }
}

#grid(align: center + horizon, columns: 2, column-gutter: 1em,
image("figures/scatter_plot.png", height: 80%),
table(columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr), inset: 0.5em, stroke: none,
table.header([], [PBE], [G#sub[0]W#sub[0]], [KI], [KIPZ], [QSGW̃]),
table.hline(),
[$E_"gap"$], [2.54], [0.56], [0.27], [0.22], [0.18],
[IP], [1.09], [0.39], [0.19], [0.21], [0.49]
))
  
]

#slide[
=== ZnO @Colonna2022
#v(-1em)
#align(center + horizon,
grid(align: center + horizon, columns: 3, column-gutter: 1em,
image("figures/ZnO_lda_cropped.png", height: 50%),
image("figures/ZnO_hse_cropped_noaxis.png", height: 50%),
image("figures/ZnO_ki_cropped_noaxis.png", height: 50%),
))
#show table.cell: it => {
  set text(size: 0.8em)
  if it.x == 5 {
    set text(fill: primary, weight: "semibold")
    it
  } else {
    it
  }
}
#table(columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1.5fr), align: center, inset: 0.5em, stroke: none,
table.header([], [LDA ], [HSE ], [GW#sub[0] ], [scGW̃ ], [KI ], [exp ]),
table.hline(),
[$E_"gap"$], [0.79], [2.79], [3.0], [3.2], [3.68], [3.60],
[$angle.l epsilon_d angle.r$], [-5.1], [-6.1], [-6.4], [-6.7], [-6.93], [-7.5 to -8.81 ],
[$Delta$], [4.15], [], [], [], [4.99], [5.3]
)
  
]

== Model systems
=== Hooke's atom@Schubert2023

#align(center + horizon, 
  image("figures/schubert_vxc_only.jpeg", height: 70%)
)

= Caveats

== Limitations

- only valid for systems with $E_"gap"$ > 0 #pause
- empty state localisation in the bulk limit #pause
- can break crystal point group symmetry

== Resonance with other efforts

- Wannier transition state method of Anisimov and Kozhevnikov@Anisimov2005
- Optimally-tuned range-separated hybrid functionals of Kronik, Pasquarello, and others@Kronik2012@Wing2021
- Ensemble DFT of Kraisler and Kronik@Kraisler2013
- Koopmans-Wannier method of Wang and co-workers@Ma2016
- Dielectric-dependent hybrid functionals of Galli and co-workers@Skone2016a
- Scaling corrections of Yang and co-workers@Li2018

= Extensions
= Non-collinear spin
== Non-collinear spin

$ rho_i (bold(r)) pause arrow.r bold(rho)_i (bold(r)) = (rho_i (bold(r)), m_i^x (bold(r)), m_i^y (bold(r)), m_i^z (bold(r))) $

#pause

e.g. for the corrective potential

$ v_i^"qKI" = - 1 / 2 integral dif bold(r) dif bold(r)' rho_i (bold(r)) f_"Hxc" (bold(r), bold(r)') rho_i (bold(r)') + (1 - f_i) integral d bold(r)' f_"Hxc" (bold(r), bold(r)') rho_i (bold(r)') $

#pause

#align(center, sym.arrow.b)

$ v_i^"qKI" = - 1 / 2 integral dif bold(r) dif bold(r)' bold(rho)_i (bold(r)) bb(F)_"Hxc" (bold(r), bold(r)') bold(rho)_i (bold(r)') sigma_0 + (1 - f_i) sum_alpha integral d bold(r)' [bb(F)_"Hxc" (bold(r), bold(r)') bold(rho)_i (bold(r)')]_alpha sigma_alpha $


#blcite(<Marrazzo2024>)

#pagebreak()

CsPbBr#sub[3] #blcite(<Marrazzo2024>)
#v(-2em)
#align(center + horizon,
image("figures/marrazzo_CsPbBr3_bands.svg", height: 50%)
)
#table(align: center, columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1.5fr), inset: 0.5em, stroke: none,
table.header([], [LDA ], [HSE ], [G#sub[0]W#sub[0] ], [scGW̃ ], [*KI*], [exp ]),
table.hline(),
[*with SOC*], [0.18], [0.78], [0.94], [1.53], [*1.78*], [1.85],
[without SOC], [1.40], [2.09], [2.56], [3.15], [3.12], [],
)

= Optical spectra
== Optical spectra

Solve the BSE, using Koopmans eigenvalues in lieu of GW

#pause

#v(-1em)
#align(center + horizon,
grid(columns: 2,
image("figures/silicon_bse_spectra.png", height: 50%),
image("figures/silicon_bse_excitons.png", height: 50%)
))

#v(-1em)

#show table.cell: it => {
  set text(size: 0.8em)
  it
}
#table(align: center + horizon, columns: (auto, 1fr, 1fr, 1fr, 1fr), inset: 0.5em, stroke: none,
table.header([silicon], [indirect gap ], [direct gap ], [first excitonic peak ], [excitonic binding energy ]),
table.hline(),
[*qKI+BSE*], [1.12], [3.31], [3.42], [0.09], 
[G#sub[0]W#sub[0]+BSE], [1.17], [3.25], [3.34], [0.09],
)

= Computational cost and scaling
== Computational cost and scaling


The vast majority of the computational cost: determining screening parameters

$
  alpha_i = (angle.l n_i|epsilon^(-1) f_"Hxc"|n_i angle.r) / (angle.l n_i|f_"Hxc"|n_i angle.r)
$

#pause

- a local measure of screening of electronic interactions #pause
- one screening parameter per orbital #pause
- must be computed #emph[ab initio] via... #pause
  - $Delta$SCF@Nguyen2018@DeGennaro2022a: embarrassingly parallel steps which each cost $cal(O)(N_"SC"^3) tilde cal(O)(N_bold(k)^3 N^3)$ #pause
  - DFPT@Colonna2018@Colonna2022: $cal(O)(N_bold(k)^2 N^3)$

#pagebreak()
#align(center + horizon,
image("figures/timings/benchmark.svg", width: 80%)
)

== Machine-learned electronic screening
#slide[
  #grid(
    columns: (1fr, 1fr),
    align: center + horizon,
    gutter: 1em,
    image(
      "figures/convergence_key.png",
      height: 5%,
    ) +  v(-1em) +
    image(
      "figures/convergence_fig.png",
      height: 55%,
    ),
    image("figures/speedup.png", height: 60%),

    [*accurate* to within $cal("O")$(10 meV) _cf._ typical band gap accuracy of $cal("O")$(100 meV)],
    [*speedup* of $cal("O")$(10) to $cal("O")$(100)],
  )

  #blcite(<Schubert2024>)
]

// == 
// #image("figures/supercell_workflow.png", width: 100%)
// 
// #image("figures/primitive_workflow.png", width: 65.5%)

#focus-slide()[
#align(center, image(width: 80%, "media/logos/koopmans_white_on_transparent.svg"))
]

== `koopmans`
#matrix-slide(alignment: horizon)[
  #image("figures/website_cropped.png")
][
  
  - automated workflows
  - `Quantum ESPRESSO` backend
  - easy installation
  - python API
  
  See `koopmans-functionals.org`
]

==
#align(center + horizon,
image("figures/supercell_workflow.png", width: 100%)
)

#matrix-slide(alignment: horizon, columns: (3fr, 2fr))[
  #image("figures/black_box_filled_square.png")
][
 
  Our goal:
  + accurate
  + robust
  + minimal input
  + fast

]

== Automated Wannierisation
#slide()[
  Koopmans functionals rely heavily on Wannier functions...
  - to initialise the minmising orbitals, _or_
  - in place of the minimising orbitals entirely

#pause

#grid(
  columns: (2fr, 2fr, 3fr),
  align: center + horizon,
  gutter: 1em,
  image("figures/proj_disentanglement_fig1a.png", height: 45%),
  image("figures/new_projs.png", height: 45%),
  image("figures/target_manifolds_fig1b.png", height: 45%),

  text("projectability-based disentanglement") + cite(<Qiao2023>),
  text("use PAOs found in pseudopotentials"),
  text("parallel transport to separate manifolds") + cite(<Qiao2023a>),
)
]

== 
#blcite(<Huber2020>)
#v(-2em)
#align(center,
  [
  #grid(columns: 3, align: horizon, column-gutter: 0.5em,
    image("media/logos/koopmans_grey_on_transparent.svg", height: 3em),
    image("figures/handshake.png", height: 2em, alt: "handshake"),
    image("media/logos/aiida.svg", height: 3em)
  )
  #pause `$ koopmans run tio2.json` #pause $arrow.r$ `$ koopmans run --engine=aiida tio2.json`
  ]
)

remote compute, parallel step execution, provenance-tracking, (requires configuration, WIP...)

#pause
#align(center, 
  image("figures/aiida-speed-up.svg", width: 70%)
)

== 
#slide()[
#show table.cell: it => {
  if it.x == 5 {
    set text(fill: primary, weight: "semibold")
    it
  } else {
    it
  }
}
#grid(columns: 3, align: horizon + center,
  [
  #set text(size: 0.45em)
  #raw(read("scripts/gaas_auto.json"), block: true, lang: "json")
  ],
  [
    #set text(size: 3em)
    #sym.arrow.r
  ],
  [
    #image("figures/Unfold_And_Interpolate_bandstructure.png", height: 60%)
    #table(columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr), inset: 0.5em, stroke: none,
    table.header([], [LDA], [HSE], [GW#sub[0]], [scGW̃ ], [KI], [exp]),
    table.hline(),
    [$E_"gap"$], [0.26], [1.28], [1.55], [1.62], [1.54], [1.55],
    [$angle.l epsilon_d angle.r$], [-14.9], [-15.6], [-17.3], [-17.6], [-17.9], [-18.9],
    [$Delta$], [12.8], [13.9], [], [], [12.7], [13.1]
    )
  ]
)
]

= Summary
== Summary
#grid(
  columns: (1fr, 2fr),
  gutter: 1em,
  image("figures/black_box_filled_square.png", width: 100%),
  text[
    Koopmans functionals...
    - impose generalised piecewise linearity condition to DFT
    - give band structures with comparable accuracy to state-of-the-art GW
    - can be used in place of GW in BSE calculation of excitons, for systems with strong SOC, ...
    - are increasingly black-box
  ],
)

== Open questions

#pause
- why does correcting _local_ charged excitations correct the description of delocalized excitations? #pause
- is there a good metric for selecting variational orbitals (_i.e._ the subspace with respect to which we enforce piecewise linearity)? #pause
- are off-diagonal corrections appropriate? What form should they take? #pause
- how to extend to metallic systems? #pause
- can we provide a formal basis for the Koopmans correction? #pause
  - GKS
  - spectral functional theory@Ferretti2014
  - ensemble DFT
  - RDMFT

== Want to find out more?
#slide()[
#show grid.cell: it => {
  let weight = "regular"
  if it.y == 1 {
    weight = "bold"
  }
  if it.x < 3 {
    set text(size: 0.8em, fill: silver.darken(50%), weight: weight)
    it
  } else {
    set text(size: 0.8em, weight: weight)
    it
  }
}
#let nm_data = read("media/mugshots/nicola_marzari.jpeg", encoding: none)
#let ms_data = read("media/mugshots/marija_stojkovic.jpg", encoding: none)
#let nc_data = read("media/mugshots/nicola_colonna.png", encoding: none)
#align(center + horizon, 
grid(columns: (1fr, 1fr, 1fr, 1fr, 1.2fr), align: center + top, inset: 0.5em,
  image-grayscale(nm_data, height: 50%),
  image-grayscale(ms_data, height: 50%),
  image-grayscale(nc_data, height: 50%),
  image("media/mugshots/junfeng_qiao.jpeg", height: 50%),
  image("media/mugshots/aleksandr_poliukhin.jpg", height: 50%),
  [Nicola Marzari], [Marija Stojkovic], [Nicola Colonna], [Junfeng Qiao], [Aleksandr Poliukhin],
  [Monday], [Monday], [Tuesday], [Poster B4.16 today!], [Thu 1000 Room C],
  [spectral theories], [band alignment for photocatalysis], [non-collinear spin], [automated Wannierisation], [electron-phonon]
) + [... or go to `koopmans-functionals.org`]
)
  
]


== Acknowledgements
#align(center + horizon, 
grid(columns: 9, column-gutter: 0.5em, align: center, row-gutter: 0.5em,
  image("media/mugshots/nicola_colonna.png", height: 40%),
  image("media/mugshots/miki_bonacci.jpg", height: 40%),
  image("media/mugshots/aleksandr_poliukhin.jpg", height: 40%),
  image("media/mugshots/marija_stojkovic.jpg", height: 40%),
  image("media/mugshots/giovanni_cistaro.jpeg", height: 40%),
  image("media/mugshots/julian_geiger.jpg", height: 40%),
  image("media/mugshots/junfeng_qiao.jpeg", height: 40%),
  image("media/mugshots/yannick_schubert.jpg", height: 40%),
  image("media/mugshots/nicola_marzari.jpeg", height: 40%),
  [Nicola Colonna], [Miki Bonacci], [Aleksandr Poliukhin], [Marija Stojkovic], [Giovanni Cistaro], [Julian Geiger], [Junfeng Qiao], [Yannick Schubert], [Nicola Marzari]
)
)

#align(
  center,
  grid(
    columns: 2,
    align: horizon + center,
    gutter: 2em,
    image("media/logos/SNF_logo_standard_web_color_pos_e.svg", height: 20%),
    image("media/logos/marvel_color_on_transparent.png", height: 20%),
  ),
)


#focus-slide()[#align(center, text(size: 2em, [Thank you!]) + linebreak() + text(size: 0.5em, style: "italic", [these slides are available at #h(0.2em) #box[#move(dy: 0.1em, image("media/logos/github-mark-white.svg", height: 1em))] `elinscott-talks`]))]


#show: appendix

#focus-slide()[#align(center, text(size: 2em, [spare slides]))]

== Frozen orbital approximation

  #v(-5em)
  #align(center + horizon, 
  grid(align: center + horizon, columns: 3, column-gutter: 2cm, row-gutter: 1cm,
  cetz.canvas({
    import cetz.draw: *
    content((1.25, 1.5), [$rho$])
    circle((0, 0), radius: 1, fill: primary, stroke: none)
    circle((2.5, 0), radius: 1, fill: primary, stroke: none)

  }),
  cetz.canvas({
    import cetz.draw: *

    content((9, 1.5), [$rho^(f_1 arrow.r 0)$])
    arc((10.75, 0), start: 0deg, stop: 360deg, radius: (1.5, 1), fill: primary, stroke: none)
    circle((8, 0), radius: 1, fill: none, stroke: (thickness: 2pt, paint: primary))
    circle((8, 0), radius: 1, fill: none, stroke: (dash: "dashed", thickness: 2pt, paint: white))
    // content((8, -1.5), [$f_1 = 0$])
  }),
  cetz.canvas({
    import cetz.draw: *

    content((17.25, 1.5), [$rho - |psi^N_1(r)|^2$])
    circle((16, 0), radius: 1, fill: none, stroke: (dash: "dashed", thickness: 2pt, paint: primary))
    circle((18.5, 0), radius: 1, fill: primary, stroke: none)
  }),
  [2-electron solution],
  [what we'd like to evaluate],
  [what we can quickly evaluate]

  ))

#matrix-slide(columns: (3fr, 2fr))[
#align(center + horizon,
  {only("1")[#image("figures/alpha_calc/fig_alpha_calc_step_0.png", height: 80%)]
  only("2")[#image("figures/alpha_calc/fig_alpha_calc_step_1.png", height: 80%)]
  only("3")[#image("figures/alpha_calc/fig_alpha_calc_step_2.png", height: 80%)]
  only("4-5")[#image("figures/alpha_calc/fig_alpha_calc_step_3.png", height: 80%)]
  only("6-7")[#image("figures/alpha_calc/fig_alpha_calc_step_4.png", height: 80%)]
  }
)
][
#only("7")[$ alpha_i = alpha_i^0 (Delta E_i - lambda_(i i)(0)) / (lambda_(i i)(alpha^0) - lambda_(i i)(0)) $
$ lambda_(i i)(alpha) = angle.l phi_i|hat(h)^"DFT" + alpha hat(v)_i^"KI"|phi_i angle.r $]
]

== Issues with extended systems

#align(center + horizon, 
  image("figures/fig_nguyen_scaling.png", width: 60%)
)

Two options: #pause _1._ use a more advanced functional#pause, or _2._ stay in the "safe" region
#blcite(<Nguyen2018>)

#slide()[
#set text(size: 0.8em)
#raw(read("scripts/gaas.json"), block: true, lang: "json")
]


== Machine-learned electronic screening

#pagebreak()

#slide[
  #align(
    center,
    grid(
      columns: 5,
      align: horizon,
      gutter: 1em,
      image("figures/orbital.emp.00191_cropped.png", height: 30%),
      $stretch(->)^("power spectrum decomposition")$,
      $vec(delim: "[", x_0, x_1, x_2, dots.v)$,
      $stretch(->)^("ridge regression")$,
      $alpha_i$,
    ),
  )

  $
    c^i_(n l m, k) & = integral dif bold(r) g_(n l) (r) Y_(l m)(theta,phi) n^i (
      bold(r) - bold(R)^i
    )
  $


  $
    p^i_(n_1 n_2 l,k_1 k_2) = pi sqrt(8 / (2l+1)) sum_m c_(n_1 l m,k_1)^(i *) c_(n_2 l m,k_2)^i
  $

  #blcite(<Schubert2024>)
]

#pagebreak()

#slide[
  #align(
    center,
    grid(
      columns: 2,
      align: horizon + center,
      gutter: 1em,
      image("figures/water.png", height: 70%),
      image("figures/CsSnI3_disordered.png", height: 70%),

      "water", "CsSnI" + sub("3"),
    ),
  )
  #blcite(<Schubert2024>)
]

The use-case

   #grid(columns: 8, column-gutter: 0.3em, row-gutter: 0.3em,
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        image("figures/CsSnI3_disordered.png", width: 100%),
        grid.cell(align: center + horizon, [...]),
        grid.cell(inset: 0.4em, align: center, fill: primary, colspan: 3, text(fill: white, "train", size: 1em, weight: "bold")),
        grid.cell(inset: 0.4em, align: center, fill: secondary, colspan: 5, text("predict", size: 1em, weight: "bold")),
  )

  #pause
  N.B. not a general model


#slide[
  #grid(
    columns: (1fr, 1fr),
    align: horizon + center,
    gutter: 1em,
    image(
      "figures/water_cls_calc_vs_pred_and_hist_bottom_panel_alphas.svg",
      height: 70%,
    ),
    image(
      "figures/CsSnI3_calc_vs_pred_and_hist_bottom_panel_alphas.svg",
      height: 70%,
    ),

    "water", "CsSnI" + sub("3"),
  )
  #blcite(<Schubert2024>)
]

#slide[
  #grid(
    columns: (1fr, 1fr),
    align: center + horizon,
    gutter: 1em,
    image(
      "figures/convergence_key.png",
      height: 5%,
    ) +  v(-1em) +
    image(
      "figures/convergence_fig.png",
      height: 55%,
    ),
    image("figures/speedup.png", height: 60%),

    [*accurate* to within $cal("O")$(10 meV) _cf._ typical band gap accuracy of $cal("O")$(100 meV)],
    [*speedup* of $cal("O")$(10) to $cal("O")$(100)],
  )

  #blcite(<Schubert2024>)
]


== Taking advantage of symmetries
To compute screening parameters via DFPT...
#algorithm(inset: 0.3em, indent: 1em, {
  import algorithmic: *
  Function("CalculateAlpha", ($n$,), {
    For($bold(q) in "BZ"$,
    {
        For($bold(k) in "BZ"$, {Comment[Linear system $A x = b$ to obtain $Delta psi_(bold(k)+bold(q),v)(bold(r))$]})
          Assign[$Delta rho^(0n)_(q)$][$sum_(bold(k)v)psi^*_(bold(k)v) (bold(r))Delta psi_(bold(k)+bold(q),v)(bold(r)) + c.c.$]
          Assign[$Pi^((r))_(0 n, bold(q))$][$angle.l Delta rho^(0 n)_(bold(q))|f_"Hxc"|rho^(0 n)_(bold(q)) angle.r$]
          Assign[$Pi^((u))_(0 n, bold(q))$][$angle.l rho^(0 n)_bold(q)|f_"Hxc"|rho^(0 n)_bold(q) angle.r$]
    })
    Return[$1 + sum_bold(q) Pi^((r))_(0 n, bold(q)) \/ sum_bold(q) Pi^((u))_(0 n, bold(q))$]
  })
})

#pagebreak()

#align(center,
  image("figures/bz-to-ibz-outer.svg", height: 80%)
)
$bold(q) in "BZ" $ $arrow.r$ $bold(q) in "IBZ"(n)$ (the symmetry of the perturbation; lower than that of the primitive cell)
#pagebreak()
#align(center,
  image("figures/bz-to-ibz-inner.svg", height: 80%)
)
$bold(k) in "BZ"$ $arrow.r$ $bold(k) in "IBZ"(bold(q))$ (can only use symmetries that leave $bold(q)$ invariant)

#align(horizon + center, image("figures/bz-to-ibz-speedup.svg", height: 100%))


== Connections with approx. self-energies

#blcite(<Ferretti2014>)#blcite(<Colonna2019>)

Orbital-density functional theory:

$ (h + alpha_i v^(K I)_i)|psi_i angle.r = lambda_i|psi_i angle.r $ $v_i^(K I)(bold(r))$ is real, local, and state-dependent #pause

cf. Green's function theory:

$ (h + Sigma_i)|psi_i angle.r = z_i|psi_i angle.r $ $Sigma_i (bold(r), bold(r)')$ is complex, non-local, and state-dependent

#slide[
Hartree-Fock self-energy in localized representation

$Sigma_x (bold(r), bold(r)') = - sum_(k sigma)^("occ") psi_(k sigma)(bold(r)) & f_H (bold(r), bold(r'))psi^*_(k sigma)(bold(r)') \
& arrow.r.double.long angle.l phi_(i sigma)|Sigma_x|phi_(j sigma') angle.r approx - angle.l phi_(i sigma)|v_H [n_(i sigma)]|phi_(i sigma)angle.r delta_(i j)delta_(sigma sigma')$

Unscreened KIPZ#sym.at Hartree ($v_"xc" arrow.r 0$; $f_"Hxc" arrow.r f_H$; $epsilon^(-1) arrow.r 1$)

$angle.l phi_(i sigma)|v^"KIPZ"_(j sigma',"xc")|phi_(j sigma') angle.r
approx {(1/2 - f_(i sigma)) angle.l n_(i sigma)|f_H|n_(i sigma) angle.r - E_H [n_(i sigma)]}
approx - angle.l phi_(i sigma)|v_H [n_(i sigma)]|phi_(i sigma)angle.r delta_(i j)delta_(sigma sigma')$

]

#slide[
Screened exchange plus Coulomb hole (COHSEX)

$ Sigma^"SEX"_"xc" (bold(s), bold(s)') = - sum_(k sigma)^"occ" psi_(k sigma)(bold(r)) psi_(k sigma)^*(bold(r)) W(bold(r), bold(r)') $

$ Sigma^"COH"_"xc" (bold(s), bold(s)') = 1/2 delta(bold(s), bold(s)'){W(bold(r), bold(r)') - f_H (bold(r), bold(r)')} $

$ arrow.r.double.long angle.l phi_(i sigma)|Sigma^"COHSEX"_"xc"|phi_(j sigma')angle.r approx {(1/2 - f_(i sigma)) angle.l n_(i sigma)|W|n_(i sigma)angle.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$

KIPZ#sym.at Hartree with RPA screening ($v_"xc" arrow.r 0$; $f_"Hxc" arrow.r f_H$; $epsilon^(-1) arrow.r "RPA"$)

$ angle.l phi_(i sigma)|v^"KIPZ"_(j sigma',"xc")|phi_(j sigma')angle.r approx{(1/2 - f_(i sigma)) angle.l n_(i sigma)|W|n_(i sigma)angle.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$
]

#slide[
  Static GWΓ#sub[xc] --- local (DFT-based) vertex corrections@Hybertsen1987@DelSole1994

  $ Sigma^(G W Gamma_"xc")_"xc"(1, 2) = i G(1, 2) W_(t-e) (1, 2) $
  
  $ W_(t-e) = (1 - f_"Hxc" chi_0)^(-1) f_H $

  $ arrow.r.double.long angle.l phi_(i sigma)|Sigma^(G W Gamma_"xc")_"xc"|phi_(j sigma')angle.r approx{(1/2 - f_(i sigma)) angle.l n_(i sigma)|W_(t-e)|n_(i sigma)angle.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$

  KIPZ#sym.at DFT ($v_"xc" arrow.r$ DFT; $f_"Hxc" arrow.r$ DFT; $epsilon^(-1) arrow.r$ DFT)

  $ angle.l phi_(i sigma)|v^"KIPZ"_(j sigma',"xc")|phi_(j sigma')angle.r approx{angle.l phi_(i sigma)|v^"DFT"_(sigma,"xc")|phi_(i sigma)angle.r + (1/2 - f_(i sigma)) angle.l n_(i sigma)|epsilon^(-1)_(t-e) f_"Hxc"|n_(i sigma)angle.r - E_H [n_(i sigma)]}delta_(i j) delta_(sigma sigma')$
]

= References
== References
#slide()[
#show bibliography: set text(0.95em)
#bibliography("references.bib", style: "nature-footnote.csl", title: none)

]