#let project(
  title: "",
  authors: (),
  time: "",
  abstract: none,
  keywords: (),
  body
) = {
  let song = ("Linux Libertine", "SimSun")
  let hei = ("Linux Libertine", "SIMHEI")
  let kai = ("Linux Libertine", "KaiTi",)
  let xbsong = "FZXiaoBiaoSong-B05"
  let code = "Consolas"
  // Moidfy the following to change the font.
  let title-font = hei
  let author-font = kai
  let body-font = song
  let heading-font = hei
  let caption-font = kai
  let header-font = kai
  let strong-font = hei
  let emph-font = kai
  let raw-font = code
  
  set document(author: authors, title: title)
  set heading(numbering: "1.1")
  set text(font: body-font, lang: "zh", region: "cn")
  set bibliography(style: "gb-7714-2015-numeric")

  show heading: it => box(width: 100%)[
    #v(0.25em)
    #set text(font: heading-font)
    #if it.numbering != none { counter(heading).display() }
    #h(0.25em)
    #it.body
]

  show heading.where(
    level: 1
  ): it => box(width: 100%)[
    #v(0.5em)
    #set align(center)
    #set heading(numbering: "一")
    #it
    #v(0.75em)
  ]

  // Title
  align(center)[
    #v(20em)
    #block(text(font: title-font, weight: 700, 2.3em, title))
    #v(5em)
  ]

  // Author information.
  pad(
    top: 0.5em,
    bottom: 0.5em,
    x: 2em,
    grid(
      columns: (1fr,) * calc.min(3, authors.len()),
      gutter: 1em,
      ..authors.map(author => align(center, text(font: author-font, author))),
    ),
  )

  // Time
  if time != "" [#align(center)[
    #v(5em)
    #set text(1em)
    #time
    ]]

  pagebreak(weak: true)

  // Table of contents
  set page(numbering: "I", number-align: center, header: align(left)[
    #set text(font: header-font)
    #title
  ])
  counter(page).update(1)

  set par(first-line-indent: 2em,leading: 1em)
  show outline.entry.where(
      level: 1
    ): it => {
      v(0.3em)
      h(-2.0em)
      set text(15pt)
      strong(it)
      }
  outline(indent: auto, depth: 2)


  // Main body
  set par(first-line-indent: 2em,leading: 1.1em)
  set enum(indent: 2em)
  set figure(gap: 0.5cm)

  show figure: it => [
    #v(12pt)
    #set text(font: caption-font)
    #it
    // #par()[#text(size: 0.0em)[#h(0.0em)]]
    #v(12pt)
  ]

  show image: it => [
    #it
    #par()[#text(size: 0.0em)[#h(0.0em)]]
  ]

  show table: it => [
    #set text(font: body-font)
    #it
  ]
  show strong: set text(font: strong-font)
  show emph: set text(font: emph-font)
  show ref: set text(red)
  show raw.where(block: true): block.with(
    width: 100%,
    fill: luma(240),
    inset: 10pt,
  )
  show raw.where(block: true): set par(leading: 0.7em)

  show raw: set text(font: (raw-font, "simsun"), size: 10pt)
  show link: underline
  show link: set text(blue)

  

  set page(numbering: "1", number-align: center, header: align(left)[
    #set text(font: header-font)
    #title
  ])

  counter(page).update(1)

  if abstract != none [
    #v(2pt)
    #h(2em) *摘要：* #abstract

    #if keywords!= () [
      *关键字：* #keywords.join("；")
    ]
    #v(2pt)
  ]

  body
}
#let problem-counter = counter("problem")
#problem-counter.step()

#let problem(body) = {
  problem-counter.step()
  set enum(numbering: "(1)")
  block(
    fill: rgb(241, 241, 255),
    inset: 8pt,
    radius: 2pt,
    width: 100%,
  )[*题目 #problem-counter.display().* #h(0.75em) #body]
}

#let solution(body) = {
  set enum(numbering: "(1)")
  block(
    inset: 8pt,
    width: 100%
  )[*解答.* #h(0.75em) #body]
}