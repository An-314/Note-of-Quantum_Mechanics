#import "@preview/scripst:1.1.1": *

#show: scripst.with(
  template: "book",
  title: "量子力学",
  author: (
    "Anzreww",
  ),
  time: "甲辰春夏于清华园",
  contents: true,
  content-depth: 3,
)

#pagebreak(weak: true)

#include "chap1.typ"

#pagebreak(weak: true)

#include "chap2.typ"

#pagebreak(weak: true)

#include "chap3.typ"

#pagebreak(weak: true)

#include "chap4.typ"

#pagebreak(weak: true)
