(TeX-add-style-hook "singlerand"
 (lambda ()
    (TeX-add-symbols
     '("code" 1))
    (TeX-run-style-hooks
     "latex2e"
     "art10"
     "article")))

