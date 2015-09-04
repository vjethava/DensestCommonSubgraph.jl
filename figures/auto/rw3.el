(TeX-add-style-hook
 "rw3"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("preview" "active" "tightpage" "border=10pt")))
   (TeX-run-style-hooks
    "latex2e"
    "minimal"
    "minimal10"
    "tikz"
    "amsmath"
    "preview")))

