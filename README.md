# ECML Poster

ETH Poster style copied from http://www.pletscher.org/

- `PosterVersion2.pdf` contains fonts in included images fixed manually using Adobe Illustrator.
- `PosterEcml15.tex` is the main file

Compilation command: `latexmk -pdf -pvc -xelatex PosterEcml15.tex`

## Special formatting

1. Uses `Garamond` and  `Cambria Math` for correct printing

		\usefonttheme{professionalfonts}
		\usepackage{xltxtra}
		\usepackage{unicode-math}
		\setmathfont{Cambria Math}
		\usepackage{fontspec}
		
		% Font specification here. 
		\setmainfont{Garamond}
		\setsansfont{Garamond}
		\setmonofont{Consolas}
		\setmathrm{Cambria Math}
		\setmathsf{Cambria Math}
		\setmathtt{Cambria Math}


2. Figures generated using TikZ `graphs` library, see: 
   1. `graphset-images` The images in this directory have to be
   compiled using `lualatex`. See `generate_graphset_images.sh` in
   this subdirectory: `lualatex -shell-escape  -interaction=nonstopmode dcs-lp.tex`

	2. `metabolic-images` This contains the images of metabolic
	reaction networks.
	
