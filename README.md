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


2. Figures generated using TikZ `graphs`  and `graphdrawing` library, see: 
   1. `graphset-images` The images in this directory have to be
   compiled using `lualatex`. See `generate_graphset_images.sh` in
   this subdirectory: `lualatex -shell-escape  -interaction=nonstopmode dcs-lp.tex`

	2. `metabolic-images` This contains the images of metabolic
	reaction networks.
	
	3. Convex hull in TikZ generated using the function:

            \newcommand{\convexpath}[2]{
            [   
                create hullnodes/.code={
                    \global\edef\namelist{#1}
                    \foreach [count=\counter] \nodename in \namelist {
                        \global\edef\numberofnodes{\counter}
                        \node at (\nodename) [draw=none,name=hullnode\counter] {};
                    }
                    \node at (hullnode\numberofnodes) [name=hullnode0,draw=none] {};
                    \pgfmathtruncatemacro\lastnumber{\numberofnodes+1}
                    \node at (hullnode1) [name=hullnode\lastnumber,draw=none] {};
                },
                create hullnodes
            ]
            ($(hullnode1)!#2!-90:(hullnode0)$)
            \foreach [
                evaluate=\currentnode as \previousnode using \currentnode-1,
                evaluate=\currentnode as \nextnode using \currentnode+1
                ] \currentnode in {1,...,\numberofnodes} {
            -- ($(hullnode\currentnode)!#2!-90:(hullnode\previousnode)$)
              let \p1 = ($(hullnode\currentnode)!#2!-90:(hullnode\previousnode) - (hullnode\currentnode)$),
                \n1 = {atan2(\y1,\x1)}, 
                \p2 = ($(hullnode\currentnode)!#2!90:(hullnode\nextnode) - (hullnode\currentnode)$),
                \n2 = {atan2(\y2,\x2)},
                \n{delta} = {-Mod(\n1-\n2,360)}
              in 
                {arc [start angle=\n1, delta angle=\n{delta}, radius=#2]}
            }
            -- cycle
            }
            % usage: \draw[violet,thick] \convexpath{1,4,3,2}{0.4cm};
	
