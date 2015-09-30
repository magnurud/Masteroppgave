set wrap
"au FileType tex map <F12> :! gnome-terminal -e "pdflatex ./main.tex" <enter>
au FileType tex map <F12> :! bibtex main.aux && gnome-terminal -e "pdflatex ./main.tex" <enter>
