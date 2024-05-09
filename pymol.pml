fetch 6ZP2
remove chain B chain C

set ray_trace_fog,0
set ray_shadows,0
unset depth_cue
bg_color white
set antialias,2
set hash_max, 300
set ray_trace_mode,  1
set ambient, 0.9


create terminals, 6ZP2 and not resid 331-528
create cterminal, terminals and resid 529-1136
create nterminal, resid 0-330 and terminals
create RBD, resid 330-529


color red, resid 18 resid 19 resid 19 resid 26 resid 67 resid 95 resid 138 resid 142 resid 190 resid 213 resid 339 resid 346 resid 371 resid 375 resid 376 resid 405 resid 408 resid 417 resid 417 resid 440 resid 446 resid 452 resid 477 resid 478 resid 484 resid 484 resid 493 resid 501 resid 505 resid 547 resid 614 resid 655 resid 679 resid 681 resid 681 resid 764 resid 796 resid 831 resid 856 resid 950 resid 954 resid 969 resid 981 resid 1027 resid 1153 resid 1176

show surface
ray 1450,1000
png raytraced_surface2.png, dpi=300
