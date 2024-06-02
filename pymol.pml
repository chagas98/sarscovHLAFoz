#fetch 6ZP2
fetch 7A94
remove chain B chain C

set ray_trace_fog,0
set ray_shadows,0
unset depth_cue
bg_color white
set antialias,2
set hash_max, 300
set ray_trace_mode,  1
set ambient, 1
set spec_reflect, 0


create terminals, 7A94 and not resid 331-528 and not chain D
create cterminal, terminals and resid 529-1136
create nterminal, resid 0-330 and terminals
create RBD, 7A94 and resid 330-529 and not chain D
create ace2, 7A94 and chain D
create junction, cterminal and resid 680-690

color deepblue, nterminal 
color limon, RBD
color skyblue, cterminal

as cartoon, nterminal or RBD or cterminal

remove terminals
remove 6ZP2

#color red, (resid 3 resid 18 resid 19 resid 19 resid 26 resid 67 resid 95 resid 138 resid 142 resid 190 resid 213 resid 339 resid 346 resid 346 resid 371 resid 375 resid 376 resid 405 resid 408 resid 417 resid 417 resid 440 resid 444 resid 444 resid 446 resid 452 resid 452 resid 460 resid 477 resid 478 resid 484 resid 484 resid 486 resid 493 resid 498 resid 501 resid 505 resid 547 resid 614 resid 655 resid 658 resid 677 resid 679 resid 681 resid 681 resid 764 resid 796 resid 831 resid 856 resid 950 resid 954 resid 969 resid 981 resid 1027 resid 1153 resid 1176) and (nterminal or RBD or cterminal)

#show surface
#set transparency, 0.7
#ray 1450,1000
png raytraced_surface2.png, dpi=300
