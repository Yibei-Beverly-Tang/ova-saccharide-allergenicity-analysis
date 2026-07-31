# PyMOL visualization for the public-data OVA structure analysis
# Run from the repository root: pymol scripts/visualize_1ova_sites.pml
reinitialize
load structures/1OVA.pdb, ova_1ova
hide everything
show cartoon, ova_1ova and chain A
color gray70, ova_1ova and chain A
set cartoon_transparency, 0.12
show sticks, ova_1ova and chain A and resn NAG
color tv_orange, ova_1ova and chain A and resn NAG
select glycation_sites, ova_1ova and chain A and resi 97+105+216+271+328+387
show sticks, glycation_sites
color marine, glycation_sites
select n_glycosylation_site, ova_1ova and chain A and resi 298
show sticks, n_glycosylation_site
color firebrick, n_glycosylation_site
select annotated_sites, ova_1ova and chain A and resi 298+97+105+216+271+328+387
label annotated_sites and name CA, "%s%s" % (resn, resi)
show surface, ova_1ova and chain A
set transparency, 0.78, ova_1ova and chain A
set label_size, 16
set label_color, black
orient ova_1ova and chain A
bg_color white
ray 1600, 1200
png outputs/1ova_annotated_sites.png, dpi=200
