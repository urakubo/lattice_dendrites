import sys, os



bound_color = [192,192,192]
mito_color  = [255,255,152]
er_color    = [179,255,179]
surfaces = {1: [bound_verts, bound_faces, bound_color],\
			2: [mito_verts, mito_faces, mito_color],\
			3: [er_verts, er_faces, er_color]}
save_uniem_annotator(annot_folder, xyzpitch, (vol_dend+vol_mito+vol_er*2).astype('uint16'), surfaces)

