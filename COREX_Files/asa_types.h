/*----------------------------------------*
 * This is the libasa header file
 *----------------------------------------*/

typedef struct asa_atom_rec_struct {
    double x, y, z, radius, area;
    char type[4];
} asa_atom_rec, *asa_atom_rec_p;

typedef struct radius_type_struct {
	char type[5];
	float size;
} radius_type, *radius_type_p;

char *asa_error();

int asa_calculate(asa_atom_rec_p atom_list, int num_atoms,
		double slice_width, double solvent_radius);

int asa_create_atom_list(asa_atom_rec_p *atom_list, int num_atoms);

void asa_destroy_atom_list(asa_atom_rec_p atom_list);

int asa_assign_radii(asa_atom_rec_p atom_list, int num_atoms,
                     radius_type_p sup_radii, int num_sup_radii);

char *asa_version();

/*----------------------------------------*/
