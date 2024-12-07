/*----------------------------------------*
 * Accessible surface area calculation
 * library
 *----------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include "libasa.h"
#include "radii.h"

char err_string[50];
asa_atom_rec_p a_pdb=(asa_atom_rec_p)NULL;

/*----------------------------------------*/

char *asa_version()
{
	err_string[0] = 0;
    return(ASA_VERSION);
}

/*----------------------------------------*/

int asa_create_atom_list(asa_atom_rec_p *atom_list, int num_atoms)
{
	err_string[0] = 0;
    *atom_list = (asa_atom_rec_p) malloc(num_atoms * sizeof(asa_atom_rec));
    if (*atom_list == NULL) {
    	strcpy(err_string, "asa_lib: malloc() failed.");
    	return(0);
    }
    return(1);
}

/*----------------------------------------*/

void asa_destroy_atom_list(asa_atom_rec_p atom_list)
{
	err_string[0] = 0;
    free(atom_list);
}

/*----------------------------------------*/

int asa_calculate(asa_atom_rec_p atom_list, int num_atoms,
		double slice_width, double solvent_radius)
{
	err_string[0] = 0;
    a_pdb = atom_list;
   
    if (!surfinit(num_atoms, solvent_radius))
    	return 0;
    
    surfarea(num_atoms, slice_width, solvent_radius);
   
    surfcleanup();
    return 1;
}

/*----------------------------------------*/

char *asa_error()
{
    return(err_string);
}

/*----------------------------------------*/

int asa_assign_radii(asa_atom_rec_p atom_list, int num_atoms,
		radius_type_p sup_radii, int num_sup_radii)
{
    int atom, radius, num_radii;

	err_string[0] = 0;
    num_radii = sizeof(radius_table) / sizeof(struct radius_type);
    for (atom = 0; atom < num_atoms; atom++) {
        for (radius = 0; radius < num_radii; radius++) {
	        if (!strncmp(atom_list[atom].type, radius_table[radius].type, 4)) {
               atom_list[atom].radius = (double) radius_table[radius].size;
               goto BreakOut;
            }
        }
    	if (sup_radii != NULL) {
        	for (radius = 0; radius < num_sup_radii; radius++) {
	       		if (!strncmp(atom_list[atom].type, sup_radii[radius].type, 4)) {
            		atom_list[atom].radius = (double) sup_radii[radius].size;
            		goto BreakOut;
        		}
        	}
        }
        strcpy(err_string, "asa_lib: Cannot assign radius of ");
        strncat(err_string, atom_list[atom].type, 4);
        strcat(err_string, ".");
        return 0;
BreakOut:
        ;
    }
    return 1;
}

/*----------------------------------------*/
