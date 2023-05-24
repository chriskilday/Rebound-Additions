/* output.h*/
#include "rebound.h"

void clean_files(char*);
void write_header(char*, int, double, double, double);
void write_sim(struct reb_simulation* const, char*, double);
void write_particles(struct reb_simulation* const, char*, double);
void write_particle(struct reb_particle, FILE*, struct reb_simulation* const, struct reb_particle);
void write_tree(struct reb_simulation* const, char*, int);
struct reb_particle get_com(struct reb_simulation* const);
void move_to_com(struct reb_simulation* const);
