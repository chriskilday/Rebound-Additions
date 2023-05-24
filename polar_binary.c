/* polar_binary.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "tree.h"
#include "output.h"
#include "communication_mpi.h"

typedef struct reb_simulation simulation;

void heartbeat(simulation* const);
void read_initial_conditions(int, char**);
void simulation_init(simulation*);
void initial_conditions(simulation*);

// Default initial conditions
int TMAX = 1e4;
double MASS_RATIO = 1.0, A = 1.0, E = 0.75;
int N_INITIAL = 1e4;
double DISC_MASS = 1e-1;
double DISC_INNER = 2;
double DISC_OUTER = 20;
int POLAR = 0;

const double BOXN = 9;
const double BOXSIZE = 75.0 / BOXN;
const double VRAND_FRAC = 0e-1;

char out_pt[512];
char out_tree[512];

int OUTPUT_INTERVAL = 100;
double OUTPUT_START = 0;
double OUTPUT_END = 1e32;

int main(int argc, char* argv[]) {
	simulation* const s = reb_create_simulation();

	read_initial_conditions(argc, argv);
	simulation_init(s);

	if (s->mpi_id == 0) {
		initial_conditions(s);
		clean_files(out_pt);
		write_header(out_pt, N_INITIAL, MASS_RATIO, A, E);
	}

	reb_communication_mpi_distribute_particles(s);
	move_to_com(s);
	MPI_Barrier(MPI_COMM_WORLD);

	reb_integrate(s, TMAX);

	MPI_Barrier(MPI_COMM_WORLD);
	reb_mpi_finalize(s);
	reb_free_simulation(s);
}

void heartbeat(simulation* const s) {
	move_to_com(s);

	if (reb_output_check(s, TMAX / 10000.))
		reb_output_timing(s, 0);

	if (s->t >= OUTPUT_START && s->t <= OUTPUT_END
			&& reb_output_check(s, OUTPUT_INTERVAL))
		write_sim(s, out_pt, MASS_RATIO);

	if (s->t == TMAX && s->mpi_id == 0)
		printf("\n");
}

void simulation_init(simulation* s) {
	s->integrator = REB_INTEGRATOR_LEAPFROG;
	s->gravity = REB_GRAVITY_TREE;
	s->boundary = REB_BOUNDARY_OPEN;
	s->heartbeat = heartbeat;
	s->dt = 3e-2;
	s->G = 1;
	s->opening_angle2 = 1.5;
	s->softening = 0.02;

	reb_configure_box(s, BOXSIZE, BOXN,BOXN,1);
	reb_mpi_init(s);
}

void initialize_binary(simulation* s) {
	/* binary */
	struct reb_particle star1 = {0}, star2 = {0};
	star1.m = 1.0;
	star2.m = star1.m * MASS_RATIO;

	double mu = star1.m*star2.m / (star1.m + star2.m);
	double vkep = sqrt(s->G*mu/A * (2/(1+E) - 1));

	if (POLAR) { star1.z = -0.5*A*(1+E); star2.z = 0.5*A*(1+E); }
	else { star1.x = -0.5 * A*(1+E); star2.x = 0.5 * A*(1+E); }
	star1.vy = -vkep; star2.vy = vkep;

	star1.hash = reb_hash("Star1");
	star2.hash = reb_hash("Star2");

	reb_add(s, star1);
	reb_add(s, star2);

	/* unary */
/*
	struct reb_particle star = {0}; star.m = 1.0;
	star.hash = reb_hash("Star1");
	reb_add(s, star);
*/
}

double MANUAL_INNER = 0;
void initialize_disc(simulation* s) {
	double inner_a = DISC_INNER * (1+E)*A, outer_a = DISC_OUTER * (1+E)*A;
	double m = DISC_MASS/N_INITIAL;
	if (MANUAL_INNER < inner_a) MANUAL_INNER = inner_a;

	for (int i = 0; i < N_INITIAL; i++) {
		struct reb_particle pt = {0};
		double a = reb_random_powerlaw(s, inner_a, outer_a, -1.5);
		double phi = reb_random_uniform(s, 0, 2*M_PI);
		double mu = (1.0 + MASS_RATIO) + (DISC_MASS - 2*m*(1/(sqrt(inner_a)) - 1/(sqrt(MANUAL_INNER)))) * 
				(pow(a, -3./2.) - pow(MANUAL_INNER, -3./2.)) / (pow(outer_a, -3./2.) - pow(MANUAL_INNER, -3./2.));
		double vkep = sqrt(s->G*mu/a);

		pt.x = a*cos(phi); pt.y = a*sin(phi);
		pt.z = a*reb_random_normal(s, 0.001);
		pt.vx = vkep*sin(phi); pt.vy = -vkep*cos(phi); pt.vz = 0;
		pt.m = m;

		double vrand = reb_random_uniform(s, 0, VRAND_FRAC * vkep);
		double min, max, rand1 = reb_random_uniform(s,0,1), rand2 = reb_random_uniform(s,0,1);
		if (rand1 < rand2) { min = rand1; max = rand2; }
		else { min = rand2; max = rand1; }

		double zrand_frac = min, xrand_frac = max-min, yrand_frac = 1-max;
		pt.vz += zrand_frac * vrand * ((rand()%2)*2 - 1);
		pt.vx += xrand_frac * vrand * ((rand()%2)*2 - 1);
		pt.vy += yrand_frac * vrand * ((rand()%2)*2 - 1);

		char hash_str[16]; sprintf(hash_str, "%d", i);
		pt.hash = reb_hash(hash_str);

		if (a >= MANUAL_INNER) reb_add(s, pt);
	}
}

void initial_conditions(simulation* s) {
	initialize_binary(s);
	initialize_disc(s);
}

void read_initial_conditions(int argc, char* argv[]) {
	sprintf(out_pt, "%s", "out_particles");
	sprintf(out_tree, "%s", "out_tree");

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-!") == 0)
			DISC_MASS = 0;
		if (strcmp(argv[i], "-pol") == 0)
			POLAR = 1;
		if (strcmp(argv[i], "-copl") == 0)
			POLAR = 0;
		if (strcmp(argv[i], "-o") == 0)
			sprintf(out_pt, "%s", argv[++i]);
		if (strcmp(argv[i], "-t") == 0)
			TMAX = strtod(argv[++i], NULL);
		if (strcmp(argv[i], "-mr") == 0)
			MASS_RATIO = strtod(argv[++i], NULL);
		if (strcmp(argv[i], "-e") == 0)
			E = strtod(argv[++i], NULL);
		if (strcmp(argv[i], "-n") == 0)
			N_INITIAL = strtol(argv[++i], NULL, 10);
		if (strcmp(argv[i], "-inner") == 0)
			MANUAL_INNER = strtol(argv[++i], NULL, 10);
		if (strcmp(argv[i], "-outer") == 0)
			DISC_OUTER = strtol(argv[++i], NULL, 10);
		if (strcmp(argv[i], "-oi") == 0)
			OUTPUT_INTERVAL = strtol(argv[++i], NULL, 10);
	}
}	
