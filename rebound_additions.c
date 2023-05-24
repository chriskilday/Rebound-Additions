/* rebound_additions.c
 * Some tools useful for running rebound with MPI
 */

// get_com and move_to_com adapted for multiple processors
struct reb_particle get_com(struct reb_simulation* const s) {
	struct reb_particle com;
	if (s->mpi_num == 1)
		com = reb_get_com(s);
	else {
		if (s->mpi_id == 0) {
			struct reb_particle coms[s->mpi_num];
			coms[s->mpi_num-1] = reb_get_com(s);

			MPI_Request request[s->mpi_num-1];
			for (int i = 0; i < s->mpi_num-1; i++) {
				MPI_Irecv(&coms[i], 1, s->mpi_particle, i+1, (i+1)*s->mpi_num, MPI_COMM_WORLD, &(request[i]));
			}

			for (int i = 0; i < s->mpi_num-1; i++) {
				MPI_Status status;
				MPI_Wait(&(request[i]), &status);
			}

			double m_tot = 0;
			for (int i = 0; i < s->mpi_num; i++)
				m_tot += coms[i].m;

			struct reb_particle local_com = {0}; local_com.m = m_tot;
			for (int i = 0; i < s->mpi_num; i++) {
				double mr = coms[i].m/m_tot;
				local_com.x += coms[i].x * mr;
				local_com.y += coms[i].y * mr;
				local_com.z += coms[i].z * mr;
			}

			for (int i = 0; i < s->mpi_num-1; i++)
				MPI_Send(&local_com, 1, s->mpi_particle, i+1, s->mpi_num+(i+1), MPI_COMM_WORLD);

			com = local_com;
		} else {
			struct reb_particle local_com = reb_get_com(s);
			MPI_Send(&local_com, 1, s->mpi_particle, 0, s->mpi_id*s->mpi_num, MPI_COMM_WORLD);

			MPI_Request request;
			MPI_Status status;
			MPI_Irecv(&local_com, 1, s->mpi_particle, 0, s->mpi_num+s->mpi_id, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);

			com = local_com;		
		}
	}

	return com;
}

void move_to_com(struct reb_simulation* const s) {
	if (s->mpi_num == 1)
		reb_move_to_com(s);
	else {
		struct reb_particle com = get_com(s);
		for (int i = 0; i < s->N; i++) {
			s->particles[i].x -= com.x;
			s->particles[i].y -= com.y;
			s->particles[i].z -= com.z;
			s->particles[i].vx -= com.vx;
			s->particles[i].vy -= com.vy;
			s->particles[i].vz -= com.vz;
		}
	}
}


// My version of particle_to_orbit
// Calculates only o.d, o.h, o.v, o.a, o.e, o.n, o.P, o.inc, o.Omega, o.omega
static double acos2(double num, double denom, double disambiguator) {
	double val;
	double cosine = num/denom;
	if (cosine > -1. && cosine < 1.) {
		val = acos(cosine);
		if (disambiguator < 0.)
			val = -val;
	} else val = (cosine <= -1.) ? M_PI : 0;

	return val;
}

struct reb_orbit particle_to_orbit(double G, struct reb_particle pt, struct reb_particle primary) {
	struct reb_orbit o;
	double mu,dx,dy,dz,dvx,dvy,dvz, vsquared, vcircsquared;
	double hx,hy,hz,vr,ex,ey,ez,nx,ny,n;

	dx = pt.x - primary.x;
	dy = pt.y - primary.y;
	dz = pt.z - primary.z;
	dvx = pt.vx - primary.vx;
	dvy = pt.vy - primary.vy;
	dvz = pt.vz - primary.vz;

	o.d = sqrt(dx*dx + dy*dy + dz*dz);
	mu = G * (pt.m + primary.m);

	vsquared = dvx*dvx + dvy*dvy + dvz*dvz;
	vcircsquared = mu / sqrt(dx*dx + dy*dy + dz*dz);

	hx = (dy*dvz - dz*dvy);
	hy = (dz*dvx - dx*dvz);
	hz = (dx*dvy - dy*dvx);
	o.h = sqrt(hx*hx + hy*hy + hz*hz);

	vr = (dx*dvx + dy*dvy + dz*dvz) / sqrt(dx*dx + dy*dy + dz*dz);
	ex = ((vsquared-vcircsquared)*dx - sqrt(dx*dx + dy*dy + dz*dz)*vr*dvx) / mu;
	ey = ((vsquared-vcircsquared)*dy - sqrt(dx*dx + dy*dy + dz*dz)*vr*dvy) / mu;
	ez = ((vsquared-vcircsquared)*dz - sqrt(dx*dx + dy*dy + dz*dz)*vr*dvz) / mu;

	o.v = sqrt(vsquared);
	o.a = -mu / (vsquared - 2*vcircsquared);
	o.e = sqrt(ex*ex + ey*ey + ez*ez);
	o.n = o.a/fabs(o.a)*sqrt(fabs(mu/(o.a*o.a*o.a)));
	o.P = 2*M_PI / o.n;

	o.inc = acos2(hz, o.h, 1);

	nx = -hy; ny = hx;
	n = sqrt(nx*nx + ny*ny);
	o.Omega = acos2(nx,n,ny);

	if (o.inc < 1e-8 || o.inc > M_PI-1e-8) {
		if (o.inc < M_PI/2) o.omega = acos2(ex, o.e, ey) - o.Omega;
		else o.omega = o.Omega - acos2(ex, o.e, ey);
	} else o.omega = acos2(nx*ex + ny*ey, n*o.e, ez);
	o.omega = reb_tools_mod2pi(o.omega);

	return o;
}
