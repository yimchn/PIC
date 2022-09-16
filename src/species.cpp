#include "species.h"

#include "utlis.h"

/*updates velocities and positions of all particles of this species*/
void Species::advance() {
    /*get the time step*/
    double dt = domain.getDt();

    /*save mesh bounds*/
    Vec2d x0 = domain.geo.x0;
    Vec2d xm = domain.geo.xm;

    /*continue while particles remain*/
    for (Particle &part : particles) {
        /*get logical coordinate of particle's position*/
        Vec2d lc = domain.XtoL(part.pos);

        /*electric field at particle position*/
        Vec3d ef_part = domain.E.gather(lc);

        /*update velocity from F=qE*/
        part.vel[0] += ef_part[0] * (dt * charge / mass);
        part.vel[1] += ef_part[1] * (dt * charge / mass);

        /*update position from v=dx/dt*/
        part.pos += part.vel * dt;

        /*did this particle leave the domain? reflect back*/
        for (int i = 0; i < 2; i++) {
            if (part.pos[i] < x0[i]) {
                part.pos[i] = 2 * x0[i] - part.pos[i];
                part.vel[i] *= -1.0;
            } else if (part.pos[i] >= xm[i]) {
                part.pos[i] = 2 * xm[i] - part.pos[i];
                part.vel[i] *= -1.0;
            }
        }
    }
}

/*compute number density*/
void Species::computeNumberDensity() {
    den.clear();
    for (Particle &part : particles) {
        Vec2d lc = domain.XtoL(part.pos);
        den.scatter(lc, part.mpw);
    }

    // divide by node volume
    den /= domain.geo.node_area;
}

/*adds a new particle, rewinding velocity by half dt*/
void Species::addParticle(Vec2d pos, Vec2d vel, double mpw) {
    // don't do anything (return) if pos outside domain bounds [x0,xd)
    if (!domain.geo.InBounds(pos)) return;

    // get particle logical coordinate
    Vec2d lc = domain.XtoL(pos);

    // evaluate electric field at particle position
    Vec3d ef_part = domain.E.gather(lc);

    // rewind velocity back by 0.5*dt*ef
    vel[0] -= charge / mass * ef_part[0] * (0.5 * domain.getDt());
    vel[1] -= charge / mass * ef_part[1] * (0.5 * domain.getDt());

    // add to list
    particles.emplace_back(pos, vel, mpw);
}

/*loads randomly distributed particles in a x1-x2 box representing num_den
 * number density*/
void Species::loadParticlesSquare(Vec2d x1, Vec2d x2, double num_den,
                                  int num_mp) {
    double square_area = (x2[0] - x1[0]) * (x2[1] - x1[1]);  // box volume
    double num_real = num_den * square_area;  // number of real particles
    double mpw = num_real / num_mp;           // macroparticle weight

    /*preallocate memory for particles*/
    particles.reserve(num_mp);

    /*load particles on an equally spaced grid*/
    for (int p = 0; p < num_mp; p++) {
        // sample random position
        Vec2d pos;
        pos[0] = x1[0] + Rnd()() * (x2[0] - x1[0]);
        pos[1] = x1[1] + Rnd()() * (x2[1] - x1[1]);

        // set initial velocity
        Vec2d vel{0, 0};  // stationary particle

        addParticle(pos, vel, mpw);  // add a new particle to the array
    }
}

/*quiet start load of num_sim[0]*num_sim[1]*num_sim[2] particles in a x1-x2 box
representing num_den number density*/
void Species::loadParticlesBoxQS(Vec2d x1, Vec2d x2, double num_den,
                                 Vec2i num_mp) {
    double square_area = (x2[0] - x1[0]) * (x2[1] - x1[1]);  // box volume
    int num_mp_tot = (num_mp[0] - 1) *
                     (num_mp[1] - 1);  // total number of simulation particles
    double num_real = num_den * square_area;  // number of real particles
    double mpw = num_real / num_mp_tot;       // macroparticle weight

    /*compute particle grid spacing*/
    double di = (x2[0] - x1[0]) / (num_mp[0] - 1);
    double dj = (x2[1] - x1[1]) / (num_mp[1] - 1);

    /*preallocate memory for particles*/
    particles.reserve(num_mp_tot);

    /*load particles on a equally spaced grid*/
    for (int i = 0; i < num_mp[0]; i++)
        for (int j = 0; j < num_mp[1]; j++) {
            double pos[2];
            pos[0] = x1[0] + i * di;
            pos[1] = x1[1] + j * dj;

            // shift particles on max faces back to the domain
            if (pos[0] == x2[0]) pos[0] -= 1e-4 * di;
            if (pos[1] == x2[1]) pos[1] -= 1e-4 * dj;

            double w = 1;  // relative weight
            if (i == 0 || i == num_mp[0] - 1) w *= 0.5;
            if (j == 0 || j == num_mp[1] - 1) w *= 0.5;

            /*add rewind*/
            double vel[2] = {0, 0};  // particle is stationary

            addParticle(pos, vel,
                        mpw * w);  // add a new particle to the array
        }
}

Species::Species(std::string name, double mass, double charge, Domain &domain)
    : name(name),
      mass(mass),
      charge(charge),
      den(domain.geo.ni, domain.geo.nj),
      domain(domain) {}

size_t Species::getNp() { return particles.size(); }

/*returns the number of real particles*/
double Species::getRealCount() {
    double mpw_sum = 0;
    for (Particle &part : particles) {
        mpw_sum += part.mpw;
    }

    return mpw_sum;
}

/* returns the species momentum*/
Vec2d Species::getMomentum() {
    Vec2d mom;
    for (Particle &part : particles) {
        mom += part.mpw * part.vel;
    }

    return mass * mom;
}

/* returns the species kinetic energy*/
double Species::getKE() {
    double ke = 0;
    for (Particle &part : particles) {
        double v2 = part.vel[0] * part.vel[0] + part.vel[1] * part.vel[1];
        ke += part.mpw * v2;
    }
    return 0.5 * mass * ke;
}