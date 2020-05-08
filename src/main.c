#include "print_particules.h"
#include "particle.h"
#include "SPH.h"
#include "derivatives.h"
#include <math.h>
#include "kernel.h"
#include "consistency.h"

//#include "crtdbg.h" // for memory leak detection; comment if you're on Linux

void script_circle_to_ellipse();
void dam_break();
void box();
void dam_break_flow();

int main() {
	// _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF); // comment if on Linux
	 // script_circle_to_ellipse();
	//script_circle_to_ellipse();
	dam_break();
    //dam_break_flow();
	return EXIT_SUCCESS;
}

void dam_break() {
	// Parameters of the problem
	double R = 0.1;
    
	
    double L = 0.5; // size of the domain: [-L,L] x [-L,L]
    double l = 0.4*L; // particle distribution on [-l,l] x [-l,l]
	double H = 1;
	double dt = 2.e-4; // physical time step
    double T = 5.; // duration of simulation
	bool gravity = 1; // 1 if we consider the gravity

	// Physical parameters
	double rho_0 = 1000.0; // initial (physical) density of water at 20°C (in kg/m^3)
    double rho_1 = 200.;
	double mu_0 = 1.0016e-3; // dynamic viscosity of water at 20°C (in N.s/m^2)
    double mu_1 = 0.086;
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 1.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 72.86e-3; // surface tension of water-air interface at 20°C (in N/m)


	// SPH parameters
	int N_x1 = 41; // flow
	int N_y1 = 41; // flow
    
	int N_x2 = ceil(N_x1*3/2); // pool
	int N_y2 = ceil(N_x2/3); // pool
    
    //initial distance between 2 particles :
    
    
    //double kh = sqrt(21) * 2 * l / 25; // min distance to have 21 particles in kernel
    double kh = 5*l/N_x1;
    printf("kh = %f\n",kh);

	int n_iter = (int)(T / dt); // number of iterations to perform
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 1.5;//1.5; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	
    int T_verlet = 4;
    double L_verlet = 1.1*(2*5.*(double)T_verlet*dt); // 􏱉L = nu*􏱑(2Vmax· C·dt)
    printf("L_verlet = %f\n",L_verlet);
    //double L_verlet = 0.;
    bool use_verlet = true;
    Verlet verlet;
    Verlet_init(&verlet,L_verlet,T_verlet,use_verlet);
    //void* verlet = NULL;
    
    
	double XSPH_epsilon = 0.5;
	Free_surface_detection surface_detection = DIVERGENCE;
    double CR = 0.9;
    double CF = 0.1;

	printf("n_iter = %d\n", n_iter);

	// Animation parameter
	double T_anim = 10; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation


	// Initialize particles on a square

	//int n_p = N_x1 * N_y1 + N_x2 * N_y2;
    int n_p = N_x1 * N_y1;
    N_x2 = 0; N_y2 = 0;
    
    double h_x1 = (0.8*L - 0.005*L) / (N_x1 - 1);
    printf("%f\n", h_x1);
    double h_y1 = (0.8*L - 0.005*L) / (N_y1 - 1);
    printf("%f\n", h_y1);
    double h_x2 = (1.2*L - 0.005*L) / (N_x2 - 1);
    printf("%f\n", h_x2);
    double h_y2 = (0.4*L - 0.005*L) / (N_y2 - 1);
    printf("%f\n", h_y2);
    
	printf("nombre de particules %d \n", n_p);
	    
	double m1 = rho_0 * h_x1 * h_y1; // dam/flow
	double m2 = rho_1 * h_x2 * h_y2; // pool
	Particle** particles = (Particle**)malloc((n_p) * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc((n_p) * sizeof(Particle_derivatives*));
	Residual** residuals = malloc((n_p) * sizeof(Residual*));
    
	//int index1 = N_x2 * N_y2;
    int index1 = 0;
    
	for (int i = 0; i < N_x2; i++) {
		for (int j = 0; j < N_y2; j++) {
			int index = i * (N_y2)+j;
			xy* pos;
			xy* v;
            pos = xy_new(-0.2*L+0.001 + i*h_x2, -L + j * h_y2 + 0.001);
			v = xy_new(0.0, 0.0); // initial velocity = 0
			particles[index] = Particle_new(index, m2, pos, v, rho_1, mu_1, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
		}
	}
	for (int i = 0; i < N_x1; i++) {
		for (int j = 0; j < N_y1; j++) {
			int index = index1 + i*N_y1 + j;
			xy* pos;
			xy* v;
            pos = xy_new(-L+0.001 + i * h_x1,-0.3*L + j * h_y1 + 0.005*L);
			v = xy_new(0.0, 0.0); // initial velocity = 0
			particles[index] = Particle_new(index, m1, pos, v, rho_0, mu_0, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
		}
	}
    int N_flow = N_x1*N_y1;
    //int N_pool = N_x2*N_y2;
    int N_pool = 0;
    
    printf("N_flow = %d, N_pool = %d\n",N_flow, N_pool);
    
	// Setup grid
    Grid* grid = Grid_new_verlet(-1.2*L, 1.2*L,-1.2*L, 1.2*L , kh, &verlet);
    //Grid* grid = Grid_new(-2.*L, 2.*L,-2.*L, 2.*L , kh);
	// Setup BOUNDARY
    Boundary* boundary = Boundary_new(-L, -L, L, -0.3*L, -0.3*L, L, CR, CF);

	// Setup setup
	Setup* setup = Setup_new(n_iter, dt, kh, &verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon, gravity);
	// Setup animation
    Animation* animation = Animation_new(N_pool, N_flow, dt_anim,grid, 1);
    //Animation* animation = NULL;
	// Simulation
	simulate_boundary(grid, particles, particles_derivatives, residuals, n_p, update_positions_seminar_5, setup, animation, boundary);
	// Free memory
	Boundary_free(boundary);
	free_particles(particles, n_p);
	free_particles_derivatives(particles_derivatives, n_p);
	free_Residuals(residuals, n_p);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}

/*
void dam_break_flow(){
	// Parameters of the problem
	double l = 0.057; // particle distribution on [-l,l] x [-l,l]
	double L = 1; // size of the domain: [-L,L] x [-L,L]
	double H = 1;
	double dt = 1.0e-4; // physical time step
	double T = 0.2; // duration of simulation
	bool gravity = 1; // 1 if we consider the gravity

	// Physical parameters
	double rho_0 = 1000.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 1.0016e-3; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 1.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 72.86e-3; // surface tension of water-air interface at 20°C (in N/m)


	// SPH parameters
	double N_x1 = 30;
	double N_y1 = 30;
	double N_x2 = 45;
	double N_y2 = 20;
	double kh = sqrt(21) * 2 * l / 25;

	//int n_per_dim = 60; // number of particles per dimension
	//double kh = sqrt(21) * 2 * l / n_per_dim; // kernel width to ensure 21 particles in the neighborhood
	int n_iter = (int)(T / dt); // number of iterations to perform
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 1.2;//1.5; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	
    
    int T_verlet = 4;
    double L_verlet = 1.1*(2*3.4*(double)T_verlet*dt); // 􏱉L = nu*􏱑(2Vmax· C·dt)
    Verlet verlet;
    Verlet_init(&verlet,L_verlet,T_verlet);
    //void* verlet = NULL;
    
    
	double XSPH_epsilon = 0.5;
	Free_surface_detection surface_detection = DIVERGENCE;
	double CR = 1.0;
	double CF = 0.0;

	printf("n_iter = %d\n", n_iter);

	// Animation parameter
	double T_anim = 10; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation


	// second amount of ball
	//int N_x1 = 25;
	//int N_x2 = 35;
	//double h_x2 = (3 * l) / (N_x2);
	//double h_y2 = l * 0.5 / (N_x2);

	// Initialize particles on a square

	int n_p = N_x1 * N_y1 + N_x2 * N_y2;
	double h_x1 = (2 * l - 0.01) / (N_x1 - 1);
	double h_y1 = (2 * l - 0.01) / (N_y1 - 1);
	double h_x2 = (3 * l - 0.01) / (N_x2 - 1);
	double h_y2 = (0.5 * l ) / (N_y2 - 1);

	//int n_p = squared(n_per_dim);// +N_x * N_y; // total number of particles
	printf("nombre de point %d", n_p);
	//double h_x = (2 * l -0.01) / (n_per_dim - N_x2 - 1); // step between neighboring particles
	//double h_y= 2 * l / (n_per_dim - 1);
	//double m1 = rho_0 * h_x*h_y;
	double m1 = rho_0 * h_x1 * h_y1;
	double m2 = rho_0 * h_x2 * h_y2;
	Particle** particles = (Particle**)malloc((n_p) * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc((n_p) * sizeof(Particle_derivatives*));
	Residual** residuals = malloc((n_p) * sizeof(Residual*));

	int index1 = N_x2 * N_y2;
	for (int i = 0; i < N_y2; i++) {
		for (int j = 0; j < N_x2; j++) {
			int index =  i * (N_x2)+j;
			xy* pos;
			xy* v;
			pos = xy_new((j)*h_x2 + 0.005, -l + i * h_y2);
			v = xy_new(0.0, 0.0); // initial velocity = 0
			particles[index] = Particle_new(index, m2, pos, v, rho_0, mu, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
		}
	}
	for (int i = 0; i < N_y1; i++) {
		for (int j = 0; j < N_x1; j++) {
			int index = index1+i * (N_x1)+j;
			xy* pos;
			xy* v;
			pos = xy_new(-4 * l + i * h_x1, j * h_y1);
			v = xy_new(0.0, 0.0); // initial velocity = 0
			particles[index] = Particle_new(index, m1, pos, v, rho_0, mu, c_0, gamma, sigma);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
		}
	}
	// Setup grid
	Grid *grid = Grid_new_verlet(-5*l, 4*l, -H, H, kh, &verlet);
	// Setup BOUNDARY
	double lb = 0.420;
	double hb = 0.440;
	double Rp = 0.001; //particle radius
	Boundary* boundary = Boundary_new(-2*l-Rp,-Rp,3*l+Rp,-Rp,-l-Rp,hb-l+Rp,CR,CF);
    

	int n_beg = N_x2 * N_y2;

	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, &verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon, gravity);
	// Setup animation
	Animation *animation = Animation_new(n_p, dt_anim, grid, 1);
	// Simulation
	simulate_boundary_flow(grid, particles, particles_derivatives, residuals, n_beg, update_positions_seminar_5, setup, animation, boundary);
	// Free memory
	Boundary_free(boundary);
	free_particles(particles, n_p);
	free_particles_derivatives(particles_derivatives, n_p);
	free_Residuals(residuals, n_p);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}
// Evolution of a 2D circle with non-zero initial velocities (no surface tension force)
// Test case from "Simulating Free Surface Flows with SPH", Monaghan (1994)

*/
