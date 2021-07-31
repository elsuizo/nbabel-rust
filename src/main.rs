#![feature(array_zip)]

use static_math::V3;
use std::fs;
use std::str::FromStr;

const fn generate_steps() -> [i64; 10000] {
    let mut result = [0i64; 10000];
    let mut i = 0;
    while i < 10000 {
        result[i] = i as i64;
        i += 1;
    }
    result
}

const STEPS: [i64; 10000] = generate_steps();

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Particle {
    position: V3<f64>,
    velocity: V3<f64>,
    acceleration: V3<f64>,
    acceleration_0: V3<f64>,
    mass: f64,
}

impl Default for Particle {
    fn default() -> Self {
        Self {
            position: V3::zeros(),
            velocity: V3::zeros(),
            acceleration: V3::zeros(),
            acceleration_0: V3::zeros(),
            mass: 0.0,
        }
    }
}

#[derive(Debug)]
struct Bodies<const N: usize> {
    particles: [Particle; N],
}

impl<const N: usize> Bodies<N> {
    pub fn new(file_path: &str) -> Self {
        let mut particles = [Particle::default(); N];
        for (index, particle) in fs::read_to_string(file_path)
            .expect(&format!("File `{}` not found!", file_path))
            .lines()
            .filter(|x| !x.is_empty())
            .map(|x| parse_row(x))
            .enumerate()
        {
            particles[index] = particle;
        }
        Bodies { particles }
    }

    pub fn update_acceleration(&mut self) {
        // for mut p in self.particles {
        //     p.acceleration = V3::zeros();
        // }
        // for i in 0..N {
        //     for j in 0..N {
        //         if (i != j) {
        //             let rij = self.particles[i].position - self.particles[j].position;
        //             let r_dot_r = rij * rij;
        //             let apre = 1.0 / f64::sqrt(r_dot_r * r_dot_r * r_dot_r);
        //             self.particles[i].acceleration -= self.particles[j].mass * apre * rij;
        //         }
        //     }
        // }
        for (mut p1, p2) in self.particles.zip(self.particles) {
            if p1 != p2 {
                let rij = p1.position - p2.position;
                let r_dot_r = rij * rij;
                let apre = 1.0 / f64::sqrt(r_dot_r * r_dot_r * r_dot_r);
                p1.acceleration -= p2.mass * apre * rij;
            }
        }
    }

    pub fn update_positions(&mut self, dt: f64) {
        self.particles.iter_mut().for_each(|p| {
            // p.acceleration_0 = p.acceleration;
            p.position += dt * p.velocity + 0.5 * dt * dt * p.acceleration_0
        });
    }
    // update velocities based on previous and new accelerations
    pub fn update_velocities(&mut self, dt: f64) {
        self.particles.iter_mut().for_each(|p| {
            p.velocity += 0.5 * dt * (p.acceleration_0 + p.acceleration);
            // p.acceleration_0 = p.acceleration;
        });
    }

    // NOTE(elsuizo:2021-07-29): wide code is better code :)
    pub fn compute_energies(&self) -> (f64, f64) {
        let mut kinetic_energy = 0.0;
        let mut potential_energy = 0.0;
        for i in 0..N {
            let temp_vel = self.particles[i].velocity * self.particles[i].velocity;
            kinetic_energy += 0.5 * self.particles[i].mass * temp_vel;
            for j in i + 1..N {
                // distance between the two stars
                let rij = self.particles[i].position - self.particles[j].position;
                let temp_rij = rij * rij;
                potential_energy -=
                    (self.particles[i].mass * self.particles[j].mass) / f64::sqrt(temp_rij);
            }
        }
        (potential_energy, kinetic_energy)
    }
}

pub fn parse_row(line: &str) -> Particle {
    let row_vec: Vec<f64> = line
        .split_whitespace()
        .skip(1)
        .filter_map(|s| f64::from_str(s).ok())
        .collect();

    Particle {
        position: V3::new_from(row_vec[1], row_vec[2], row_vec[3]),
        velocity: V3::new_from(row_vec[4], row_vec[5], row_vec[6]),
        acceleration: V3::zeros(),
        acceleration_0: V3::zeros(),
        mass: row_vec[0],
    }
}

fn main() {
    let mut bodies: Bodies<128> = Bodies::new("input128");
    // start time, end time and simulation step
    let mut t = 0.0;
    let t_end = 10.0;
    let dt = 0.001;
    let mut k = 0;

    bodies.update_acceleration();
    let (mut potential_energy, mut kinetic_energy) = bodies.compute_energies();
    let total_energy = potential_energy + kinetic_energy;
    for _ in STEPS {
        // update positions based on velocities and accelerations
        bodies.update_positions(dt);

        // get new accelerations
        bodies.update_acceleration();

        // update velocities
        bodies.update_velocities(dt);
    }

    let r = bodies.compute_energies();
    let ekin = r.1;
    let epot = r.0;
    println!(
        "Final dE/E = {}",
        ((ekin + epot) - total_energy) / total_energy
    );
}
