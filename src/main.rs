// TODO: write a local macro within propagate to make it more readable

use core::num;
use std::{ops::Index, thread::sleep};

use anyhow::Result;
use plotters::{prelude::*, style::HSLColor};
use rayon::prelude::*;

const OUTPUT_FILE: &str = "plot.png";

// Poisson's equation convergence parameters
const EPSILON: f64 = 1e-3;
const MAX_ITERATIONS: u32 = 100;

enum ChartType {
    Pressure,
    Velocity,
}

struct Paramters {
    pub rho: f64,
    pub gravity: f64,
    pub mu: f64,
    pub velocity: f64,
}
impl Paramters {
    /// Fluid properites for water at 25 degrees
    pub fn Water(velocity: f64) -> Self {
        Self {
            rho: 1000.0,
            gravity: -9.8,
            mu: 8.90e-4,
            velocity,
        }
    }
}

fn main() -> Result<()> {
    let example = Fluid::new(Paramters::Water(10.0), 4.0, 400, 2.0, 200);
    let n_interations = 1000;
    for _ in 0..n_interations {}

    Ok(())
}

struct Grid {
    time: f64,
    velocity: Vec<(f64, f64)>,
    pressure: Vec<f64>,
}

struct Fluid {
    params: Paramters,
    length: f64,
    x_divisions: usize,
    width: f64,
    y_divisions: usize,

    grid: Vec<Grid>,
}
impl Fluid {
    pub fn new(
        fluid_params: Paramters,
        length: f64,
        x_divisions: usize,
        width: f64,
        y_divisions: usize,
    ) -> Self {
        let num_cells = x_divisions * y_divisions;
        Self {
            params: fluid_params,
            length,
            x_divisions,
            width,
            y_divisions,
            grid: vec![Grid {
                time: 0.0,
                pressure: vec![0.0; num_cells],
                velocity: vec![(0.0, 0.0); num_cells],
            }],
        }
    }

    pub fn propagate(&mut self, delta_t: f64) {
        let num_cells = self.x_divisions * self.y_divisions;
        let prev_grid = self.grid.last().expect("Grid vec is empty");

        // previous (i) values
        let pressure = &prev_grid.pressure;
        let velocity = &prev_grid.velocity;

        // current (i + 1) values
        let mut new_pressure: Vec<f64> = Vec::with_capacity(num_cells);
        let mut new_velocity: Vec<(f64, f64)> = Vec::with_capacity(num_cells);

        // 0. Calculate unifroms
        let delta_x = self.length / self.x_divisions as f64;
        let delta_y = self.width / self.y_divisions as f64;
        let rho = self.params.rho;
        let v_inlet = self.params.velocity;

        // 1. Solve Poisson's euqation at every cell
        new_pressure
            .par_iter_mut()
            .enumerate()
            .for_each(|(index, p)| {
                // get x and y coords from buffer index
                let (i, j) = (index % self.x_divisions, index / self.x_divisions);

                // inlet boundary
                if i == 0 {
                    *p = pressure[self.index(1, j)] - (rho * v_inlet * delta_x) / 2.0;
                }
                // outlet boundary
                else if i == self.x_divisions - 1 {
                    *p = 0.0;
                }
                // wall condition (copy press from adjacent cell to force 0 pressure gradient)
                // top wall
                else if j == 0 {
                    *p = pressure[self.index(i, j + 1)]
                }
                // bottom wall
                else if j == self.y_divisions - 1 {
                    *p = pressure[self.index(i, j - 1)]
                }
                // normal case
                else {
                    for _ in 0..MAX_ITERATIONS {
                        // see [method](method.md) for equation
                        *p = 0.25
                            * (pressure[self.index(i, j - 1)]
                                + pressure[self.index(i, j + 1)]
                                + pressure[self.index(i - 1, j)]
                                + pressure[self.index(i + 1, j)]
                                - (rho * delta_x * delta_y / delta_t
                                    * (2.0 * delta_x * velocity[self.index(i, j + 1)].0
                                        - velocity[self.index(i, j - 1)].0
                                        + 2.0 * delta_y * velocity[self.index(i + 1, j)].1
                                        - velocity[self.index(i - 1, j)].1)))
                    }
                }
            });

        // 2. Apply momentum euqation for x, then y directions
        new_velocity
            .par_iter_mut()
            .enumerate()
            .for_each(|(index, v)| {
                // get x and y coords from buffer index
                let (i, j) = (index % self.x_divisions, index / self.x_divisions);

                // inlet boundary
                if i == 0 {
                    *v = (v_inlet, 0.0);
                }
                // outlet (extrapolate interior velocities)
                else if i == self.x_divisions - 1 {
                    *v = velocity[self.index(i - 1, j)]
                }
                // no-slip wall condition
                else if j == 0 || j == self.y_divisions - 1 {
                    *v = (0.0, 0.0);
                } else {
                    let nu = self.params.mu / rho;
                    // u velocity (in horizontal x-direction)
                    v.0 = velocity[self.index(i, j)].0
                        - delta_t / rho
                            * (pressure[self.index(i + 1, j)] - pressure[self.index(i, j)])
                        + delta_t
                            * nu
                            * ((velocity[self.index(i + 1, j)].0
                                - 2.0 * velocity[self.index(i, j)].0
                                + velocity[self.index(i - 1, j)].0)
                                / delta_x.powi(2)
                                + (velocity[self.index(i, j + 1)].0
                                    - 2.0 * velocity[self.index(i, j)].0
                                    + velocity[self.index(i, j - 1)].0)
                                    / delta_x.powi(2));

                    // v velocity (in vertical y-direction)
                    v.1 = velocity[self.index(i, j)].1
                        - delta_t / rho
                            * (pressure[self.index(i + 1, j)] - pressure[self.index(i, j)])
                        + delta_t
                            * nu
                            * ((velocity[self.index(i + 1, j)].1
                                - 2.0 * velocity[self.index(i, j)].1
                                + velocity[self.index(i - 1, j)].1)
                                / delta_x.powi(2)
                                + (velocity[self.index(i, j + 1)].1
                                    - 2.0 * velocity[self.index(i, j)].1
                                    + velocity[self.index(i, j - 1)].1)
                                    / delta_x.powi(2));
                }
            });
    }

    pub fn draw_chart(&self, chart_type: ChartType, path: &str) -> Result<()> {
        let root = BitMapBackend::new(path, (800, 400)).into_drawing_area();
        root.fill(&WHITE)?;

        // create chart
        let mut chart = ChartBuilder::on(&root)
            .margin(20)
            .x_label_area_size(10)
            .y_label_area_size(10)
            .build_cartesian_2d(0.0..self.length, (-self.width / 2.0)..(self.width / 2.0))?;

        // configure chart
        let ctx = chart
            .configure_mesh()
            .disable_x_mesh()
            .disable_y_mesh()
            .draw()?;

        let plotting_area = chart.plotting_area();

        match chart_type {
            ChartType::Pressure => todo!(),
            ChartType::Velocity => {
                let min_v = 0.0;
                let max_v = 1.0;

                let color = |u: (f64, f64)| {};

                self.grid
                    .last()
                    .unwrap()
                    .velocity
                    .iter()
                    .enumerate()
                    .for_each(|(index, v)| {
                        let (i, j) = (index % self.x_divisions, index / self.x_divisions);
                        let mag = (v.0.powi(2) + v.1.powi(2)).sqrt();

                        plotting_area
                            .draw_pixel((i as f64, j as f64), &map_color(mag, 0.0, 10.0))
                            .expect("Failed to write pixel");
                    })
            }
        }

        Ok(())
    }

    pub fn prograte_n_steps(&mut self, n: usize, timestep: f64) {
        for _ in 0..n {
            self.propagate(timestep)
        }
    }

    /// Calculate padded buffer index given the cartesian coords.
    fn index(&self, i: usize, j: usize) -> usize {
        j * self.x_divisions + i
    }
}

fn map_color(mag: f64, max: f64, min: f64) -> HSLColor {
    let normalized_velocity = (mag - min) / (max - min);
    // Scale the velocity to the hue range (blue to red)
    let hue = 240.0 - normalized_velocity * 240.0;
    HSLColor(hue, 1.0, 1.0)
}
