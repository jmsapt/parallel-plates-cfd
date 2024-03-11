use std::{fmt, ops::Index};

use anyhow::Result;
use plotters::{prelude::*, style::HSLColor};
use rayon::prelude::*;

// Poisson's equation convergence parameters
const EPSILON: f64 = 1e-9;
const MAX_ITERATIONS: u32 = 100;

pub enum ChartType {
    Pressure,
    Velocity,
}

pub struct Paramaters {
    pub rho: f64,
    pub gravity: f64,
    pub mu: f64,
    pub velocity: f64,
}
impl Paramaters {
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

#[derive(Debug)]
pub struct Grid {
    pub time: f64,
    pub velocity: Vec<(f64, f64)>,
    pub pressure: Vec<f64>,
}

pub struct Fluid {
    params: Paramaters,
    length: f64,
    x_divisions: usize,
    width: f64,
    y_divisions: usize,

    pub grid: Vec<Grid>,
}
impl Fluid {
    pub fn new(
        fluid_params: Paramaters,
        length: f64,
        x_divisions: usize,
        width: f64,
        y_divisions: usize,
    ) -> Self {
        let num_cells = x_divisions * y_divisions;

        let mut velocity = vec![(2.0, 0.0); num_cells];
        for j in 0..y_divisions {
            velocity[j * x_divisions] = (fluid_params.velocity, 0.0);
        }
        Self {
            params: fluid_params,
            length,
            x_divisions,
            width,
            y_divisions,
            grid: vec![Grid {
                time: 0.0,
                pressure: vec![0.0; num_cells],
                velocity,
            }],
        }
    }

    pub fn propagate(&mut self, delta_t: f64) {
        // let num_cells = self.x_divisions * self.y_divisions;
        let prev_grid = self.grid.last().expect("Grid vec is empty");

        // previous (i) values
        let pressure = &prev_grid.pressure;
        let velocity = &prev_grid.velocity;

        // current (i + 1) values
        // let mut new_pressure: Vec<f64> = Vec::with_capacity(num_cells);
        // let mut new_velocity: Vec<(f64, f64)> = Vec::with_capacity(num_cells);
        let mut new_pressure = pressure.clone();
        let mut new_velocity = velocity.clone();

        // 0. Calculate unifroms
        let time = prev_grid.time + delta_t;
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
                        let p_init = *p;
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
                                        - velocity[self.index(i - 1, j)].1)));

                        if (p_init - *p).abs() / p_init < EPSILON {
                            break;
                        }
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

        // 3. Push new velocities and pressures
        self.grid.push(Grid {
            time,
            velocity: new_velocity,
            pressure: new_pressure,
        })
    }

    pub fn draw_chart(&self, chart_type: ChartType, path: &str) -> Result<()> {
        let root = BitMapBackend::new(path, (1000, 500)).into_drawing_area();
        root.fill(&WHITE)?;

        // create chart
        let mut chart = ChartBuilder::on(&root)
            .margin(40)
            .x_label_area_size(20)
            .y_label_area_size(20)
            .build_cartesian_2d(0.0..self.length, (-self.width / 2.0)..(self.width / 2.0))?;

        // configure chart
        chart
            .configure_mesh()
            .disable_x_mesh()
            .disable_y_mesh()
            .draw()?;

        let plotting_area = chart.plotting_area();

        match chart_type {
            ChartType::Pressure => todo!(),
            ChartType::Velocity => {
                let min = 0.0;
                let max = 3.0;

                let grid = self.grid.last().unwrap();

                let (x_range, y_range) = plotting_area.get_pixel_range();

                for x in x_range.start..x_range.end {
                    for y in y_range.start..y_range.end {
                        let i = (self.x_divisions as i32 * (x - x_range.start))
                            / (x_range.end - x_range.start);
                        let j = (self.y_divisions as i32 * (y - y_range.start))
                            / (y_range.end - y_range.start);
                        let index = self.index(i as usize, j as usize);

                        let mag = (grid.velocity[index].0.powi(2) + grid.velocity[index].1.powi(2))
                            .sqrt();
                        root.draw_pixel((x, y), &HSLColor((mag - min) / (max - min), 1.0, 0.5))?;
                        if mag > 0.5 {
                            // if mag < 10.0 {
                            //     println!("mag {mag:.6}")
                            // }
                        }
                    }
                }
            }
        }
        root.present()?;
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

fn map_color(mag: f64, min: f64, max: f64) -> HSLColor {
    let normalized_velocity = (mag - min) / (max - min);
    // dbg!(normalized_velocity, mag);
    // Scale the velocity to the hue range (blue to red)
    let hue = 240.0 - normalized_velocity * 240.0;
    HSLColor(hue, 1.0, 0.5)
}

#[cfg(test)]
mod Test {
    use super::*;
    #[test]
    fn test_map_color() {
        assert_eq!(map_color(10.0, 0.0, 10.0), HSLColor(0.0, 1.0, 1.0));
        assert_eq!(map_color(0.0, 0.0, 10.0), HSLColor(240.0, 1.0, 1.0));
    }
}
