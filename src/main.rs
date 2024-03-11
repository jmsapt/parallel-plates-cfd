// TODO: write a local macro within propagate to make it more readable
mod sim;
use anyhow::Result;

const OUTPUT_FILE: &str = "plot.png";
use sim::{ChartType, Fluid, Paramaters};

fn main() -> Result<()> {
    let mut example = Fluid::new(Paramaters::Water(2.0), 0.2,200, 2e-3, 100);
    let n_iterations = 20_000;
    for _ in 0..n_iterations {
        example.propagate(0.005);
    }

    example.draw_chart(ChartType::Velocity, "plates.png")?;

    let max = example
        .grid
        .last()
        .unwrap()
        .velocity
        .iter()
        .fold((0_f64, 0_f64), |max, x| {
            if (x.0.powi(2) + x.1.powi(2)).sqrt() > (max.0.powi(2) + max.1.powi(2)).sqrt() {
                *x
            } else {
                max
            }
        });

    println!("Max velocity: {max:?}");
    Ok(())
}
