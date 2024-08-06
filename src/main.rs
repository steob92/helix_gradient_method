mod utils;
use utils::{snells_law, line};
mod radiator;
use radiator::{Radiator, Surface};
use gnuplot::*;
use gnuplot::{Caption, Color, Figure, PointSize, PointSymbol};


fn s_curve(theta : f64, t : f64, n : f64) -> f64 {
    let theta_rad = theta.to_radians();
    let n_air = 1.0003;
    let mut inner =  n_air * theta_rad.sin() / n;
    inner = theta_rad.cos()*(inner).asin().tan();
    t * (theta_rad.sin() - inner)
}



fn main() -> Result<(), &'static str>{

    // Run some initial test
    let mut rad = Radiator::new();
    let refractive_index = vec![1.15, 0., 0.];
    let thickness = vec![12.0, 0., 0.];
    let surf_front = vec![10.0, 0., 0.];
    let surf_back = vec![20.0, 0., 0.];
    

    rad.set_offset(0.);
    rad.set_refractive_index(&refractive_index)?;
    rad.set_thickness(&thickness)?;
    rad.set_surface(Surface::Front(&surf_front))?;
    rad.set_surface(Surface::Back(&surf_back))?;
    rad.set_step_size(0.001)?;

    let ang = (45_f64).to_radians();
    rad.theta0 = ang;
    rad.x0 = 0.;
    rad.y0 = -surf_front[0] * ang.tan();

    let x_inter = rad.get_intercept(0., 20., ang, &surf_front);
    let y_inter = line(x_inter, rad.theta0.tan(), rad.y0);
    println!("Intercept: {:0.2}, {:0.2}", x_inter, y_inter);
    let (x_prop, y_prop) = rad.propagate()?;
    println!("Screen intercept: {:0.2}, {:0.2}", x_prop, y_prop);


    // let mut x_prop_vec = vec![];
    let mut y_prop_vec: Vec<f64> = vec![];
    let mut y_no_radiator: Vec<f64> = vec![];
    let mut y_plot: Vec<f64> = vec![];

    let mut angles: Vec<f64> = vec![];
    let n_points: usize = 50;
    let ang_min = -50_f64;
    let ang_max = 50_f64;
    let del_angle = (ang_max - ang_min) / (n_points as f64);
    

    for i in 0..n_points{
        let angle = ang_min + (i as f64) * del_angle;
        angles.push(angle);

        let ang = (angle).to_radians();
        rad.theta0 = ang;
        rad.x0 = 0.;
        rad.y0 = -surf_front[0] * ang.tan();

        y_no_radiator.push(rad.y0 + rad.screen * (ang).tan());

        let (_, y) = rad.propagate()?;
        y_prop_vec.push(y);

        y_plot.push(y_no_radiator[i] - y_prop_vec[i]);
    }
    

    println!("{:#?}", y_prop_vec);

    let s1 = angles.iter().map(|&x| {
        s_curve(x, surf_back[0] - surf_front[0], refractive_index[0])
    }).collect::<Vec<f64>>();

    let s2 = angles.iter().map(|&x| {
        s_curve(x, (surf_back[0] - surf_front[0]) / x.to_radians().cos(), refractive_index[0])
    }).collect::<Vec<f64>>();


    let mut fg = Figure::new();
    fg.axes2d()
        .points(&angles, &y_plot, &[Caption("Ray Tracing"), Color("black"), PointSize(2.), PointSymbol('O')])
        .lines(&angles, &s1, &[Caption("t'=t"), Color("orange")])
        .lines(&angles, &s2, &[Caption("t' = t / cos($\\theta$)"), Color("red"), LineWidth(4.0), LineStyle(DashType::Dot)])
        .set_x_ticks(Some((Auto, 1)), &[], &[])
        .set_y_ticks(Some((Auto, 1)), &[], &[])
        .set_x_label(
            "Displacement [mm]",
            &[],
        )
        .set_y_label("Incodince Angle [degrees]", &[])
        .set_x_grid(true)
        .set_y_grid(true)
        .set_x_minor_grid(true)
        .set_y_minor_grid(true);
        // .set_x_range(Fix(1.), Fix(100.))
        // .set_y_range(Fix(2e-4), Fix(0.2))
        // .set_x_log(Some(10.0))
        // .set_y_log(Some(10.0));
    let _ = fg.save_to_png("test.png", 1100, 600);

    Ok(())
}
