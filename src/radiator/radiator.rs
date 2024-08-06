#[macro_use]
use peroxide::fuga::*;
use anyhow::Result;
use peroxide::newton;

use std::io::Error;
use std::f64::consts::{FRAC_PI_2, PI};
use crate::utils::{dist, get_angle, get_norm, line, parabola, snells_law, N_AIR};

pub enum Surface <'a>{
    Front(&'a Vec<f64>),
    Back(&'a Vec<f64>),  
}

#[derive(Debug)]
pub struct Radiator{
    // # Initial Conditions
    // ## Start location
    pub x0 : f64,
    pub y0 : f64,
    pub y_shift :  f64,
    pub x_shift :  f64,
    // ## incodence angle
    pub theta0 : f64,
    // ## Step size
    dx : f64,
    // ## Imaging plane location
    pub screen : f64,
    // # Define a 12 mm region
    radiator_world_xmin : f64,
    radiator_world_xmax : f64,
    // # Evolution of refractive index in y
    n_y : Vec<f64>,
    // # Front aerogel surface
    surf_front : Vec<f64>,
    // # Back areogel surface
    surf_back : Vec<f64>,
    // # Thickness of the aerogel
    thickness : Vec<f64>,
    // # Radiator offset
    rad_offset_x : f64,
    delta_min : f64,
}


impl Radiator{

    pub fn new() -> Self{
        Radiator{
            x0 : 0.,
            y0 : 0.,
            y_shift :  0.,
            x_shift :  0.,
            // ## incodence angle
            theta0 : 0.,
            // ## Step size
            dx : 0.001,
            // ## Imaging plane location
            screen : 30.,
            // # Define a 12 mm region
            radiator_world_xmin : 7.,
            radiator_world_xmax : 30.,
            // # Evolution of refractive index in y
            // n_y : vec![1.15855784e+00, -5.83882672e-05,  1.09534706e-02],
            n_y : vec![1.15, 0.0, 0.0],
            // # Front aerogel surface
            surf_front : vec![10., 0., 0.],
            // # Back areogel surface
            surf_back : vec![20., 0., 0.],
            // # Thickness of the aerogel
            thickness : vec![12., 0., 0.,],
            // # Radiator offset
            rad_offset_x : 0.0,
            delta_min : 1e-3,
        }
    }


    pub fn set_step_size(&mut self, dx : f64) -> Result<(), &'static str> {
        if dx <= 0.{
            Err("Step size must be a positive integer")
        } else{
            self.dx = dx;
            Ok(())   
        }
    }

    pub fn set_refractive_index(&mut self, n_y: &Vec<f64>) -> Result<(), &'static str>{
        if n_y.len() != 3 {
            Err("Insufficient parameters for refractive index")
        } else{
            for i in 0..self.n_y.len(){
                self.n_y[i] = n_y[i];
            }
            Ok(())
        }
    }

    pub fn set_surface(&mut self, surf : Surface) -> Result<(), &'static str> {

        match surf {
            Surface::Front(parms) =>{
                if parms.len() != 3{
                    Err("Insufficient Parameters for front surface")
                } else{
                    for i in 0..self.surf_front.len(){
                        self.surf_front[i] = parms[i];
                    }
                    Ok(())
                }
            },
            Surface::Back(parms) =>{
                if parms.len() != 3{
                    Err("Insufficient parameters for back surface")
                } else{
                    for i in 0..self.surf_front.len(){
                        self.surf_back[i] = parms[i];
                    }
                    Ok(())
                }
            }
        }
    }


    pub fn set_thickness(&mut self, parms: &Vec<f64>) -> Result<(), &'static str> {
        
        if parms.len() != 3{
            Err("Insufficient parameters for thickness")
        } else{
            for i in 0..self.thickness.len(){
                self.thickness[i] = parms[i];
            }
            Ok(())
        }
    }


    pub fn set_offset(&mut self, offset: f64) {
        self.rad_offset_x = offset;
    }


    
    pub fn get_intercept(&self, xmin : f64, xmax : f64, theta : f64, surf : &Vec<f64>) -> f64 {
        // Doing something stupid until I find a better solution
        // Brute force line search
        let n_points = ((xmax - xmin) / self.delta_min) as usize;
        let mut dist_min = 9999.;
        let mut x_min = 0.;

        (0..n_points).for_each(|i|{
            let x = i as f64 * self.delta_min + xmin;
            let d = dist(x, theta, self.y0,  surf, self.rad_offset_x);
            // println!("{}, {}", x, d);
            if d < dist_min {
                // println!("New min: {}, {}, {}", x, d, dist_min);
                dist_min = d;
                x_min = x;
            }
        });
        x_min
    }


    pub fn get_refractive_index(&self, y : f64) -> f64{
        parabola(y, self.y_shift, &self.n_y)
    }


    pub fn is_in_radiator(&self, x : f64, y : f64) -> bool{
        let x_lower = parabola(y, self.y_shift, &self.surf_front);
        let x_upper = parabola(y, self.y_shift, &self.surf_back);
        x > x_lower && x < x_upper
        
    }

    pub fn propagate(&self) -> Result<(f64, f64), &'static str> {

        // # determine the number of points in the "world"
        let npoints = ((self.radiator_world_xmax - self.radiator_world_xmin) / self.dx) as usize;

        // # Store the y, n and angles
        // y_points = np.zeros(npoints)
        // x_points = np.arange(npoints) * self.dx + self.radiator_world_xmin
        // theta_points = np.zeros(npoints)
        // n_path = np.zeros(npoints)


        // # Get the intersection point between the ray and the material
        let theta = self.theta0;
        let x_inter = self.get_intercept(self.rad_offset_x, self.rad_offset_x + self.surf_back[0], theta, &self.surf_front);
        let y_inter = line(x_inter,self.theta0.tan(), self.y0);

        let norm = get_norm(y_inter,
            self.surf_front[0], 
            self.surf_front[1],
            self.surf_front[1]
        );
        let norm_inv = norm.iter().map(|x| -x).collect::<Vec<_>>();

        let angle_norm = FRAC_PI_2 - get_angle(&vec![x_inter - self.x0, y_inter - self.y0], &norm_inv);
        let theta_surf1 = get_angle(&vec![1.,0.], &norm);


        // println!("1st Interface-> Intersection Point: [{x_inter:0.2},{y_inter:0.2}] mm");
        // println!("1st Interface-> Angle between ray and surface: {:0.2} degrees", angle_norm.to_degrees());
        // println!("1st Interface-> Angle between surface and coordinate system: {:0.2} degrees", theta_surf1.to_degrees());
        // println!("1st Interface-> Angle between normal to surface and coordinate system: {:0.2} degrees", (theta_surf1 - FRAC_PI_2).to_degrees());

        let mut n_current = self.get_refractive_index(y_inter);
        let theta_prop = snells_law(N_AIR, n_current, angle_norm).unwrap();

        // println!("1st Interface-> Snells Law: {:0.2} degrees", theta_prop.to_degrees());

        let mut angle_prop = theta_surf1 - FRAC_PI_2 + theta_prop;
        // println!("1st Interface-> Angle of propagation: {:0.2} degrees", angle_prop.to_degrees());
        
        let mut x = x_inter;
        let mut y = y_inter;
        'in_radiator : loop {
            // take a step
            y = y + self.dx * angle_prop.tan();
            x += self.dx;

            // Get refractive index of the next step
            let n_next = self.get_refractive_index(y);

            // Snell's Law for the next step
            angle_prop = snells_law(n_current, n_next, angle_prop).unwrap();

            // Update refractive index
            n_current = n_next;

            if !self.is_in_radiator(x, y){
                break 'in_radiator;
            }
        }


        // # Get the normal vector at the point of intersection
        let norm_exit = get_norm(y,
            self.surf_back[0], 
            self.surf_back[1],
            self.surf_back[1]
        );
        let norm_exit_inv = norm_exit.iter().map(|x| -x).collect::<Vec<_>>();


        // # Angle between surface and coordinates
        let theta_surf2 = get_angle(&vec![1.,0.], &norm_exit); 
        let angle_norm_surf2 = theta_surf2 - FRAC_PI_2 - angle_prop;

        // println!("2nd Interface-> Angle between ray and surface 2: {:0.2} degrees", angle_norm_surf2.to_degrees());
        // println!("2nd Interface-> Angle between surface 2 and coordinate system: {:0.2} degrees", theta_surf2.to_degrees());

        let theta_exit = snells_law(
            n_current,
            N_AIR, 
            angle_norm_surf2
        ).unwrap();
        // println!("2nd Interface-> Snell's law: {:0.3}", theta_exit.to_degrees());
        let theta_exit_prop = FRAC_PI_2 - theta_surf1 - theta_exit;
        // println!("2nd Interface-> Propagation Angle: {:0.3} degrees", theta_exit_prop.to_degrees());

        let x_exit = parabola(
            y, 
            self.y_shift, 
            &self.surf_back
        ) + self.rad_offset_x;

        let x_prop = self.screen;
        let y_prop = y + (x_prop - x_exit) * theta_exit_prop.tan();

        // println!("Distance to travel {:0.5}", x_prop - x_exit);
        // println!("Intersection on screen at: ({:0.3}, {:0.3})", x_prop, y_prop );


        Ok((x_prop, y_prop))

    }





}